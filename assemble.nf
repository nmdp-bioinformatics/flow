#!/usr/bin/env nextflow

/*

    flow  Consensus assembly, variant calling, and allele interpretation workflow.
    Copyright (c) 2014 National Marrow Donor Program (NMDP)

    This library is free software; you can redistribute it and/or modify it
    under the terms of the GNU Lesser General Public License as published
    by the Free Software Foundation; either version 3 of the License, or (at
    your option) any later version.

    This library is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; with out even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
    License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this library;  if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA.

    > http://www.gnu.org/licenses/lgpl.html

*/

params.experiment = "tutorial"
params.reference = "tutorial/ref/chr6-ex.fa"

raw = "${params.experiment}/raw"
intermediate = "${params.experiment}/intermediate"
fine = "${params.experiment}/final"

raw1 = "${raw}/**_R1*{fastq,fq,fastq.gz,fq.gz}"
raw2 = "${raw}/**_R2*{fastq,fq,fastq.gz,fq.gz}"
reads1 = Channel.fromPath(raw1).ifEmpty { error "cannot find any reads matching ${raw1}" }.map { path -> tuple(sample(path), path) }
reads2 = Channel.fromPath(raw2).ifEmpty { error "cannot find any reads matching ${raw2}" }.map { path -> tuple(sample(path), path) }

alignmentReadPairs = Channel.create()
readPairs = reads1.phase(reads2).map{ read1, read2 -> [ read1[0], read1[1], read2[1] ] }.tap(alignmentReadPairs)

process fastqToSsake {
  input:
    set s, file(r1), file(r2) from readPairs
  output:
    set s, file("${s}.ssake.fa.gz") into ssakeFasta

  """
  ngs-fastq-to-ssake -1 ${r1} -2 ${r2} -o ${s}.ssake.fa.gz --insert-size 500
  """
}

process reformat {
  input:
    set s, file(f) from ssakeFasta
  output:
    set s, file("${s}.ssake.reformatted.fa.gz") into reformatted

  """
  zcat $f | sed -e 's/a/A/g' -e 's/c/C/g' -e 's/g/G/g' -e 's/t/T/g' -e 's/n/N/g' -e 's/^>.*/>reformatted:500/' | gzip -c > ${s}.ssake.reformatted.fa.gz
  """
}

process ssake {
  input:
    set s, file(f) from reformatted
  output:
    set s, file("${s}.ssake.d/${s}.contigs") into contigs

  """
  gunzip -c ${f} > ${f}.tmp
  mkdir ${s}.ssake.d
  SSAKE -f ${f}.tmp -b ${s}.ssake.d/${s} -w 1 -h 1 -p 1 -m 50 -o 30 -c 1 -e 0.90 -k 4 -a 0.1 -x 20
  rm ${f}.tmp
  """
}

ref = file("${params.reference}")

process alignContigs {
  input:
    set s, file(f) from contigs
  output:
    set s, file("${s}.contigs.bwa.sorted.bam") into contigsBam
    set s, file("${s}.contigs.bwa.sorted.vcf.gz") into contigsVcf

  """
  bwa bwasw -b 1 ${ref} ${f} | samtools view -hub - | samtools sort -l 0 -T ${s}.contigs.bwa.tmp -o ${s}.contigs.bwa.sorted.bam
  samtools index ${s}.contigs.bwa.sorted.bam
  samtools mpileup -RB -C 0 -Q 0 -f ${ref} ${s}.contigs.bwa.sorted.bam -v -u | gzip -c > ${s}.contigs.bwa.sorted.vcf.gz
  """
}

contigsBam.subscribe() {
  println "aligned contigs $it"
}

contigsVcf.subscribe() {
  println "contigs vcf $it"
}

process interleave {
  input:
    set s, file(r1), file(r2) from alignmentReadPairs
  output:
    set s, file("${s}.paired.fq.gz") into interleavedReads

  """
  ngs-interleave-fastq -1 ${r1} -2 ${r2} -p ${s}.paired.fq.gz -u ${s}.unpaired.fq.gz
  """
}

process alignReads {
  input:
    set s, file(r) from interleavedReads
  output:
    set s, file("${s}.reads.bwa.sorted.bam") into readsBam
    set s, file("${s}.reads.bwa.sorted.vcf.gz") into readsVcf

  """
  bwa mem ${ref} -a -p ${r} | samtools view -hub - | samtools sort -l 0 -T ${s}.reads.bwa.tmp -o ${s}.reads.bwa.sorted.bam
  samtools index ${s}.reads.bwa.sorted.bam
  samtools mpileup -f ${ref} ${s}.reads.bwa.sorted.bam -v -u | gzip -c > ${s}.reads.bwa.sorted.vcf.gz
  """
}

readsBam.subscribe() {
  println "aligned reads $it"
}

readsVcf.subscribe() {
  println "reads vcf $it"
}

def sample(Path path) {
  def name = path.getFileName().toString()
  int start = Math.max(0, name.lastIndexOf('/'))
  return name.substring(start, name.indexOf("_R"))
}