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

params.sample = "example"
params.r1 = "${params.sample}_R1.fq.gz"
params.r2 = "${params.sample}_R2.fq.gz"

sample = Channel.from(params.sample)
fastq = Channel.from(file(params.r1), file(params.r2))

process fastqToSsake {    
    input:
        val s from sample
        file r1 from fastq
        file r2 from fastq
    output:
        file "${s}.ssake.fa.gz" into ssakeFasta

    """
    ngs-fastq-to-ssake -1 ${r1} -2 ${r2} -o ${s}.ssake.fa.gz --insert-size 500
    """
}

sample = Channel.from(params.sample)

process reformat {
    input:
        val s from sample
        file f from ssakeFasta
    output:
        file "${s}.ssake.reformatted.fa.gz" into reformatted

    """
    gzcat $f | sed -e 's/a/A/g' -e 's/c/C/g' -e 's/g/G/g' -e 's/t/T/g' -e 's/n/N/g' -e 's/^>.*/>reformatted:500/' | gzip -c > ${s}.ssake.reformatted.fa.gz
    """
}

sample = Channel.from(params.sample)

process ssake {
    input:
        val s from sample
        file f from reformatted
    output:
        file "${s}.ssake.d" into assembled

    """
    gunzip -c ${f} > ${f}.tmp
    mkdir ${s}.ssake.d
    SSAKE -f ${f}.tmp -b ${s}.ssake.d/${s} -w 1 -h 1 -p 1 -m 50 -o 30 -c 1 -e 0.90 -k 4 -a 0.1 -x 20
    rm ${f}.tmp
    """
}

assembled.subscribe() {
    println "assembled $it"
}
