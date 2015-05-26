flow
====

Consensus assembly and variant calling workflow.


###Hacking flow

Install [Nextflow](http://www.nextflow.io/)
```bash
$ curl -fsSL get.nextflow.io | bash 

      N E X T F L O W
      Version 0.13.5 build 2985
      last modified 07-05-2015 12:30 UTC (07:30 CDT)
      http://nextflow.io

Nextflow installation completed.
```

Run workflow locally (requires workflow dependencies to be installed locally)
```bash
$ ./nextflow run main.nf 
N E X T F L O W  ~  version 0.13.5
Launching main.nf
[warm up] executor > local
[83/654161] Submitted process > interleave (example)
[6b/88ef80] Submitted process > fastqToSsake (example)
[50/c3d689] Submitted process > reformat (example)
[ee/21d368] Submitted process > alignReads (example)
[47/a1c58f] Submitted process > ssake (example)
Copying example.reads.bwa.sorted.bam (readsBam) into: tutorial/final
Copying example.reads.bwa.sorted.vcf.gz (readsVcf) into: tutorial/final
[a4/e1a618] Submitted process > alignContigs (example)
Copying example.contigs (contigsFasta) into: tutorial/final
Copying example.contigs.bwa.sorted.bam (contigsBam) into: tutorial/final
Copying example.contigs.bwa.sorted.vcf.gz (configVcf) into: tutorial/final
```

Run workflow locally using Docker image (requires only Docker to be installed locally)
```bash
$ ./nextflow run main.nf -with-docker nmdpbioinformatics/docker-flow:1.0-snapshot-2
N E X T F L O W  ~  version 0.13.5
Launching main.nf
[warm up] executor > local
[ad/d7bd89] Submitted process > fastqToSsake (example)
[a7/b6cf05] Submitted process > interleave (example)
[ee/86ea16] Submitted process > reformat (example)
[ec/3370df] Submitted process > alignReads (example)
[3c/f4bbf7] Submitted process > ssake (example)
Copying example.reads.bwa.sorted.bam (readsBam) into: tutorial/final
Copying example.reads.bwa.sorted.vcf.gz (readsVcf) into: tutorial/final
[1b/fbdffd] Submitted process > alignContigs (example)
Copying example.contigs (contigsFasta) into: tutorial/final
Copying example.contigs.bwa.sorted.bam (contigsBam) into: tutorial/final
Copying example.contigs.bwa.sorted.vcf.gz (configVcf) into: tutorial/final
```

Use [SLURM](https://computing.llnl.gov/linux/slurm/) executor, with or without Docker
```bash
$ echo "process.executor = 'slurm'" > nextflow.config
$ ./nextflow run main.nf -with-docker nmdpbioinformatics/docker-flow:1.0-snapshot-2
N E X T F L O W  ~  version 0.13.5
[warm up] executor > slurm
[22/882358] Submitted process > fastqToSsake (example)
[16/d3ac82] Submitted process > interleave (example)
[b7/9cbbaf] Submitted process > reformat (example)
[43/2f652d] Submitted process > alignReads (example)
[ef/270f79] Submitted process > ssake (example)
Copying example.reads.bwa.sorted.bam (readsBam) into: tutorial/final
Copying example.reads.bwa.sorted.vcf.gz (readsVcf) into: tutorial/final
[18/55682c] Submitted process > alignContigs (example)
Copying example.contigs (contigsFasta) into: tutorial/final
Copying example.contigs.bwa.sorted.bam (contigsBam) into: tutorial/final
Copying example.contigs.bwa.sorted.vcf.gz (configVcf) into: tutorial/final
```

Note that when using the SLURM executor, the Nextflow working directory must be mounted on a volume shared to all SLURM nodes.