flow
====

Consensus assembly, variant calling, and allele interpretation workflow.


###Hacking ngs

Install [Nextflow](http://www.nextflow.io/)
```bash
$ curl -fsSL get.nextflow.io | bash 

      N E X T F L O W
      Version 0.12.5 build 2800
      last modified 18-03-2015 09:11 UTC (04:11 CDT)
      http://nextflow.io

Nextflow installation completed.
```

Run assembly
```bash
$ ./nextflow assemble.nf 
N E X T F L O W  ~  version 0.12.5
[warm up] executor > local
[94/cfa72d] Submitted process > fastqToSsake (1)
[bb/73c436] Submitted process > interleave (1)
[48/019fbb] Submitted process > alignReads (1)
[2e/491ce1] Submitted process > reformat (1)
[f6/3c6d74] Submitted process > ssake (1)
reads vcf [example, work/48/019fbbf84290466e2223373c9a3af3/example.reads.bwa.sorted.vcf.gz]
aligned reads [example, work/48/019fbbf84290466e2223373c9a3af3/example.reads.bwa.sorted.bam]
[2e/a62ff3] Submitted process > alignContigs (1)
aligned contigs [example, work/2e/a62ff3adc7c7d73845fa8c58c47f5f/example.contigs.bwa.sorted.bam]
contigs vcf [example, work/2e/a62ff3adc7c7d73845fa8c58c47f5f/example.contigs.bwa.sorted.vcf.gz]
```

Use [SLURM](https://computing.llnl.gov/linux/slurm/) executor
```bash
$ echo "process.executor = 'slurm'" > netflow.config
$ ./nextflow assemble.nf 
N E X T F L O W  ~  version 0.12.5
[warm up] executor > slurm
[ab/8b26c9] Submitted process > interleave (1)
[51/19b1c5] Submitted process > fastqToSsake (1)
[38/f27ca8] Submitted process > alignReads (1)
aligned reads [example, work/38/f27ca8187337e548ec11c9d5003e99/example.reads.bwa.sorted.bam]
reads vcf [example, work/38/f27ca8187337e548ec11c9d5003e99/example.reads.bwa.sorted.vcf.gz]
[d9/0ef973] Submitted process > reformat (1)
[61/4fae5f] Submitted process > ssake (1)
[81/e19095] Submitted process > alignContigs (1)
aligned contigs [example, work/81/e19095c5509793edf9d2047ffc4ecf/example.contigs.bwa.sorted.bam]
contigs vcf [example, work/81/e19095c5509793edf9d2047ffc4ecf/example.contigs.bwa.sorted.vcf.gz]
```

Note that when using the SLURM executor, the Nextflow working directory must be mounted on a volume shared to all SLURM nodes.