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
[22/8312f0] Submitted process > interleave (1)
[9d/3cfe76] Submitted process > fastqToSsake (1)
[77/c60d3b] Submitted process > alignReads (1)
[90/395c41] Submitted process > reformat (1)
[ad/3eae51] Submitted process > ssake (1)
[13/6d2035] Submitted process > copyReadsBam (1)
[73/c3c8ba] Submitted process > copyReadsVcf (1)
[ab/32323b] Submitted process > alignContigs (1)
[82/ca7732] Submitted process > copyContigsFasta (1)
[fd/6e470f] Submitted process > copyContigsBam (1)
[21/cc12ce] Submitted process > copyContigsVcf (1)
```

Run workflow locally using Docker image (requires only Docker to be installed locally)
```bash
$ ./nextflow run main.nf -with-docker nmdpbioinformatics/docker-flow:1.0-snapshot-2
N E X T F L O W  ~  version 0.13.5
Launching main.nf
[warm up] executor > local
[71/ab4304] Submitted process > fastqToSsake (1)
[fb/b1bf04] Submitted process > interleave (1)
[0f/f5049d] Submitted process > reformat (1)
[a0/68b5e2] Submitted process > alignReads (1)
[40/b8f680] Submitted process > ssake (1)
[f7/459a41] Submitted process > copyReadsBam (1)
[16/b40f02] Submitted process > copyReadsVcf (1)
[37/7f3f64] Submitted process > alignContigs (1)
[35/96759b] Submitted process > copyContigsBam (1)
[de/12b07b] Submitted process > copyContigsFasta (1)
[75/0deafc] Submitted process > copyContigsVcf (1)
```

Use [SLURM](https://computing.llnl.gov/linux/slurm/) executor, with or without Docker
```bash
$ echo "process.executor = 'slurm'" > nextflow.config
$ ./nextflow run main.nf -with-docker nmdpbioinformatics/docker-flow:1.0-snapshot-2
N E X T F L O W  ~  version 0.13.5
[warm up] executor > slurm
[be/99c179] Submitted process > interleave (1)
[d4/f1f89d] Submitted process > fastqToSsake (1)
[9f/011a75] Submitted process > alignReads (1)
[ff/2f03d9] Submitted process > reformat (1)
[60/996a78] Submitted process > ssake (1)
[86/acdaab] Submitted process > copyReadsBam (1)
[ac/2897ed] Submitted process > copyReadsVcf (1)
[75/7f9983] Submitted process > alignContigs (1)
[16/664353] Submitted process > copyContigsBam (1)
[8e/82be94] Submitted process > copyContigsFasta (1)
[27/a1c267] Submitted process > copyContigsVcf (1)
```

Note that when using the SLURM executor, the Nextflow working directory must be mounted on a volume shared to all SLURM nodes.