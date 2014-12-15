flow
====

Consensus assembly, variant calling, and allele interpretation workflow.


###Hacking ngs

Install [Nextflow](http://www.nextflow.io/)
```bash
$ curl -fsSL get.nextflow.io | bash 
```

Run assembly
```bash
$ ./nextflow run assemble.nf
N E X T F L O W  ~  version 0.11.3
[warm up] executor > local
[0f/289ed2] Submitted process > fastqToSsake (1)
[29/2ccb7f] Submitted process > reformat (1)
[27/f277f6] Submitted process > ssake (1)
assembled /tmp/work/27/f277f6d0e48ae00f491108d7e63ec5/example.ssake.d
```

Use [SLURM](https://computing.llnl.gov/linux/slurm/) executor
```bash
$ echo "process.executor = 'slurm'" > netflow.config
$ ./nextflow run assemble.nf
N E X T F L O W  ~  version 0.11.3
[warm up] executor > slurm
[64/2d8202] Submitted process > fastqToSsake (1)
[c8/2ee769] Submitted process > reformat (1)
[f4/6aa524] Submitted process > ssake (1)
assembled /tmp/work/f4/6aa524b5a20bcf385d11f88e37cdf8/example.ssake.d
```