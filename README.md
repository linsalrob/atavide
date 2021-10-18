[![Edwards Lab](https://img.shields.io/badge/Bioinformatics-EdwardsLab-03A9F4)](https://edwards.flinders.edu.au)
[![DOI](https://www.zenodo.org/badge/403921714.svg)](https://www.zenodo.org/badge/latestdoi/403921714)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![GitHub language count](https://img.shields.io/github/languages/count/linsalrob/atavide)


# atavide

A simple metagenome processing pipeline that is designed to clean the data, assemble it into MAGs, and do something interesting with it.

It is definitely a work in progress, but you can run it with the following command 

```bash
snakemake --configfile config/atavide.yaml -s workflow/atavide.snakefile --profile slurm
```

But you will need a [slurm profile](https://fame.flinders.edu.au/blog/2021/08/02/snakemake-profiles-updated) to make this work!


# Installation and getting going

Atavide is mostly self-installing. We rely on conda for most of what we do. However, if you wish to run super-focus, please [download the appropriate super-focus database](https://github.com/metageni/SUPER-FOCUS/issues/66) [hint: probably version 2]
and set the `SUPERFOCUS_DB` directory to [point to the location of those files](https://github.com/metageni/SUPER-FOCUS#database).

Kraken2, singlem, and focus have their own datasets downloaded on install time.

To get atavide, clone it from git:

```
git clone git@github.com:linsalrob/atavide.git
ATAVIDE_DIR=$PWD/atavide/
```

Then, if your data is in a directory, `$DATA_DIR` with some fastq files in a directory `$DATA_DIR/fastq` you can use:

```
cd $DATA_DIR
snakemake --configfile $ATAVIDE_DIR/atavide.yaml -s $ATAVIDE_DIR/workflow/atavide.snakefile --profile slurm
```




