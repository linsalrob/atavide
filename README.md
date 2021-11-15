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


1. Clone this repository from GitHub: `git clone https://github.com/linsalrob/atavide.git`
2. Set the location of the repository: `export ATAVIDE_DIR=$PWD/atavide/`
2. Install the [appropriate super-focus database](https://github.com/metageni/SUPER-FOCUS/issues/66) [hint: probably version 2] and set the `SUPERFOCUS_DB` directory to [point to the location of those files](https://github.com/metageni/SUPER-FOCUS#database).
3. Have a directory of fastq files with both `_R1_` and `_R2_` files in a data directory: `$DATA_DIR/fastq` 
4. Run atavide: `cd $DATA_DIR && snakemake --configfile $ATAVIDE_DIR/atavide.yaml -s $ATAVIDE_DIR/workflow/atavide.snakefile --profile slurm`


# Current processing steps:


1. Assemble each pair of reads (`_R1_` and `_R2_` separately)
2. Merge those contigs
3. Map all reads to those contigs and identify unassembled reads
4. Assemble all of the unassembled reads
5. Merge the contigs from the unassembled reads with the initial contigs to get a final assembly
6. Map all the reads back to the final assembly
7. Run pairwise correlations between the number of reads mapped to each contig to identify similar contigs
8. Generate atavide bins from the pairwise mapping.
9. Generate concoct bins from the final contigs and mapped reads
10. Generate metabat bins from the final contigs and mapped reads
11. Run read-based annotations:
    - focus
    - super-focus
    - kraken
    - singlem
12. Combine all the read based annotations with the contig-read mappings to get contig level annotations






