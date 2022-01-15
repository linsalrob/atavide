[![Edwards Lab](https://img.shields.io/badge/Bioinformatics-EdwardsLab-03A9F4)](https://edwards.flinders.edu.au)
[![DOI](https://www.zenodo.org/badge/403921714.svg)](https://www.zenodo.org/badge/latestdoi/403921714)
[![DOI](https://img.shields.io/badge/DOI-WorkflowHub-yellowgreen)](https://doi.org/10.48546/WORKFLOWHUB.WORKFLOW.241.1)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![GitHub language count](https://img.shields.io/github/languages/count/linsalrob/atavide)


# atavide

`atavide` is a simple, yet complete workflow for metagenomics data analysis, including QC/QA, optional host removal, assembly and cross-assembly, and individual read based annotations. We have also built in some advanced analytics including tools to assign annotations from reads to contigs, and to generate metagenome-assembled genomes in several different ways, giving you the power to explore your data!

`atavide` is 100% snakemake and conda, so you only need to install the snakemake workflow, and then everything else will be installed with conda.


It is definitely a work in progress, but you can run it with the following command 

```bash
snakemake --configfile config/atavide.yaml -s workflow/atavide.snakefile --profile slurm
```

But you will need a [slurm profile](https://fame.flinders.edu.au/blog/2021/08/02/snakemake-profiles-updated) to make this work!


# Installation and getting going


1. Clone this repository from GitHub: `git clone https://github.com/linsalrob/atavide.git`
2. Set the location of the repository: `export ATAVIDE_DIR=$PWD/atavide/`
3. Install a few python dependencies. You probably already have most of these, but the one that trips up is `pysam`. We're working on getting `conda` configured properly to do this automatically. `pip install -r $ATAVIDE_DIR/requirements.txt`
4. Install the [appropriate super-focus database](https://github.com/metageni/SUPER-FOCUS/issues/66) [hint: probably version 2] and set the `SUPERFOCUS_DB` directory to [point to the location of those files](https://github.com/metageni/SUPER-FOCUS#database).
5. Copy the [NCBI taxonomy](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/) (You really just need the [taxdump.tar.gz](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz) file), and set the `NCBI_TAXONOMY` environment variable to point to the location of those files.
6. Have a directory of fastq files with both `_R1_` and `_R2_` files in a data directory: `$DATA_DIR/fastq` 
7. Run atavide: `cd $DATA_DIR && snakemake --configfile $ATAVIDE_DIR/atavide.yaml -s $ATAVIDE_DIR/workflow/atavide.snakefile --profile slurm`


# Current processing steps:

### Steps:
1. QC/QA with [prinseq++](https://github.com/Adrian-Cantu/PRINSEQ-plus-plus)
2. optional host removal using bowtie2 and samtools, [as described previously](https://edwards.flinders.edu.au/command-line-deconseq/). To enable this, you need to provide a path to the host db and a host db.

### Metagenome assembly
1. pairwise assembly of each sample using [megahit](https://github.com/voutcn/megahit)
2. extraction of all reads that do not assemble using samtools flags
3. assembly of all unassembled reads using [megahit](https://github.com/voutcn/megahit)
4. compilation of _all_ contigs into a single unified set using [Flye](https://github.com/fenderglass/Flye)
5. comparison of reads -> contigs to generate coverage

### MAG creation
1. [metabat](https://bitbucket.org/berkeleylab/metabat/src/master/)
2. [concoct](https://github.com/BinPro/CONCOCT)
3. Pairwise comparisons using [turbocor](https://github.com/dcjones/turbocor) followed by clustering

### Read-based annotations
1. [Kraken2](https://ccb.jhu.edu/software/kraken2/)
2. [singlem](https://github.com/wwood/singlem)
3. [SUPER-focus](https://github.com/metageni/SUPER-FOCUS)
4. [FOCUS](https://github.com/metageni/FOCUS)

Want something else added to the suite? File an issue on GitHub and we'll add it ASAP!


