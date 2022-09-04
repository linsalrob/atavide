"""

A metagenome assembly pipeline written by Rob Edwards, Jan 2020-Jan 2022.

The idea is to start with pairs of reads and assemble those, and then 
combine all the assembled contigs and dereplicate using mmseqs, then 
map reads using bowtie2, and take the unmapped reads and reassemble.

Repeat that with the site data sets and then again with the sample
datasets.

To run on the deepthought use this command:

snakemake --configfile config/process_metagenomes.yaml -s workflow/process_metagenomes.snakefile --profile slurm

to run this on anthill use this command:
snakemake --configfile config/process_metagenomes.yaml -s workflow/process_metagenomes.snakefile --profile slurm

"""

import os
import sys
import socket
import re


READDIR = config['directories']['Reads']
PSEQDIR = config['directories']['prinseq']
STATS   = config['directories']['statistics']
TMPDIR  = config['directories']['temp_directory']
RBADIR  = config['directories']['read_based_annotations']
METABAT = os.path.join(config['directories']['binning'], 'metabat')
CONCOCT = os.path.join(config['directories']['binning'], 'concoct')
ATAVIDE_BINNING = os.path.join(config['directories']['binning'], 'atavide')

# assembly directories
ASSDIR  = os.path.join(config['directories']['assemblies'], "assembly.1")
CRMDIR  = os.path.join(config['directories']['assemblies'], "reads.contigs.1")
UNASSM  = os.path.join(config['directories']['assemblies'], "unassembled_reads")
REASSM  = os.path.join(config['directories']['assemblies'], "reassembled_reads")
CCMO    = os.path.join(config['directories']['assemblies'], "final.combined_contigs")
RMRD    = os.path.join(config['directories']['assemblies'], "reads_vs_final_assemblies")

# how much memory do we have
LARGE_MEM = config['parameters']['large_mem']

SAMPLE_ID=re.sub('\W+','', config['sample_id'])


# what is the directory of atavide.snakefile.
# We need to add that to the pythonpath and also
# use it for scripts

ATAVIDE_DIR = workflow.basedir
# append to pythonpath
sys.path.append(ATAVIDE_DIR)

# Where is the taxonomy data
TAXON = None
if 'ncbi_taxonomy' in config['directories']:
    TAXON = config['directories']['ncbi_taxonomy']
elif 'NCBI_TAXONOMY' in os.environ:
    TAXON = os.environ['NCBI_TAXONOMY']
else:
    sys.stderr.write("FATAL: Please set the location of your NCBI taxonomy data using the NCBI_TAXONOMY environment variable\n")
    sys.exit(1)


# Where is the superfocus database?
SFDB = None
if 'superfocus_db' in config['directories']:
    SFDB = config['directories']['superfocus_db']
    os.environ['SUPERFOCUS_DB'] = SFDB
elif 'SUPERFOCUS_DB' in os.environ:
    SFDB = os.environ['SUPERFOCUS_DB']
else:
    sys.stderr.write("FATAL: Please set the location of your superfocus databases using the SUPERFOCUS_DB environment variable\n")
    sys.exit(1)
# how many reads do we want to subsample for superfocus?
superfocus_subsample_reads = 0
if 'superfocus_reads' in config['parameters']:
    superfocus_subsample_reads = config['parameters']['superfocus_reads']


# do we want to do host removal?

# we use an additional PSEQDIR variables here. PSEQDIR is the regular output from prinseq
# if we are NOT doing host removal PSEQDIR AND PSEQDIR_TWO are the same
# if we ARE doing host removal, they are not the same and so we have to make the 
# second path
if 'host_dbpath' in config['directories'] and config['directories']['host_dbpath']:
    if not 'host_dbname' in config['options']:
        sys.stderr.write(f"ERROR: You have set host_dbpath as {config['directories']['host_dbpath']} but not defined the db_name\n")
        sys.exit(0)
   
    bt2l = os.path.join(config['directories']['host_dbpath'], f"{config['options']['host_dbname']}.1.bt2l")
    bt2r = os.path.join(config['directories']['host_dbpath'], f"{config['options']['host_dbname']}.1.bt2")

    if not (os.path.exists(bt2l) or os.path.exists(bt2r)):
        sys.stderr.write(f"Error: don't seem to be able to find a bowtie2 index for {config['options']['host_dbname']}\n")
        sys.stderr.write(f"\tWe looked for either of\n\t{bt2r}  or\n\t{bt2l}\n")
        sys.exit(0)

    PSEQDIR_TWO = f"{PSEQDIR}_after_hostremoval"

    PSEQDIR = f"{PSEQDIR}_before_hostremoval"

    include: "rules/deconseq.snakefile"
else:
    PSEQDIR_TWO = PSEQDIR

os.makedirs(PSEQDIR_TWO, exist_ok=True)

# A Snakemake regular expression matching the forward mate FASTQ files.
# the comma after SAMPLES is important!
SAMPLES,EXTENSIONS, = glob_wildcards(os.path.join(READDIR, '{sample}_R1{extn}'))
if len(SAMPLES) == 0:
    sys.stderr.write(f"We did not find any fastq files in {SAMPLES}. Is this the right read dir?\n")
    sys.exit(0)
if len(set(EXTENSIONS)) != 1:
    sys.stderr.write("FATAL: You have more than one type of file extension\n\t")
    sys.stderr.write("\n\t".join(set(EXTENSIONS)))
    sys.stderr.write("\nWe don't know how to handle these\n")
    sys.exit(0)

FQEXTN = EXTENSIONS[0]
PATTERN_R1 = '{sample}_R1' + FQEXTN
PATTERN_R2 = '{sample}_R2' + FQEXTN

# read the rules for running different pieces and parts of the code
include: "rules/qc_qa.snakefile"
include: "rules/focus.snakefile"
include: "rules/superfocus.snakefile"
include: "rules/round1_assembly.snakefile"
include: "rules/compress_outputs.snakefile"
include: "rules/round2_assembly.snakefile"
include: "rules/statistics.snakefile"
include: "rules/binning.snakefile"
include: "rules/kraken_taxonomy.snakefile"
include: "rules/kraken_rarefaction.snakefile"
include: "rules/singlem.snakefile"
include: "rules/combine_read_annotations.snakefile"
include: "rules/atavide_clusters.snakefile"
include: "rules/fqchk.snakefile"

localrules: all

rule all:
    input:
        expand(
            [
                os.path.join(PSEQDIR_TWO, "{sample}_good_out_R1.fastq.gz"),
                os.path.join(RBADIR, "{sample}", "focus", "output_All_levels.csv.zip"),
                os.path.join(RBADIR, "{sample}", "superfocus", "output_all_levels_and_function.xls.zip"),
                os.path.join(RBADIR, "{sample}", "superfocus", "{sample}_good_out_R1.taxonomy.zip"),
                os.path.join(RBADIR, "{sample}", "kraken", "{sample}.report.tsv.zip"),
                os.path.join(RBADIR, "{sample}", "kraken", "{sample}.output.tsv.zip"),
                os.path.join(RBADIR, "{sample}", "kraken", "{sample}.taxonomy.tsv.zip"), 
                os.path.join(RBADIR, "{sample}", "singlem", "singlem_otu_table.tsv.zip"),
                os.path.join(RMRD, "{sample}." + SAMPLE_ID + ".assembly.bam.bai"),
                os.path.join(RBADIR, "{sample}", "{sample}_contig_taxonomy.comparison.tsv.zip")
            ],
               sample=SAMPLES),
        os.path.join(REASSM, "merged_contigs.fa"),
        os.path.join(CCMO, "flye.log"),
        os.path.join(STATS, "post_qc_statistics.tsv"),
        os.path.join(STATS, "initial_read_statistics.tsv"),
        os.path.join(STATS, "final_assembly.txt.zip"),
        os.path.join(STATS, "sample_coverage.tsv.zip"),
        os.path.join(STATS, "sample_rpkm.tsv"),
        os.path.join(RBADIR, "superfocus_functions.tsv.gz"),
        os.path.join(STATS, "av_quality_scores_by_position.tsv"),
        os.path.join(STATS, "kraken_species_rarefaction.tsv"),
        os.path.join(STATS, "kraken_species.tsv"),
        os.path.join(STATS, "kraken_phyla.tsv"),
        os.path.join(STATS, "kraken_families.tsv"),
        os.path.join(STATS, "kraken_genera.tsv"),
        os.path.join(METABAT, "metabat_depth"),
        os.path.join(METABAT, "metabat_bins/metabat_bins.1.fa"),
        os.path.join(CONCOCT, "concoct_output/clustering_gt1000.csv"),
        os.path.join(CONCOCT, "concoct_bins/1.fa"),
        os.path.join(ATAVIDE_BINNING, "stats", "sample_coverage.h5"),
        os.path.join(ATAVIDE_BINNING, "stats", "atavide_clusters.json"),
        os.path.join(ATAVIDE_BINNING, "bins"),


