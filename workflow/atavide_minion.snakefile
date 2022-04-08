"""

A metagenome assembly pipeline written by Rob Edwards, Jan 2020-Mar 2022.

This version is a trimmed down version that should accept minion files.

The goal is to have your samples as:

    fastq/barcode01.fastq
    fastq/barcode02.fastq
    fastq/barcode03.fastq
    ...

and then we will process them all. Note at the moment we combine all these
in the assemblies, but that may not be what you really want.

NOTE: Need to implement .gz reading in superfocus/focus

To run on the deepthought use this command:

snakemake --configfile config/process_metagenomes.yaml -s workflow/atavide_minion.snakefile --profile slurm


"""

import os
import sys
import socket
import re


READDIR = config['directories']['Reads']
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



# do we want to do host removal?

# we use an additional PSEQDIR variables here. PSEQDIR is the regular output from prinseq
# if we are NOT doing host removal PSEQDIR AND PSEQDIR_TWO are the same
# if we ARE doing host removal, they are not the same and so we have to make the 
# second path
if 'host_dbpath' in config['directories'] and config['directories']['host_dbpath']:
    if not 'host_dbname' in config['options']:
        sys.stderr.write(f"ERROR: You have set host_dbpath as {config['directories']['host_dbpath']} but not defined the db_name\n")
        sys.exit(0)
   
    if not os.path.exists(
        os.path.join(config['directories']['host_dbpath'], f"{config['options']['host_dbname']}.1.bt2l") or
        os.path.join(config['directories']['host_dbpath'], f"{config['options']['host_dbname']}.1.bt2")
    ):
        sys.stderr.write(f"Error: don't seem to be able to find a bowtie2 index for {config['options']['host_dbname']}\n")
        sys.exit(0)

    PSEQDIR_TWO = f"after_hostremoval"
    PSEQDIR = READDIR
    include: "rules/deconseq_minion.snakefile"
    os.makedirs(PSEQDIR_TWO, exist_ok=True)
else:
    PSEQDIR_TWO = READDIR 


# A Snakemake regular expression matching the forward mate FASTQ files.
# the comma after SAMPLES is important!
SAMPLES, = glob_wildcards(os.path.join(READDIR, '{sample}.fastq'))
if len(SAMPLES) == 0:
    sys.stderr.write(f"We did not find any fastq files in {SAMPLES}. Is this the right read dir?\n")
    sys.exit(0)

FASTQ = '{sample}.fastq'

# read the rules for running different pieces and parts of the code
include: "rules/atavide_clusters.snakefile"
include: "rules/binning.snakefile"
include: "rules/combine_read_annotations.snakefile"
include: "rules/compress_outputs.snakefile"
include: "rules/focus_minion.snakefile"
include: "rules/fqchk_minion.snakefile"
include: "rules/kraken_rarefaction_minion.snakefile"
include: "rules/kraken_taxonomy_minion.snakefile"
include: "rules/qc_qa_minion.snakefile"
include: "rules/round1_assembly_minion.snakefile"
include: "rules/round2_assembly_minion.snakefile"
include: "rules/singlem_minion.snakefile"
include: "rules/statistics.snakefile"
include: "rules/superfocus_minion.snakefile"



rule all:
    input:
        expand(
            [
                os.path.join(RBADIR, "{sample}", "focus", "output_All_levels.csv.zip"),
                os.path.join(RBADIR, "{sample}", "superfocus", "output_all_levels_and_function.xls.zip"),
                os.path.join(RBADIR, "{sample}", "superfocus", "{sample}.taxonomy.zip"),
                os.path.join(RBADIR, "{sample}", "kraken", "{sample}.report.tsv.zip"),
                os.path.join(RBADIR, "{sample}", "kraken", "{sample}.output.tsv.zip"),
                os.path.join(RBADIR, "{sample}", "kraken", "{sample}.taxonomy.tsv"), 
                os.path.join(RBADIR, "{sample}", "singlem", "singlem_otu_table.tsv.zip"),
                os.path.join(RMRD, "{sample}." + SAMPLE_ID + ".assembly.bam.bai"),
                os.path.join(RBADIR, "{sample}", "{sample}_contig_taxonomy.comparison.tsv")
            ],
               sample=SAMPLES),
        os.path.join(REASSM, "merged_contigs.fa"),
        os.path.join(CCMO, "flye.log"),
        os.path.join(STATS, "final_assembly.txt.zip"),
        os.path.join(STATS, "sample_coverage.tsv.zip"),
        os.path.join(STATS, "sample_rpkm.tsv"),
        os.path.join(RBADIR, "superfocus_functions.tsv.gz"),
        os.path.join(STATS, "av_quality_scores_by_position.tsv"),
        os.path.join(STATS, "kraken_species_rarefaction.tsv"),
        os.path.join(STATS, "kraken_families.tsv"),
        os.path.join(METABAT, "metabat_depth"),
        os.path.join(METABAT, "metabat_bins/metabat_bins.1.fa"),
        os.path.join(CONCOCT, "concoct_output/clustering_gt1000.csv"),
        os.path.join(CONCOCT, "concoct_bins/1.fa"),
        os.path.join(ATAVIDE_BINNING, "stats", "sample_coverage.h5"),
        os.path.join(ATAVIDE_BINNING, "stats", "atavide_clusters.json"),
        os.path.join(ATAVIDE_BINNING, "bins"),


