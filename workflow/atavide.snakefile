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


READDIR = config['directories']['Reads']
PSEQDIR = config['directories']['prinseq']
STATS   = config['directories']['statistics']
RBADIR  = config['directories']['read_based_annotations']
METABAT = config['directories']['metabat']
CONCOCT = config['directories']['concoct']
TAXON   = config['directories']['ncbi_taxonomy']


# atavide binning
ASSDIR  = os.path.join(config['directories']['atavide_binning'], "assembly.1")
CRMDIR  = os.path.join(config['directories']['atavide_binning'], "reads.contigs.1")
UNASSM  = os.path.join(config['directories']['atavide_binning'], "unassembled_reads")
REASSM  = os.path.join(config['directories']['atavide_binning'], "reassembled_reads")
CCMO    = os.path.join(config['directories']['atavide_binning'], "final.combined_contigs")
RMRD    = os.path.join(config['directories']['atavide_binning'], "reads_vs_final_assemblies")

# what is the directory of atavide.snakefile.
# We need to add that to the pythonpath and also
# use it for scripts

ATAVIDE_DIR = workflow.basedir
# append to pythonpath
sys.path.append(ATAVIDE_DIR)

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

    PSEQDIR_TWO = f"{PSEQDIR}_after_hostremoval"

    PSEQDIR = f"{PSEQDIR}_before_hostremoval"

    include: "rules/deconseq.snakefile"
else:
    PSEQDIR_TWO = PSEQDIR


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
include: "rules/read_annotations.snakefile"
include: "rules/round1_assembly.snakefile"
include: "rules/compress_outputs.snakefile"
include: "rules/round2_assembly.snakefile"
include: "rules/statistics.snakefile"
include: "rules/binning.snakefile"
include: "rules/kraken_taxonomy.snakefile"
include: "rules/singlem.snakefile"
include: "rules/combine_read_annotations.snakefile"

rule all:
    input:
        # these rules are from rules/read_annotations.snakefile
        expand(
            [
                os.path.join(PSEQDIR_TWO, "{sample}_good_out_R1.fastq"),
                os.path.join(RBADIR, "{sample}", "focus", "output_All_levels.csv.zip"),
                os.path.join(RBADIR, "{sample}", "superfocus", "output_all_levels_and_function.xls.zip"),
                os.path.join(RBADIR, "{sample}", "superfocus", "{sample}_good_out.taxonomy"),
                os.path.join(RBADIR, "{sample}", "kraken", "{sample}.report.tsv.zip"),
                os.path.join(RBADIR, "{sample}", "kraken", "{sample}.output.tsv.zip"),
                os.path.join(RBADIR, "{sample}", "kraken", "{sample}.taxonomy.tsv"), 
                os.path.join(RBADIR, "{sample}", "singlem", "singlem_otu_table.tsv.zip"),
                os.path.join(RMRD, "{sample}.final_contigs.bam.bai"),
                os.path.join(RBADIR, "{sample}", "{sample}_contig_taxonomy.comparison.tsv")
            ],
               sample=SAMPLES),
        os.path.join(RBADIR, "superfocus_functions.tsv.zip"),
        os.path.join(REASSM, "merged_contigs.fa"),
        os.path.join(CCMO, "flye.log"),
        os.path.join(STATS, "final_assembly.txt.zip"),
        os.path.join(STATS, "sample_coverage.tsv.zip"),
        os.path.join(STATS, "sample_coverage.h5"),
        os.path.join(METABAT, "metabat_depth"),
        os.path.join(METABAT, "metabat_bins/metabat_bins.1.fa"),
        os.path.join(CONCOCT, "concoct_output/clustering_gt1000.csv"),
        os.path.join(CONCOCT, "concoct_bins/1.fa")


