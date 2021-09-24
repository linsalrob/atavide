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

#configfile: "config/process_metagenomes.yaml"

READDIR = config['directories']['Reads']
ASSDIR  = config['directories']['round1_assembly_output']
CRMDIR  = config['directories']['round1_contig_read_mapping']
UNASSM  = config['directories']['round2_unassembled_reads']
REASSM  = config['directories']['round2_assembly_output']
CCMO    = config['directories']['combined_contig_merging']
RMRD    = config['directories']['reads_vs_final_assemblies']
PSEQDIR = config['directories']['prinseq']
STATS   = config['directories']['statistics']
RBADIR  = config['directories']['read_based_annotations']

# what is the directory of process_metagenomes.snakefile.
# We need to add that to the pythonpath and also
# use it for scripts

PMSDIR = workflow.basedir
# append to pythonpath
sys.path.append(PMSDIR)

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
SAMPLES,EXTENSIONS, = glob_wildcards(os.path.join(READDIR, '{sample}_R1.{extn}'))
if len(SAMPLES) == 0:
    sys.stderr.write(f"We did not find any fastq files in {SAMPLES}. Is this the right read dir?\n")
    sys.exit(0)
if len(set(EXTENSIONS)) != 1:
    sys.stderr.write("FATAL: You have more than one type of file extension\n\t")
    sys.stderr.write("\n\t".join(set(EXTENSIONS)))
    sys.stderr.write("\nWe don't know how to handle these\n")
    sys.exit(0)

FQEXTN = EXTENSIONS[0]
PATTERN_R1 = '{sample}_R1.' + FQEXTN
PATTERN_R2 = '{sample}_R2.' + FQEXTN

# read the rules for running kraken, focus, superfocus, singlem, etc. 
include: "rules/read_annotations.snakefile"

rule all:
    input:
        # these rules are from rules/read_annotations.snakefile
        expand(
            [
                os.path.join(PSEQDIR_TWO, "{sample}_good_out_R1.fastq"),
                os.path.join(RBADIR, "{sample}", "focus", "output_All_levels.csv"),
                os.path.join(RBADIR, "{sample}", "superfocus", "output_all_levels_and_function.xls"), ## TODO note we can not use superfocus at the moment until the environmental variable stuff is incorporated
                os.path.join(RBADIR, "{sample}", "kraken", "{sample}.report.tsv"),
                os.path.join(RBADIR, "{sample}", "kraken", "{sample}.output.tsv"),
                os.path.join(RBADIR, "{sample}", "singlem", "singlem_otu_table.tsv"),
                os.path.join(RMRD, "{sample}.final_contigs.bam.bai")
            ],
               sample=SAMPLES),
        os.path.join(RBADIR, "superfocus_functions.tsv"),
        os.path.join(RBADIR, "superfocus_taxonomy.tsv"),
        os.path.join(REASSM, "merged_contigs.fa"),
        os.path.join(CCMO, "flye.log"),
        os.path.join(STATS, "final_assembly.txt"),
        os.path.join(STATS, "sample_coverage.tsv")



rule prinseq:
    input:
        r1 = os.path.join(READDIR, PATTERN_R1),
        r2 = os.path.join(READDIR, PATTERN_R2)
    output:
        r1 = os.path.join(PSEQDIR, "{sample}_good_out_R1.fastq"),
        r2 = os.path.join(PSEQDIR, "{sample}_good_out_R2.fastq"),
        s1 = os.path.join(PSEQDIR, "{sample}_single_out_R1.fastq"),
        s2 = os.path.join(PSEQDIR, "{sample}_single_out_R2.fastq"),
        b1 = temporary(os.path.join(PSEQDIR, "{sample}_bad_out_R1.fastq")),
        b2 = temporary(os.path.join(PSEQDIR, "{sample}_bad_out_R2.fastq"))
    conda: "envs/prinseq.yaml"
    params:
        o = os.path.join(PSEQDIR, "{sample}")
    shell:
        """
            prinseq++ -min_len 60 -min_qual_mean 25 -ns_max_n 1 -derep 1 \
                    -out_format 0 -trim_tail_left 5 -trim_tail_right 5 \
                    -ns_max_n 5  -trim_qual_type min -trim_qual_left 30 \
                    -trim_qual_right 30 -trim_qual_window 10 \
                    -threads {threads} \
                    -out_name {params.o} \
                    -fastq {input.r1} \
                    -fastq2 {input.r2};
        """

# Here, we need to add an option to integrate rules/deconseq.snakefile
# so that we can move from {sample}_good_out_R1.fastq to 
# os.path.join(outdir, '{sample}_singletons_' + hostname + '.unmapped.fastq')


rule megahit_assemble:
    input:
        r1 = os.path.join(PSEQDIR_TWO, "{sample}_good_out_R1.fastq"),
        r2 = os.path.join(PSEQDIR_TWO, "{sample}_good_out_R2.fastq"),
        s1 = os.path.join(PSEQDIR_TWO, "{sample}_single_out_R1.fastq"),
        s2 = os.path.join(PSEQDIR_TWO, "{sample}_single_out_R2.fastq")
    output:
        os.path.join(ASSDIR, "{sample}/final.contigs.fa"),
        os.path.join(ASSDIR, "{sample}/log"),
        temporary(os.path.join(ASSDIR, "{sample}/checkpoints.txt")),
        temporary(os.path.join(ASSDIR, "{sample}/done")),
        temporary(directory(os.path.join(ASSDIR, "{sample}/intermediate_contigs"))),
        temporary(os.path.join(ASSDIR, "{sample}/options.json"))
    params:
        odir = directory(os.path.join(ASSDIR, '{sample}'))
    resources:
        mem_mb=64000,
        cpus=16,
        time=7200
    conda:
        "envs/megahit.yaml"
    shell:
        """
        rmdir {params.odir} ; 
        megahit -1 {input.r1} -2 {input.r2} -r {input.s1} -r {input.s2} \
                -o {params.odir} -t {resources.cpus}
        """

rule combine_contigs:
    """
    Here we take all contigs produced by megahit and combine them into a single (redundant)
    fasta file. This also creates a file called round1_contigs.ids that has the original
    contig ids
    """

    input:
        expand(os.path.join(ASSDIR, "{sample}/final.contigs.fa"), sample=SAMPLES)
    output:
        contigs = os.path.join(ASSDIR, "round1_contigs.fa"),
        ids = os.path.join(ASSDIR, "round1_contigs.ids")
    params:
        sct = os.path.join(PMSDIR, "scripts/renumber_merge_fasta.py")
    shell:
        """
        python3 {params.sct} -f {input} -o {output.contigs} -i {output.ids} -v
        """


rule index_contigs:
    """
    build the bowtie2 index of the contigs
    """
    input:
        os.path.join(ASSDIR, "round1_contigs.fa")
    params:
        baseoutput = os.path.join(CRMDIR, "round1_contigs")
    output:
        idx1 = temporary(os.path.join(CRMDIR, "round1_contigs.1.bt2l")),
        idx2 = temporary(os.path.join(CRMDIR, "round1_contigs.2.bt2l")),
        idx3 = temporary(os.path.join(CRMDIR, "round1_contigs.3.bt2l")),
        idx4 = temporary(os.path.join(CRMDIR, "round1_contigs.4.bt2l")),
        ridx1 = temporary(os.path.join(CRMDIR, "round1_contigs.rev.1.bt2l")),
        ridx2 = temporary(os.path.join(CRMDIR, "round1_contigs.rev.2.bt2l"))
    resources:
        mem_mb=64000,
        cpus=8
    conda:
        "envs/bowtie.yaml"
    shell:
        # note that we build a large index by default, because
        # at some point we will end up doing that, and so this
        # always makes a large index
        """
        mkdir -p {CRMDIR} && \
        bowtie2-build --threads {resources.cpus} --large-index {input} {params.baseoutput}
        """
        

rule map_reads:
    """
    Here we map the original reads back to the contigs
    """
    input:
        idx1 = os.path.join(CRMDIR, "round1_contigs.1.bt2l"),
        idx2 = os.path.join(CRMDIR, "round1_contigs.2.bt2l"),
        idx3 = os.path.join(CRMDIR, "round1_contigs.3.bt2l"),
        idx4 = os.path.join(CRMDIR, "round1_contigs.4.bt2l"),
        ridx1 = os.path.join(CRMDIR, "round1_contigs.rev.1.bt2l"),
        ridx2 = os.path.join(CRMDIR, "round1_contigs.rev.2.bt2l"),
        r1 = os.path.join(PSEQDIR_TWO, "{sample}_good_out_R1.fastq"),
        r2 = os.path.join(PSEQDIR_TWO, "{sample}_good_out_R2.fastq"),
        s1 = os.path.join(PSEQDIR_TWO, "{sample}_single_out_R1.fastq"),
        s2 = os.path.join(PSEQDIR_TWO, "{sample}_single_out_R2.fastq")
    params:
        contigs = os.path.join(CRMDIR, "round1_contigs")
    output:
        os.path.join(CRMDIR, "{sample}.contigs.bam")
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "envs/bowtie.yaml"
    shell:
        """
        bowtie2 --mm -x {params.contigs} -1 {input.r1} -2 {input.r2} \
        -U {input.s1} -U {input.s2} --threads {resources.cpus} | \
        samtools view -bh | samtools sort -o {output} -
        """

"""
Using samtools to extract unmapped reads from the bam files
A samtools flag of -f 77 means the read is paird, neither the read nor mate is mapped, and the it is the first read in the pair, while a flag of -f 141 means the same except it is the second mate in the pair.
Then we use two flags: -f 4 (the read is unmapped) and -F 1 (the read is not paired) to find the single reads that are not mapped

See: https://edwards.sdsu.edu/research/command-line-deconseq/
"""

rule umapped_left_reads:
    """
    get left reads
    """
    input:
        os.path.join(CRMDIR, "{sample}.contigs.bam")
    output:
        os.path.join(UNASSM, "{sample}.unassembled.R1.fastq")
    conda:
        "envs/bowtie.yaml"
    resources:
        mem_mb=32000,
        cpus=8
    shell:
        """
        samtools view -@ {resources.cpus} -h {input} | 
                awk 'BEGIN {{FS="\t"; OFS="\t"}} 
                {{if (/^@/ && substr($2, 3, 1)==":") {{print}} 
                else if (and($2, 0x1) && and($2, 0x40) && 
                (and($2, 0x4) || and($2, 0x8))) {{print}}}}' \
                | samtools bam2fq -@ {resources.cpus} > {output}
        """

rule umapped_right_reads:
    """
    get right reads
    """
    input:
        os.path.join(CRMDIR, "{sample}.contigs.bam")
    output:
        os.path.join(UNASSM, "{sample}.unassembled.R2.fastq")
    conda:
        "envs/bowtie.yaml"
    resources:
        mem_mb=32000,
        cpus=8
    shell:
        """
        samtools view -@ {resources.cpus} -h {input} | 
                awk 'BEGIN {{FS="\t"; OFS="\t"}} 
                {{if (/^@/ && substr($2, 3, 1)==":") {{print}} 
                else if (and($2, 0x1) && and($2, 0x80) && 
                (and($2, 0x4) || and($2, 0x8))) {{print}}}}' \
                | samtools bam2fq -@ {resources.cpus} > {output}
        """

rule umapped_single_reads:
    """
    get singletons
    """
    input:
        os.path.join(CRMDIR, "{sample}.contigs.bam")
    output:
        os.path.join(UNASSM, "{sample}.unassembled.singles.fastq")
    conda:
        "envs/bowtie.yaml"
    resources:
        mem_mb=32000,
        cpus=8
    shell:
            "samtools fastq -@ {resources.cpus} -f 4 -F 1 {input} > {output}"

"""
We concatanate the unassembled reads into separate R1/R2/s files so
we can assmble them all together.

We also do this in 3 separate threads to take advantage of parallelization
and to make the command easier
"""

rule concatenate_R1_unassembled:
    """
    Start with R1 reads
    """
    input:
        expand(os.path.join(UNASSM, "{sample}.unassembled.R1.fastq"), sample=SAMPLES)
    output:
        os.path.join(UNASSM, "R1.unassembled.fastq")
    shell:
        "cat {input} > {output}"

rule concatenate_R2_unassembled:
    """
    Concat R2 reads
    """
    input:
        expand(os.path.join(UNASSM, "{sample}.unassembled.R2.fastq"), sample=SAMPLES)
    output:
        os.path.join(UNASSM, "R2.unassembled.fastq")
    shell:
        "cat {input} > {output}"

rule concatenate_single_unassembled:
    """
    Concate singletons
    """
    input:
        expand(os.path.join(UNASSM, "{sample}.unassembled.singles.fastq"),sample=SAMPLES)
    output:
        os.path.join(UNASSM, "single.unassembled.fastq")
    shell:
        "cat {input} > {output}"

rule assemble_unassembled:
    """
    assemble the unassembled reads
    """
    input:
        r1 = os.path.join(UNASSM, "R1.unassembled.fastq"),
        r2 = os.path.join(UNASSM, "R2.unassembled.fastq"),
        s0 = os.path.join(UNASSM, "single.unassembled.fastq")
    output:
        os.path.join(REASSM, "final.contigs.fa"),
        os.path.join(REASSM, "log"),
        temporary(os.path.join(REASSM, "checkpoints.txt")),
        temporary(os.path.join(REASSM, "done")),
        temporary(directory(os.path.join(REASSM, "intermediate_contigs"))),
        temporary(os.path.join(REASSM, "options.json"))
    params:
        odir = os.path.join(REASSM)
    conda:
        "envs/megahit.yaml"
    resources:
        mem_mb=128000,
        cpus=32,
        time=7200
    shell:
        """
        rmdir {params.odir};
        megahit -1 {input.r1} -2 {input.r2} -r {input.s0} -o {params.odir} -t {resources.cpus}
        """

"""
Combine all the contigs and use metaflye to merge the subassmblies
"""


rule concatenate_all_assemblies:
    """
    Again we take all contigs produced by megahit and combine them into a single (redundant)
    fasta file. This also creates a file called merged_contigs.ids that has the 
    contig ids. 
    """

    input:
        expand(os.path.join(ASSDIR, "{sample}/final.contigs.fa"), sample=SAMPLES),
        os.path.join(REASSM, "final.contigs.fa")
    output:
        contigs = os.path.join(REASSM, "merged_contigs.fa"),
        ids = os.path.join(REASSM, "merged_contigs.ids") 
    params:
        sct = os.path.join(PMSDIR, "scripts/renumber_merge_fasta.py")
    shell:
        """
        python3 {params.sct} -f {input} -o {output.contigs} -i {output.ids} -v
        """

rule merge_assemblies_with_flye:
    """
    Run flye on all the merged contigs fromt the last round to merge everything one more time
    """
    input:
        contigs = os.path.join(REASSM, "merged_contigs.fa")
    output:
        os.path.join(CCMO, "assembly.fasta"),
        os.path.join(CCMO, "assembly_graph.gfa"),
        os.path.join(CCMO, "assembly_graph.gv"),
        os.path.join(CCMO, "assembly_info.txt"),
        os.path.join(CCMO, "flye.log"),
    conda:
        "envs/flye.yaml"
    resources:
        mem_mb=512000,
        cpus=32
    shell:
        """
        flye --meta --subassemblies {input.contigs} -o {CCMO} --threads {resources.cpus}
        """

rule final_assembly_stats:
    input:
        os.path.join(CCMO, "assembly.fasta")
    output:
        os.path.join(STATS, "final_assembly.txt")
    params:
        sct = os.path.join(PMSDIR, "scripts/countfasta.py")
    shell:
        """
        python3 {params.sct} -f {input} > {output}
        """


rule index_final_contigs:
    """
    build the bowtie2 index of the final assembly
    """
    input:
        os.path.join(CCMO, "assembly.fasta")
    params:
        baseoutput = os.path.join(RMRD, "final_contigs")
    output:
        idx1 = temporary( os.path.join(RMRD, "final_contigs" + ".1.bt2l")),
        idx2 = temporary( os.path.join(RMRD, "final_contigs" + ".2.bt2l")),
        idx3 = temporary( os.path.join(RMRD, "final_contigs" + ".3.bt2l")),
        idx4 = temporary( os.path.join(RMRD, "final_contigs" + ".4.bt2l")),
        rid1 = temporary( os.path.join(RMRD, "final_contigs" + ".rev.1.bt2l")),
        rid2 = temporary( os.path.join(RMRD, "final_contigs" + ".rev.2.bt2l"))
    resources:
        mem_mb=64000,
        cpus=16
    conda:
        "envs/bowtie.yaml"
    shell:
        # note that we build a large index by default, because
        # at some point we will end up doing that, and so this
        # always makes a large index
        """
        mkdir -p {RMRD} && \
        bowtie2-build --threads {resources.cpus} --large-index {input} {params.baseoutput}
        """
        

rule map_reads_to_final:
    """
    Here we map the original reads back to the contigs
    """
    input:
        idx1 = os.path.join(RMRD, "final_contigs" + ".1.bt2l"),
        idx2 = os.path.join(RMRD, "final_contigs" + ".2.bt2l"),
        idx3 = os.path.join(RMRD, "final_contigs" + ".3.bt2l"),
        idx4 = os.path.join(RMRD, "final_contigs" + ".4.bt2l"),
        rid1 = os.path.join(RMRD, "final_contigs" + ".rev.1.bt2l"),
        rid2 = os.path.join(RMRD, "final_contigs" + ".rev.2.bt2l"),
        r1 = os.path.join(PSEQDIR_TWO, "{sample}_good_out_R1.fastq"),
        r2 = os.path.join(PSEQDIR_TWO, "{sample}_good_out_R2.fastq"),
        s1 = os.path.join(PSEQDIR_TWO, "{sample}_single_out_R1.fastq"),
        s2 = os.path.join(PSEQDIR_TWO, "{sample}_single_out_R2.fastq")
    params:
        contigs = os.path.join(RMRD, "final_contigs")
    output:
        os.path.join(RMRD, "{sample}.final_contigs.bam")
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "envs/bowtie.yaml"
    shell:
        """
        bowtie2 --mm -x {params.contigs} -1 {input.r1} -2 {input.r2} \
        -U {input.s1} -U {input.s2} --threads {resources.cpus} | \
        samtools view -bh | samtools sort -o {output} -
        """

rule bai_final_bams:
    input:
        os.path.join(RMRD, "{sample}.final_contigs.bam")
    output:
        os.path.join(RMRD, "{sample}.final_contigs.bam.bai")
    conda:
        "envs/bowtie.yaml"
    shell:
        "samtools index {input}"

rule count_mapped_reads:
    input:
        bam = os.path.join(RMRD, "{sample}.final_contigs.bam"),
        bai = os.path.join(RMRD, "{sample}.final_contigs.bam.bai")
    output:
        os.path.join(RMRD, "{sample}_contig_hits.tsv")
    resources:
        cpus=8
    conda:
        "envs/bowtie.yaml"
    shell:
        """
        samtools idxstats -@ {resources.cpus} {input.bam} | cut -f 1,3 > {output}
        """
rule make_table:
    input:
        expand(os.path.join(RMRD, "{smpl}_contig_hits.tsv"), smpl=SAMPLES)
    output:
        os.path.join(STATS, "sample_coverage.tsv")
    resources:
        mem_mb=64000
    params:
        sct = os.path.join(PMSDIR, "scripts/joinlists.pl")
    shell:
        """
        perl {params.sct} -t {input} > {output}
        """

rule count_contig_lengths:
    input:
        os.path.join(CCMO, "assembly.fasta")
    output:
        os.path.join(STATS, "sequence_lengths.tsv")
    params:
        sct = os.path.join(PMSDIR, "scripts/countfasta.py")
    shell:
        """
        python3 {params.sct} -f {input} -ln > {output}
        """

