"""

Snakemake rule to run several read-based annotation systems
including kraken, focus, super-focus and singlem.

For focus/superfocus we need a directory so we make hard links

For kraken and singlem

"""


import os
import sys

# we will skip the superfocus taxonomy if this is not defined
taxonomy_database = None

# This is the bit Mike hates, but hey, at least I provide an option for you to set the directory
if "TAXONOMY_DB" in os.environ and os.path.exists(os.environ['TAXONOMY_DB']):
    taxonomy_database = os.environ['TAXONOMY_DB']
elif os.path.exists("/home/edwa0468/ncbi/taxonomy/taxonomy.sqlite3"):
     taxonomy_database = "/home/edwa0468/ncbi/taxonomy/taxonomy.sqlite3"
else:
    sys.stderr.write("WARNING: We can not find a taxonomy database, therefore we will skip the superfocus taxonomy step\n")
    sys.stderr.write("You can probably download this database or use Rob's ... ask him!\n")

if taxonomy_database:
    include: "superfocus_taxonomy.snakefile"
else:
    include: "skip_superfocus_taxonomy.snakefile"


if 'SUPERFOCUS_DB' not in os.environ:
    sys.stderr.write("FATAL: Please set the location of your superfocus databases using the SUPERFOCUS_DB environment variable\n")
    sys.exit(1)



# just check there is something to actually do!
if len(SAMPLES) == 0:
    sys.stderr.write("FATAL: We could not detect any samples at all.\n")
    sys.stderr.write(f"Do you have a directory called {PSEQDIR_TWO} with some fastq files in it?\n")
    sys.exit()


rule read_annotation_all:
    """
    This rule is probably not used unless you use this as a standalone snakemake script!
    """
    input:
        expand(
            [
                os.path.join(RBADIR, "{sample}", "focus", "output_All_levels.csv"),
                os.path.join(RBADIR, "{sample}", "superfocus", "output_all_levels_and_function.xls"),
                os.path.join(RBADIR, "{sample}", "kraken", "{sample}.report.tsv"),
                os.path.join(RBADIR, "{sample}", "singlem", "singlem_otu_table.tsv"),
                os.path.join(RBADIR, "{sample}", "kraken", "{sample}.output.tsv"),
            ],
               sample=SAMPLES),
        os.path.join(RBADIR, "superfocus_taxonomy.tsv.zip")



# focus/super-focus uses directories, and so if we point it at the whole prinseq
# directory it will process all the data in there.

# therefore, we make a new directory and link just the file we need. Doesn't consume
# any more space!

rule link_psq_good:
    input:
        r1 = os.path.join(PSEQDIR_TWO, "{sample}_good_out_R1.fastq")
    output:
        d = directory(os.path.join(RBADIR, "{sample}", "prinseq_good")),
        f = os.path.join(RBADIR, "{sample}", "prinseq_good", "{sample}.good_out_R1.fastq")
    params:
        d = os.path.join(RBADIR, "{sample}", "prinseq_good")
    shell:
        """
        mkdir -p {params.d} && ln {input.r1} {output.f}
        """

rule run_focus:
    """
    Run focus on the directory with just the prinseq good R1 data
    """
    input:
        os.path.join(RBADIR, "{sample}", "prinseq_good")
    output:
        d = directory(os.path.join(RBADIR, "{sample}", "focus")),
        a = os.path.join(RBADIR, "{sample}", "focus", "output_All_levels.csv"),
        f = temporary(os.path.join(RBADIR, "{sample}", "focus", "output_Family_tabular.csv")),
        k = temporary(os.path.join(RBADIR, "{sample}", "focus", "output_Kingdom_tabular.csv")),
        p = temporary(os.path.join(RBADIR, "{sample}", "focus", "output_Phylum_tabular.csv")),
        s = temporary(os.path.join(RBADIR, "{sample}", "focus", "output_Strain_tabular.csv")),
        c = temporary(os.path.join(RBADIR, "{sample}", "focus", "output_Class_tabular.csv")),
        g = temporary(os.path.join(RBADIR, "{sample}", "focus", "output_Genus_tabular.csv")),
        o = temporary(os.path.join(RBADIR, "{sample}", "focus", "output_Order_tabular.csv")),
        sp = temporary(os.path.join(RBADIR, "{sample}", "focus", "output_Species_tabular.csv"))
    resources:
        cpus=8,
        mem_mb=16000
    conda:
        "../envs/focus.yaml"
    shell:
        """
        focus -q {input} -o {output.d} -t {resources.cpus}
        """


rule run_superfocus:
    input:
        os.path.join(RBADIR, "{sample}", "prinseq_good")
    output:
        d = directory(os.path.join(RBADIR, "{sample}", "superfocus")),
        a = os.path.join(RBADIR, "{sample}", "superfocus", "output_all_levels_and_function.xls"),
        m8 = os.path.join(RBADIR, "{sample}", "superfocus", "{sample}.good_out_R1.fastq_alignments.m8"),
        l1 = temporary(os.path.join(RBADIR, "{sample}", "superfocus", "output_subsystem_level_1.xls")),
        l2 = temporary(os.path.join(RBADIR, "{sample}", "superfocus", "output_subsystem_level_2.xls")),
        l3 = temporary(os.path.join(RBADIR, "{sample}", "superfocus", "output_subsystem_level_3.xls")),
    resources:
        cpus=16,
        mem_mb=32000,
        time=7200
    conda:
        "../envs/superfocus.yaml"
    shell:
        """
        superfocus -q {input} -dir {output.d} -a diamond -t {resources.cpus} -n 0
        """

rule merge_sf_outputs:
    input:
        expand(os.path.join(RBADIR, "{smps}", "superfocus", "output_all_levels_and_function.xls"), smps=SAMPLES)
    output:
        os.path.join(RBADIR, "superfocus_functions.tsv")
    shell:
        """
        perl /home/edwa0468/GitHubs/atavide/workflow/scripts/joinsuperfocus.pl {input}  > {output}
        """


rule run_kraken:
    input:
        r1 = os.path.join(PSEQDIR_TWO, "{sample}_good_out_R1.fastq.gz"),
        r2 = os.path.join(PSEQDIR_TWO, "{sample}_good_out_R2.fastq.gz")
    output:
        rt = os.path.join(RBADIR, "{sample}", "kraken", "{sample}.report.tsv"),
        ot = os.path.join(RBADIR, "{sample}", "kraken", "{sample}.output.tsv")
    resources:
        cpus=8,
        mem_mb=400000
    conda:
        "../envs/kraken.yaml"
    shell:
        """
        kraken2 --report {output.rt} \
                --output {output.ot} \
                --threads {resources.cpus} \
                {input.r1}
        """

rule run_singlem:
    input:
        r1 = os.path.join(PSEQDIR_TWO, "{sample}_good_out_R1.fastq.gz"),
        r2 = os.path.join(PSEQDIR_TWO, "{sample}_good_out_R2.fastq.gz")
    output:
        d = directory(os.path.join(RBADIR, "{sample}", "singlem")),
        otu = os.path.join(RBADIR, "{sample}", "singlem", "singlem_otu_table.tsv")
    conda:
        "../envs/singlem.yaml"
    resources:
        cpus=8,
        mem_mb=40000
    shell:
        """
        mkdir --parents {output.d};
        singlem pipe --forward {input.r1} --reverse {input.r2} --otu_table {output.otu} --output_extras --threads {resources.cpus}
        """

