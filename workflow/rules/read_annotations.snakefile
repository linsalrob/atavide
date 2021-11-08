"""

Snakemake rule to run several read-based annotation systems
including kraken, focus, super-focus 

For focus/superfocus we need a directory so we make hard links


"""


import os
import sys

include: "superfocus_taxonomy.snakefile"


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
    params:
        sct = os.path.join(ATAVIDE_DIR, "scripts/joinsuperfocus.pl")
    shell:
        """
        perl {params.sct} {input}  > {output}
        """

