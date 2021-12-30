"""

Snakemake rule to run several read-based annotation systems
including kraken, focus, super-focus 

For focus/superfocus we need a directory so we make hard links


"""


import os
import sys

include: "superfocus_taxonomy.snakefile"



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
                os.path.join(RBADIR, "focus", "output_All_levels.csv"),
                os.path.join(RBADIR, "superfocus", "output_all_levels_and_function.xls"),
            ],
               sample=SAMPLES),




rule run_focus:
    """
    Run focus on the directory with just the prinseq good R1 data
    Run focus on the directory with all the prinseq good data.
    This will make a single focus output file with all data.
    If we want, we can separate those by samples
    """
    input:
        os.path.join(PSEQDIR_TWO)
    output:
        d = directory(os.path.join(RBADIR, "focus")),
        a = temporary(os.path.join(RBADIR, "focus", "output_All_levels.csv")),
        f = temporary(os.path.join(RBADIR, "focus", "output_Family_tabular.csv")),
        k = temporary(os.path.join(RBADIR, "focus", "output_Kingdom_tabular.csv")),
        p = temporary(os.path.join(RBADIR, "focus", "output_Phylum_tabular.csv")),
        s = temporary(os.path.join(RBADIR, "focus", "output_Strain_tabular.csv")),
        c = temporary(os.path.join(RBADIR, "focus", "output_Class_tabular.csv")),
        g = temporary(os.path.join(RBADIR, "focus", "output_Genus_tabular.csv")),
        o = temporary(os.path.join(RBADIR, "focus", "output_Order_tabular.csv")),
        sp = temporary(os.path.join(RBADIR, "focus", "output_Species_tabular.csv"))
    threads: 8
    resources:
        mem_mb=16000,
        load_superfocus=25
    conda:
        "../envs/focus.yaml"
    shell:
        """
        focus -q {input} -o {output.d} -t {threads}
        """


rule run_superfocus:
    input:
        os.path.join(PSEQDIR_TWO)
    output:
        d = directory(os.path.join(RBADIR, "superfocus")),
        a = temporary(os.path.join(RBADIR, "superfocus", "output_all_levels_and_function.xls")),
        l1 = temporary(os.path.join(RBADIR, "superfocus", "output_subsystem_level_1.xls")),
        l2 = temporary(os.path.join(RBADIR, "superfocus", "output_subsystem_level_2.xls")),
        l3 = temporary(os.path.join(RBADIR, "superfocus", "output_subsystem_level_3.xls")),
    threads: 16
    resources:
        mem_mb=32000,
        load_superfocus=25,
        time=7200
    conda:
        "../envs/superfocus.yaml"
    shell:
        """
        superfocus -q {input} -dir {output.d} -a diamond -t {threads} -n 0
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



rule zip_sf_output:
    input:
        os.path.join(RBADIR, "superfocus_functions.tsv"),
    output:
        os.path.join(RBADIR, "superfocus_functions.tsv.zip"),
    shell:
        """
        for F in {input}; do zip $F.zip $F; done
        """

rule zip_compress_focus:
    input:
        os.path.join(RBADIR, "focus", "output_All_levels.csv"),
        os.path.join(RBADIR, "focus", "output_Family_tabular.csv"),
        os.path.join(RBADIR, "focus", "output_Kingdom_tabular.csv"),
        os.path.join(RBADIR, "focus", "output_Phylum_tabular.csv"),
        os.path.join(RBADIR, "focus", "output_Strain_tabular.csv"),
        os.path.join(RBADIR, "focus", "output_Class_tabular.csv"),
        os.path.join(RBADIR, "focus", "output_Genus_tabular.csv"),
        os.path.join(RBADIR, "focus", "output_Order_tabular.csv"),
        os.path.join(RBADIR, "focus", "output_Species_tabular.csv")
    output:
        os.path.join(RBADIR, "focus", "output_All_levels.csv.zip"),
        os.path.join(RBADIR, "focus", "output_Family_tabular.csv.zip"),
        os.path.join(RBADIR, "focus", "output_Kingdom_tabular.csv.zip"),
        os.path.join(RBADIR, "focus", "output_Phylum_tabular.csv.zip"),
        os.path.join(RBADIR, "focus", "output_Strain_tabular.csv.zip"),
        os.path.join(RBADIR, "focus", "output_Class_tabular.csv.zip"),
        os.path.join(RBADIR, "focus", "output_Genus_tabular.csv.zip"),
        os.path.join(RBADIR, "focus", "output_Order_tabular.csv.zip"),
        os.path.join(RBADIR, "focus", "output_Species_tabular.csv.zip")
    shell:
        """
        for F in {input}; do zip $F.zip $F; done
        """



rule zip_compress_superfocus:
    input:
        os.path.join(RBADIR, "superfocus", "output_all_levels_and_function.xls"),
        os.path.join(RBADIR, "superfocus", "output_subsystem_level_1.xls"),
        os.path.join(RBADIR, "superfocus", "output_subsystem_level_2.xls"),
        os.path.join(RBADIR, "superfocus", "output_subsystem_level_3.xls"),
    output:
        os.path.join(RBADIR, "superfocus", "output_all_levels_and_function.xls.zip"),
        os.path.join(RBADIR, "superfocus", "output_subsystem_level_1.xls.zip"),
        os.path.join(RBADIR, "superfocus", "output_subsystem_level_2.xls.zip"),
        os.path.join(RBADIR, "superfocus", "output_subsystem_level_3.xls.zip"),
    shell:
        """
        for F in {input}; do zip $F.zip $F; done
        """





