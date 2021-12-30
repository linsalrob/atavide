"""

Snakemake rule to run several read-based annotation systems
including kraken, focus, super-focus 

For focus/superfocus we need a directory so we make hard links


"""


import os
import sys



rule run_focus:
    """
    Run focus on the directory with just the prinseq good R1 data
    Run focus on the directory with all the prinseq good data.
    This will make a single focus output file with all data.
    If we want, we can separate those by samples
    """
    input:
        r1 = os.path.join(PSEQDIR_TWO, "{sample}_good_out_R1.fastq"),
        r2 = os.path.join(PSEQDIR_TWO, "{sample}_good_out_R2.fastq")
    output:
        d = directory(os.path.join(RBADIR, "{sample}", "focus")),
        a = temporary(os.path.join(RBADIR, "{sample}", "focus", "output_All_levels.csv")),
        f = temporary(os.path.join(RBADIR, "{sample}", "focus", "output_Family_tabular.csv")),
        k = temporary(os.path.join(RBADIR, "{sample}", "focus", "output_Kingdom_tabular.csv")),
        p = temporary(os.path.join(RBADIR, "{sample}", "focus", "output_Phylum_tabular.csv")),
        s = temporary(os.path.join(RBADIR, "{sample}", "focus", "output_Strain_tabular.csv")),
        c = temporary(os.path.join(RBADIR, "{sample}", "focus", "output_Class_tabular.csv")),
        g = temporary(os.path.join(RBADIR, "{sample}", "focus", "output_Genus_tabular.csv")),
        o = temporary(os.path.join(RBADIR, "{sample}", "focus", "output_Order_tabular.csv")),
        sp = temporary(os.path.join(RBADIR,"{sample}",  "focus", "output_Species_tabular.csv"))
    threads: 8
    resources:
        mem_mb=16000,
        load_superfocus=25
    conda:
        "../envs/focus.yaml"
    shell:
        """
        mkdir -p {output.d};
        focus -q {input.r1} -q {input.r2} -o {output.d} -t {threads}
        """



rule zip_compress_focus:
    input:
        os.path.join(RBADIR, "{sample}", "focus", "output_All_levels.csv"),
        os.path.join(RBADIR, "{sample}", "focus", "output_Family_tabular.csv"),
        os.path.join(RBADIR, "{sample}", "focus", "output_Kingdom_tabular.csv"),
        os.path.join(RBADIR, "{sample}", "focus", "output_Phylum_tabular.csv"),
        os.path.join(RBADIR, "{sample}", "focus", "output_Strain_tabular.csv"),
        os.path.join(RBADIR, "{sample}", "focus", "output_Class_tabular.csv"),
        os.path.join(RBADIR, "{sample}", "focus", "output_Genus_tabular.csv"),
        os.path.join(RBADIR, "{sample}", "focus", "output_Order_tabular.csv"),
        os.path.join(RBADIR, "{sample}", "focus", "output_Species_tabular.csv")
    output:
        os.path.join(RBADIR, "{sample}", "focus", "output_All_levels.csv.zip"),
        os.path.join(RBADIR, "{sample}", "focus", "output_Family_tabular.csv.zip"),
        os.path.join(RBADIR, "{sample}", "focus", "output_Kingdom_tabular.csv.zip"),
        os.path.join(RBADIR, "{sample}", "focus", "output_Phylum_tabular.csv.zip"),
        os.path.join(RBADIR, "{sample}", "focus", "output_Strain_tabular.csv.zip"),
        os.path.join(RBADIR, "{sample}", "focus", "output_Class_tabular.csv.zip"),
        os.path.join(RBADIR, "{sample}", "focus", "output_Genus_tabular.csv.zip"),
        os.path.join(RBADIR, "{sample}", "focus", "output_Order_tabular.csv.zip"),
        os.path.join(RBADIR, "{sample}", "focus", "output_Species_tabular.csv.zip")
    shell:
        """
        for F in {input}; do zip $F.zip $F; done
        """




