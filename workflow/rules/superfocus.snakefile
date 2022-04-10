"""

snakemake rule just for superfocus

"""


import os
import sys


rule run_superfocus:
    # This is the superfocus command for diamond
    #    superfocus -q {input.r1} -q {input.r2} -dir {output.d} -a diamond -t {threads} -n 0 -tmp $(mktemp -d -p {params.TMPDIR})
    # but be aware that diamond does not play well run in parallel like snakemake uses.
    input:
        r1 = os.path.join(PSEQDIR_TWO, "{sample}_good_out_R1.fastq"),
        r2 = os.path.join(PSEQDIR_TWO, "{sample}_good_out_R2.fastq")
    output:
        d = directory(os.path.join(RBADIR, "{sample}", "superfocus")),
        a = temporary(os.path.join(RBADIR, "{sample}", "superfocus", "output_all_levels_and_function.xls")),
        l1 = temporary(os.path.join(RBADIR, "{sample}", "superfocus", "output_subsystem_level_1.xls")),
        l2 = temporary(os.path.join(RBADIR, "{sample}", "superfocus", "output_subsystem_level_2.xls")),
        l3 = temporary(os.path.join(RBADIR, "{sample}", "superfocus", "output_subsystem_level_3.xls")),
        m81 = os.path.join(RBADIR, "{sample}", "superfocus", "{sample}_good_out_R1.fastq_alignments.m8"),
        m82 = os.path.join(RBADIR, "{sample}", "superfocus", "{sample}_good_out_R2.fastq_alignments.m8"),
    threads: 16
    params: 
        TMPDIR=TMPDIR
    resources:
        mem_mb=64000,
        load_superfocus=25,
        time=7200
    conda:
        "../envs/superfocus.yaml"
    shell:
        """
        superfocus -q {input.r1} -q {input.r2} -dir {output.d} -a mmseqs2 -t {threads} -n 0 -tmp $(mktemp -d -p {params.TMPDIR})
        """

rule merge_sf_outputs:
    input:
        expand(os.path.join(RBADIR, "{smps}", "superfocus", "output_all_levels_and_function.xls"), smps=SAMPLES)
    output:
        os.path.join(RBADIR, "superfocus_functions.tsv.gz")
    params:
        sct = os.path.join(ATAVIDE_DIR, "scripts/joinsuperfocus.pl"),
        out = os.path.join(RBADIR, "superfocus_functions.tsv")
    shell:
        """
        perl {params.sct} {input}  > {params.out} && gzip {params.out}
        """




"""
In the rules below we use the m8 files from superfocus
to crate taxonomy tables.
"""

rule superfocus_taxonomy:
    input:
        m8 = os.path.join(RBADIR, "{sample}", "superfocus", "{sample}_good_out_R1.fastq_alignments.m8"),
    output:
        temporary(os.path.join(RBADIR, "{sample}", "superfocus", "{sample}_good_out_R1.taxonomy"))
    params:
        t = TAXON
    threads: 4
    resources:
        mem_mb=16000,
        load_superfocus=25
    conda:
        "../envs/taxonkit.yaml"
    shell:
        """
        perl -F"\\t" -lane 'if ($F[1] =~ /fig\|(\d+)\.\d+/ && $1 != 6666666) \
           {{print "$F[0]\\t$1"}}' {input.m8} | \
        taxonkit lineage -j {threads} -i 2 --data-dir {params.t} | \
        taxonkit reformat -j {threads} --data-dir {params.t} -i 3 \
          -f "Root; d__{{k}}; p__{{p}}; c__{{c}}; o__{{o}}; f__{{f}}; g__{{g}}" \
          --fill-miss-rank | \
        cut -f 1,2,4 > {output}
        """




rule zip_compress_superfocus:
    input:
        os.path.join(RBADIR, "{sample}", "superfocus", "output_all_levels_and_function.xls"),
        os.path.join(RBADIR, "{sample}", "superfocus", "output_subsystem_level_1.xls"),
        os.path.join(RBADIR, "{sample}", "superfocus", "output_subsystem_level_2.xls"),
        os.path.join(RBADIR, "{sample}", "superfocus", "output_subsystem_level_3.xls"),
        os.path.join(RBADIR, "{sample}", "superfocus", "{sample}_good_out_R1.taxonomy"),
    output:
        os.path.join(RBADIR, "{sample}", "superfocus", "output_all_levels_and_function.xls.zip"),
        os.path.join(RBADIR, "{sample}", "superfocus", "output_subsystem_level_1.xls.zip"),
        os.path.join(RBADIR, "{sample}", "superfocus", "output_subsystem_level_2.xls.zip"),
        os.path.join(RBADIR, "{sample}", "superfocus", "output_subsystem_level_3.xls.zip"),
        os.path.join(RBADIR, "{sample}", "superfocus", "{sample}_good_out_R1.taxonomy.zip"),
    shell:
        """
        for F in {input}; do zip $F.zip $F; done
        """



