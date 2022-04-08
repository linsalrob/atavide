
rule run_singlem:
    input:
        r1 = os.path.join(PSEQDIR_TWO, FASTQ)
    output:
        d = directory(os.path.join(RBADIR, "{sample}", "singlem")),
        otu = os.path.join(RBADIR, "{sample}", "singlem", "singlem_otu_table.tsv")
    conda:
        "../envs/singlem.yaml"
    threads: 8
    resources:
        mem_mb=40000
    shell:
        """
        mkdir --parents {output.d};
        singlem pipe --forward {input.r1} --otu_table {output.otu} --output_extras --threads {threads}
        """


