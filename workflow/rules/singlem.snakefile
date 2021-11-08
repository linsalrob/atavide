
rule run_singlem:
    input:
        r1 = os.path.join(PSEQDIR_TWO, "{sample}_good_out_R1.fastq"),
        r2 = os.path.join(PSEQDIR_TWO, "{sample}_good_out_R2.fastq")
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


