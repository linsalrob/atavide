



rule run_kraken:
    input:
        r1 = os.path.join(PSEQDIR_TWO, "{sample}_good_out_R1.fastq"),
        r2 = os.path.join(PSEQDIR_TWO, "{sample}_good_out_R2.fastq")
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



rule kraken_taxonomy:
    input:
        os.path.join(RBADIR, "{sample}", "kraken", "{sample}.output.tsv")
    output:
        os.path.join(RBADIR, "{sample}", "kraken", "{sample}.taxonomy.tsv")
    resources:
        cpus=4,
        mem_mb=40000
    params:
        t = TAXON
    conda:
        "../envs/taxonkit.yaml"
    shell:
        """
        grep -v ^U {input} | \
        taxonkit lineage -i 3 --data-dir {params.t} | \
        taxonkit reformat --data-dir {params.t} -i 6 -f \
        "Root; d__{{k}}; p__{{p}}; c__{{c}}; o__{{o}}; f__{{f}}; g__{{g}}" \
        --fill-miss-rank | cut -f 2,3,7 taxonout2 > {output}
        """
