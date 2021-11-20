



rule run_kraken:
    input:
        r1 = os.path.join(PSEQDIR_TWO, "{sample}_good_out_R1.fastq"),
        r2 = os.path.join(PSEQDIR_TWO, "{sample}_good_out_R2.fastq")
    output:
        rt = os.path.join(RBADIR, "{sample}", "kraken", "{sample}.report.tsv"),
        ot = os.path.join(RBADIR, "{sample}", "kraken", "{sample}.output.tsv")
    threads: 8
    resources:
        mem_mb=250000
    conda:
        "../envs/kraken.yaml"
    shell:
        """
        kraken2 --report {output.rt} \
                --output {output.ot} \
                --threads {threads} \
                {input.r1}
        """



rule kraken_taxonomy:
    input:
        os.path.join(RBADIR, "{sample}", "kraken", "{sample}.output.tsv")
    output:
        os.path.join(RBADIR, "{sample}", "kraken", "{sample}.taxonomy.tsv")
    threads: 4
    resources:
        mem_mb=40000
    params:
        t = TAXON
    conda:
        "../envs/taxonkit.yaml"
    shell:
        """
        grep -v ^U {input} | \
        taxonkit lineage -j {threads} -i 3 --data-dir {params.t} | \
        taxonkit reformat -j {threads} --data-dir {params.t} -i 6 -f \
        "Root; d__{{k}}; p__{{p}}; c__{{c}}; o__{{o}}; f__{{f}}; g__{{g}}" \
        --fill-miss-rank | cut -f 2,3,7 > {output}
        """
