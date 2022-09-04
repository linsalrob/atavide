

rule merge_read_annotations:
    input:
        bam = os.path.join(RMRD, "{sample}." + SAMPLE_ID + ".assembly.bam"),
        bai = os.path.join(RMRD, "{sample}." + SAMPLE_ID + ".assembly.bam.bai"),
        singlem = os.path.join(RBADIR, "{sample}", "singlem", "singlem_otu_table.tsv"),
        kraken = os.path.join(RBADIR, "{sample}", "kraken", "{sample}.taxonomy.tsv"),
    output:
        temporary(os.path.join(RBADIR, "{sample}", "{sample}_contig_taxonomy.comparison.tsv")),
        temporary(os.path.join(RBADIR, "{sample}", "{sample}_contig_taxonomy.kraken.allmatches.tsv")),
        temporary(os.path.join(RBADIR, "{sample}", "{sample}_contig_taxonomy.kraken.besthits.tsv")),
        temporary(os.path.join(RBADIR, "{sample}", "{sample}_contig_taxonomy.singlem.allmatches.tsv")),
        temporary(os.path.join(RBADIR, "{sample}", "{sample}_contig_taxonomy.singlem.besthits.tsv"))
    conda:
        "../envs/pysam.yaml"
    resources:
        mem_mb=8000
    params:
        out = os.path.join(RBADIR, "{sample}", "{sample}_contig_taxonomy"),
        sct = os.path.join(ATAVIDE_DIR, "scripts/map_reads_to_contigs.py")
    shell:
        """
        python {params.sct} --bam {input.bam} \
            --singlem {input.singlem} \
            --kraken {input.kraken} \
            --output {params.out} \
            --verbose
        """
