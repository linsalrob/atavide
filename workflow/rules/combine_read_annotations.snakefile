

rule merge_read_annotations:
    input:
        bam = os.path.join(RMRD, "{sample}.final_contigs.bam"),
        bai = os.path.join(RMRD, "{sample}.final_contigs.bam.bai"),
        singlem = os.path.join(RBADIR, "{sample}", "singlem", "singlem_otu_table.tsv"),
        kraken = os.path.join(RBADIR, "{sample}", "kraken", "{sample}.taxonomy.tsv")
        superfocus = os.path.join(RBADIR, "{sample}", "superfocus", "{sample}_good_out.taxonomy")
    output:
        os.path.join(RBADIR, "{sample}", "{sample}_contig_taxonomy.comparison.tsv"),
        os.path.join(RBADIR, "{sample}", "{sample}_contig_taxonomy.kraken.allmatches.tsv"),
        os.path.join(RBADIR, "{sample}", "{sample}_contig_taxonomy.kraken.besthits.tsv"),
        os.path.join(RBADIR, "{sample}", "{sample}_contig_taxonomy.singlem.allmatches.tsv"),
        os.path.join(RBADIR, "{sample}", "{sample}_contig_taxonomy.singlem.besthits.tsv")
    conda:
        "../envs/pysam.yaml"
    params:
        out = os.path.join(RBADIR, "{sample}", "{sample}_contig_taxonomy"),
        sct = os.path.join(ATAVIDE_DIR, "scripts/map_reads_to_contigs.py")
    shell:
        """
        python3 {params.sct} --bam {input.bam} \
            --singlem {input.singlem} \
            --kraken {input.kraken} \
            --superfocus {input.superfocus} \
            --output {params.out} \
            --verbose
        """
