#####################################################################
#                                                                   #
# Compress all the outputs                                          #
#                                                                   #
# For files that we are likely to use on a linux machine            #
# (e.g. *fastq, *fasta) we use gzip to compress the output.         # 
#                                                                   #
# For .tsv, .txt, and .xls we use zip to compress the outputs       #
# so they are easier to expand on a desktop                         #
#                                                                   #
#####################################################################



rule compress_kraken:
    input:
        os.path.join(RBADIR, "{sample}", "kraken", "{sample}.report.tsv"),
        os.path.join(RBADIR, "{sample}", "kraken", "{sample}.output.tsv"),
        os.path.join(RBADIR, "{sample}", "kraken", "{sample}.taxonomy.tsv"),
        os.path.join(RBADIR, "{sample}", "singlem", "singlem_otu_table.tsv"),
        os.path.join(RBADIR, "{sample}", "{sample}_contig_taxonomy.comparison.tsv"),
        os.path.join(RBADIR, "{sample}", "{sample}_contig_taxonomy.kraken.allmatches.tsv"),
        os.path.join(RBADIR, "{sample}", "{sample}_contig_taxonomy.kraken.besthits.tsv"),
        os.path.join(RBADIR, "{sample}", "{sample}_contig_taxonomy.singlem.allmatches.tsv"),
        os.path.join(RBADIR, "{sample}", "{sample}_contig_taxonomy.singlem.besthits.tsv")
    output:
        os.path.join(RBADIR, "{sample}", "kraken", "{sample}.report.tsv.zip"),
        os.path.join(RBADIR, "{sample}", "kraken", "{sample}.output.tsv.zip"),
        os.path.join(RBADIR, "{sample}", "kraken", "{sample}.taxonomy.tsv.zip"),
        os.path.join(RBADIR, "{sample}", "singlem", "singlem_otu_table.tsv.zip"),
        os.path.join(RBADIR, "{sample}", "{sample}_contig_taxonomy.comparison.tsv.zip"),
        os.path.join(RBADIR, "{sample}", "{sample}_contig_taxonomy.kraken.allmatches.tsv.zip"),
        os.path.join(RBADIR, "{sample}", "{sample}_contig_taxonomy.kraken.besthits.tsv.zip"),
        os.path.join(RBADIR, "{sample}", "{sample}_contig_taxonomy.singlem.allmatches.tsv.zip"),
        os.path.join(RBADIR, "{sample}", "{sample}_contig_taxonomy.singlem.besthits.tsv.zip")
    shell:
        """
        for F in {input}; do zip $F.zip $F; done
        """




rule compress_wildcard_text:
    input:
        os.path.join(RMRD, "{sample}_contig_hits.tsv"),
    output:
        os.path.join(RMRD, "{sample}_contig_hits.tsv.zip"),
    shell:
        """
        for F in {input}; do zip $F.zip $F; done
        """

rule compress_text_files:
    input:
        os.path.join(STATS, "final_assembly.txt"),
        os.path.join(STATS, "sample_coverage.tsv"),
        os.path.join(STATS, "sequence_lengths.tsv")
    output:
        os.path.join(STATS, "final_assembly.txt.zip"),
        os.path.join(STATS, "sample_coverage.tsv.zip"),
        os.path.join(STATS, "sequence_lengths.tsv.zip")
    shell:
        """
        for F in {input}; do zip $F.zip $F; done
        """


