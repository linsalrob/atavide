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


rule compress_prinseq_output:
    """
    Note: we do this indirectly, with a stdout redirect because
    elsewhere we set a link to this file, and gzip
    breaks if there is >1 link to the file.

    We have also set these files to be temporary, so once 
    we have compressed them, they should be deleted!

    Same thing as gzip, but in two steps
    """
    input:
        os.path.join(PSEQDIR, "{sample}_good_out_R1.fastq"),
        os.path.join(PSEQDIR, "{sample}_good_out_R2.fastq"),
        os.path.join(PSEQDIR, "{sample}_single_out_R1.fastq"),
        os.path.join(PSEQDIR, "{sample}_single_out_R2.fastq")
    output:
        os.path.join(PSEQDIR, "{sample}_good_out_R1.fastq.gz"),
        os.path.join(PSEQDIR, "{sample}_good_out_R2.fastq.gz"),
        os.path.join(PSEQDIR, "{sample}_single_out_R1.fastq.gz"),
        os.path.join(PSEQDIR, "{sample}_single_out_R2.fastq.gz")
    shell:
        """
        for F in {input}; do gzip -c $F > $F.gz; done
        """


rule compress_unassembled_reads:
    input:
        os.path.join(UNASSM, "{sample}.unassembled.R1.fastq"),
        os.path.join(UNASSM, "{sample}.unassembled.R2.fastq"),
        os.path.join(UNASSM, "{sample}.unassembled.singles.fastq")
    output:
        os.path.join(UNASSM, "{sample}.unassembled.R1.fastq.gz"),
        os.path.join(UNASSM, "{sample}.unassembled.R2.fastq.gz"),
        os.path.join(UNASSM, "{sample}.unassembled.singles.fastq.gz")
    shell:
        """
        for F in {input}; do gzip $F; done
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

rule zip_sf_output:
    input:
        os.path.join(RBADIR, "superfocus_functions.tsv"),
    output:
        os.path.join(RBADIR, "superfocus_functions.tsv.zip"),
    shell:
        """
        for F in {input}; do zip $F.zip $F; done
        """

rule zip_compress_read_annotation_outputs:
    input:
        os.path.join(RBADIR, "{sample}", "focus", "output_All_levels.csv"),
        os.path.join(RBADIR, "{sample}", "superfocus", "output_all_levels_and_function.xls"),
        os.path.join(RBADIR, "{sample}", "kraken", "{sample}.report.tsv"),
        os.path.join(RBADIR, "{sample}", "kraken", "{sample}.output.tsv"),
        otu = os.path.join(RBADIR, "{sample}", "singlem", "singlem_otu_table.tsv")
    output:
        os.path.join(RBADIR, "{sample}", "focus", "output_All_levels.csv.zip"),
        os.path.join(RBADIR, "{sample}", "superfocus", "output_all_levels_and_function.xls.zip"),
        os.path.join(RBADIR, "{sample}", "kraken", "{sample}.report.tsv.zip"),
        os.path.join(RBADIR, "{sample}", "kraken", "{sample}.output.tsv.zip"),
        otu = os.path.join(RBADIR, "{sample}", "singlem", "singlem_otu_table.tsv.zip")
    shell:
        """
        for F in {input}; do zip $F.zip $F; done
        """

rule gzip_compress_read_annotation_outputs:
    input:
       os.path.join(RBADIR, "{sample}", "superfocus", "{sample}.good_out_R1.fastq_alignments.m8"),
    output:
       os.path.join(RBADIR, "{sample}", "superfocus", "{sample}.good_out_R1.fastq_alignments.m8.gz"),
    shell:
        """
        for F in {input}; do gzip $F; done
        """




