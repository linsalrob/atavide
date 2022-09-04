"""

Rules for quality control and quality assurance
"""

rule prinseq:
    input:
        r1 = os.path.join(READDIR, PATTERN_R1),
        r2 = os.path.join(READDIR, PATTERN_R2)
    output:
        r1 = temporary(os.path.join(PSEQDIR, "{sample}_good_out_R1.fastq")),
        r2 = temporary(os.path.join(PSEQDIR, "{sample}_good_out_R2.fastq")),
        s1 = temporary(os.path.join(PSEQDIR, "{sample}_single_out_R1.fastq")),
        s2 = temporary(os.path.join(PSEQDIR, "{sample}_single_out_R2.fastq")),
    conda: "../envs/prinseq.yaml"
    params:
        o = os.path.join(PSEQDIR, "{sample}")
    shell:
        """
            prinseq++ -min_len 60 -min_qual_mean 25 -ns_max_n 1 -derep 1 \
                    -out_format 0 -trim_tail_left 5 -trim_tail_right 5 \
                    -ns_max_n 5  -trim_qual_type min -trim_qual_left 30 \
                    -trim_qual_right 30 -trim_qual_window 10 \
                    -threads {threads} \
                    -out_name {params.o} \
                    -out_bad /dev/null \
                    -out_bad2 /dev/null \
                    -fastq {input.r1} \
                    -fastq2 {input.r2};
        """

rule compress_prinseq:
    """
    Here we compress the prinseq files, but we also have as input
    the output of some of the other rules to make sure we don't 
    compress too early!
    """
    input:
        os.path.join(PSEQDIR, "{sample}_good_out_R1.fastq"),
        os.path.join(PSEQDIR, "{sample}_good_out_R2.fastq"),
        os.path.join(PSEQDIR, "{sample}_single_out_R1.fastq"),
        os.path.join(PSEQDIR, "{sample}_single_out_R2.fastq"),
        os.path.join(RMRD, "{sample}_rpkm.tsv"),
        os.path.join(RBADIR, "{sample}", "focus", "output_All_levels.csv.zip"),
        os.path.join(STATS, "kraken_species_rarefaction.tsv"),
        os.path.join(RBADIR, "{sample}", "superfocus", "output_all_levels_and_function.xls.zip"),
        os.path.join(RBADIR, "{sample}", "singlem", "singlem_otu_table.tsv"),
        os.path.join(STATS, "av_quality_scores_by_position.tsv")
    output:
        r1 = os.path.join(PSEQDIR, "{sample}_good_out_R1.fastq.gz"),
        r2 = os.path.join(PSEQDIR, "{sample}_good_out_R2.fastq.gz"),
        s1 = os.path.join(PSEQDIR, "{sample}_single_out_R1.fastq.gz"),
        s2 = os.path.join(PSEQDIR, "{sample}_single_out_R2.fastq.gz"),
    shell:
        """
        for F in {input}; do gzip -c $F > $F.gz; done
        """



