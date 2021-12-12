
# Use seqtk to summarize the fastq quality scores. We make a table that has reads
# and the average quality score at each position.



rule fqchk_r1_sequences:
    input:
        r1 = os.path.join(PSEQDIR_TWO, "{sample}_good_out_R1.fastq")
    output:
        os.path.join(RBADIR, "{sample}", "statistics", "{sample}_good_out_R1.fqchk.tsv")
    threads: 8
    resources:
        mem_mb=25000
    params:
        d =  os.path.join(RBADIR, "{sample}", "statistics")
    conda:
        "../envs/seqtk.yaml"
    shell:
        """
        mkdir -p {params.d};
        seqtk fqchk {input.r1} > {output}
        """


rule get_fqchk_cols:
    input:
        os.path.join(RBADIR, "{sample}", "statistics", "{sample}_good_out_R1.fqchk.tsv")
    output:
        temporary(os.path.join(RBADIR, "{sample}", "statistics", "{sample}_good_out_R1.fqchk.tsv.tmp"))
    params:
        s = "{sample}"
    shell:
        """
        cut -f 1,8 {input} | tail -n +2 | grep -v ALL | sed -e 's/avgQ/{params.s}/' > {output}
        """


rule join_fqchk_cols:
    input:
        expand(os.path.join(RBADIR, "{smp}", "statistics", "{smp}_good_out_R1.fqchk.tsv.tmp"), smp=SAMPLES)
    output:
        os.path.join(STATS, "av_quality_scores_by_position.tsv")
    params:
        sct = os.path.join(ATAVIDE_DIR, "scripts/joinlists.pl")
    shell:
        """
        perl {params.sct} -h {input} > {output}
        """
