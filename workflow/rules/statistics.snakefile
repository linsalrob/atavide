"""
The rules that generate statistics output
"""


rule input_read_stats:
    """
    Count the statistics for the initial sequences
    """
    input:
        fqdir = READDIR
    output:
        stats = os.path.join(STATS, "initial_read_statistics.tsv")
    script:
        "../scripts/countfastq.py"

rule after_qc_stats:
    """
    Count the statistics after complete QC
    """
    input:
        fqdir = PSEQDIR_TWO
    output:
        stats = os.path.join(STATS, "post_qc_statistics.tsv")
    script:
        "../scripts/countfastq.py"

rule final_assembly_stats:
    input:
        os.path.join(config['directories']['assemblies'], f"{SAMPLE_ID}_assembly.fasta")
    output:
        os.path.join(STATS, "final_assembly.txt")
    params:
        sct = os.path.join(ATAVIDE_DIR, "scripts/countfasta.py")
    shell:
        """
        python3 {params.sct} -f {input} > {output}
        """

rule make_table:
    input:
        expand(os.path.join(RMRD, "{smpl}_contig_hits.tsv"), smpl=SAMPLES)
    output:
        os.path.join(STATS, "sample_coverage.tsv")
    resources:
        mem_mb=64000
    params:
        sct = os.path.join(ATAVIDE_DIR, "scripts/joinlists.pl")
    shell:
        """
        perl {params.sct} -t {input} > {output}
        """

rule rpkm:
    input:
        expand(os.path.join(RMRD, "{smpl}_rpkm.tsv"), smpl=SAMPLES)
    output:
         os.path.join(STATS, "sample_rpkm.tsv")
    resources:
        mem_mb=64000
    params:
        sct = os.path.join(ATAVIDE_DIR, "scripts/joinlists.pl")
    shell:
        """
        perl {params.sct} -t {input} > {output}
        """


rule count_contig_lengths:
    input:
        os.path.join(config['directories']['assemblies'], f"{SAMPLE_ID}_assembly.fasta")
    output:
        os.path.join(STATS, "sequence_lengths.tsv")
    params:
        sct = os.path.join(ATAVIDE_DIR, "scripts/countfasta.py")
    shell:
        """
        python3 {params.sct} -f {input} -ln > {output}
        """

rule sf_taxonomy:
    input:
        expand(os.path.join(RBADIR, "{smpl}", "superfocus", "{smpl}_R1_taxonomy_best_hits.tsv"), smpl=SAMPLES)
    output:
        os.path.join(STATS, "superfocus_best_hits_taxonomy.tsv")
    resources:
        mem_mb=64000
    params:
        sct = os.path.join(ATAVIDE_DIR, "scripts/joinlists.pl")
    shell:
        """
        perl {params.sct} -t {input} > {output}
        """


rule merge_sf_outputs:
    input:
        expand(os.path.join(RBADIR, "{smps}", "superfocus", "output_all_levels_and_function.xls"), smps=SAMPLES)
    output:
        os.path.join(STATS, "superfocus_functions.tsv.gz")
    params:
        sct = os.path.join(ATAVIDE_DIR, "scripts/joinsuperfocus.pl"),
        out = os.path.join(STATS, "superfocus_functions.tsv")
    shell:
        """
        perl {params.sct} {input}  > {params.out} && gzip {params.out}
        """






