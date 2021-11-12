
rule make_h5_table:
    input:
        os.path.join(STATS, "sample_coverage.tsv")
    output:
        h5 = os.path.join(STATS, "sample_coverage.h5"),
        idx = os.path.join(STATS, "sample_coverage.contig_index.tsv")
    resources:
        mem_mb=64000
    conda:
        "../envs/h5py.yaml"
    params:
        sct = os.path.join(ATAVIDE_DIR, "scripts/matrix_to_h5.py")
    shell:
        """
        python3 {params.sct} --file {input} --output {output.h5} --header --dataset contigs --indexfile {output.idx}
        """

rule run_turbocor:
    input:
        os.path.join(STATS, "sample_coverage.h5")
    output:
        pearson = os.path.join(STATS, "sample_coverage.pearson"),
        tmp = temporary(os.path.join(STATS, "sample_coverage.corr"))
    conda:
        "../envs/turbocor.yaml"
    shell:
        """
        turbocor compute {input} {output.tmp} --dataset contigs;
        turbocor topk 1000000 {output.tmp} > {output.pearson}
        """

rule generate_clusters:
    input:
        os.path.join(STATS, "sample_coverage.pearson")
    output:
        os.path.join(STATS, "atavide_clusters.json")
    params:
        sct = os.path.join(ATAVIDE_DIR, "scripts/correlation_clustering.py"),
        pearson_threshold = 0.99
    shell:
        """
        python3 {params.sct} --file {input} --output {output} \
            --separator , --threshold {params.pearson_threshold}
        """




