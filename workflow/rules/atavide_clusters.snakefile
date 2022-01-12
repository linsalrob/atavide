
rule make_h5_table:
    input:
        os.path.join(STATS, "sample_coverage.tsv")
    output:
        h5 = os.path.join(ATAVIDE_BINNING, "stats", "sample_coverage.h5"),
        idx = os.path.join(ATAVIDE_BINNING, "stats", "sample_coverage.contig_index.tsv")
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
        os.path.join(ATAVIDE_BINNING, "stats", "sample_coverage.h5")
    output:
        pearson = os.path.join(ATAVIDE_BINNING, "stats", "sample_coverage.pearson"),
        tmp = os.path.join(ATAVIDE_BINNING, "stats", "sample_coverage.corr")
    conda:
        "../envs/turbocor.yaml"
    shell:
        """
        turbocor compute {input} {output.tmp} --dataset contigs;
        turbocor topk 1000000 {output.tmp} > {output.pearson}
        """

rule generate_clusters:
    input:
        os.path.join(ATAVIDE_BINNING, "stats", "sample_coverage.pearson")
    output:
        os.path.join(ATAVIDE_BINNING, "stats", "atavide_clusters.json")
    resources:
        mem_mb=64000
    params:
        sct = os.path.join(ATAVIDE_DIR, "scripts/correlation_clustering.py"),
        pearson_threshold = 0.99
    shell:
        """
        python3 {params.sct} --file {input} --output {output} \
            --separator , --threshold {params.pearson_threshold}
        """


rule extract_fasta_sequences:
    input:
        fa = os.path.join(config['directories']['assemblies'], f"{SAMPLE_ID}_assembly.fasta"),
        cl = os.path.join(ATAVIDE_BINNING, "stats", "atavide_clusters.json"),
        idx = os.path.join(ATAVIDE_BINNING, "stats", "sample_coverage.contig_index.tsv")
    output:
        directory = directory(os.path.join(ATAVIDE_BINNING, "bins"))
    params:
        sct = os.path.join(ATAVIDE_DIR, "scripts/correlation_clusters_to_fasta.py"),
        pearson_threshold = 0.97
    shell:
        """
        python3 {params.sct} --fasta {input.fa} --clusters {input.cl} \
        --directory {output.directory} --idx {input.idx}
        """


