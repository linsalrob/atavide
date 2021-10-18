"""
This set of rules handles the binning of reads using concoct and metabat

"""

import os
import sys
import socket



rule metabat_contigs_depth:
    input:
        bam = expand(os.path.join(RMRD, "{sample}.final_contigs.bam"), sample=SAMPLES),
        contigs = os.path.join(CCMO, "assembly.fasta")
    output:
        depth=os.path.join(METABAT, "metabat_depth")
    params:
        basename=os.path.join(METABAT, "metabat_bins")
    conda:
        "envs/metabat2.yaml"
    shell:
        """
        mkdir --parents {params.basename};
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bam}
        """

rule metabat_bins:
    input:
        contigs=os.path.join(CCMO, "assembly.fasta"),
        depth= os.path.join(METABAT, "metabat_depth")
    output:
        o=os.path.join(METABAT, "metabat_bins/metabat_bins.1.fa")
    params:
        base=os.path.join(METABAT, "metabat_bins/metabat_bins")
    conda:
        "envs/metabat2.yaml"
    resources:
        mem_mb=20000,
        cpus=16
    shell:
        """
        touch {output}
        metabat2 -i {input.contigs} -a {input.depth} -m 1500 -o {params.base} -t {resources.cpus}
        """

rule concoct:
    input:
       contigs= os.path.join(CCMO, "assembly.fasta"),
       tsv=os.path.join(STATS, "sample_coverage.tsv")
    output:
       odir=os.path.join(CONCOCT, "concoct_output/clustering_gt1000.csv")
    params:
       od=os.path.join(CONCOCT, "concoct_output")
    conda:
        "envs/concoct.yaml"
    resources:
        mem_mb=20000,
        cpus=16
    shell:
       """
       concoct --composition_file {input.contigs} --coverage_file {input.tsv} -b {params.od} -t {resources.cpus}
       """

rule extract_concoct_bins:
    input:
       contigs= os.path.join(CCMO, "assembly.fasta"),
       odir=os.path.join(CONCOCT, "concoct_output/clustering_gt1000.csv")
    output:
       os.path.join(CONCOCT, "concoct_bins/1.fa")
    params:
       base=os.path.join(CONCOCT, "concoct_bins")
    conda:
        "envs/concoct.yaml"
    shell:
       """
       touch {output}
       extract_fasta_bins.py {input.contigs} {input.odir} --output_path {params.base}
       """      
