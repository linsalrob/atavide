"""
The rules that generate statistics output
"""


rule final_assembly_stats:
    input:
        os.path.join(CCMO, "assembly.fasta")
    output:
        os.path.join(STATS, "final_assembly.txt")
    params:
        sct = os.path.join(PMSDIR, "scripts/countfasta.py")
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
        sct = os.path.join(PMSDIR, "scripts/joinlists.pl")
    shell:
        """
        perl {params.sct} -t {input} > {output}
        """

rule make_h5_table:
    input:
        expand(os.path.join(RMRD, "{smpl}_contig_hits.tsv"), smpl=SAMPLES)
    output:
        os.path.join(STATS, "sample_coverage.h5")
    resources:
        mem_mb=64000
    conda:
        "../envs/h5py.yaml"
    shell:
        """
         python3 ~/GitHubs/EdwardsLab/h5py/files_to_h5.py -f {input} -o {output} -s
         """

rule count_contig_lengths:
    input:
        os.path.join(CCMO, "assembly.fasta")
    output:
        os.path.join(STATS, "sequence_lengths.tsv")
    params:
        sct = os.path.join(PMSDIR, "scripts/countfasta.py")
    shell:
        """
        python3 {params.sct} -f {input} -ln > {output}
        """


