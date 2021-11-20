"""
In the rules below we use the m8 files from superfocus
to crate taxonomy tables.
"""

rule superfocus_taxonomy:
    input:
        m8 = os.path.join(RBADIR, "{sample}", "superfocus", "{sample}.good_out_R1.fastq_alignments.m8"),
    output:
        os.path.join(RBADIR, "{sample}", "superfocus", "{sample}_good_out.taxonomy")
    params:
        t = TAXON
    threads: 4
    resources:
        mem_mb=16000
    conda:
        "../envs/taxonkit.yaml"
    shell:
        """
        perl -F"\\t" -lane 'if ($F[1] =~ /fig\|(\d+)\.\d+/ && $1 != 6666666) \
           {{print "$F[0]\\t$1"}}' {input.m8} | \
        taxonkit lineage -j {threads} -i 2 --data-dir {params.t} | \
        taxonkit reformat -j {threads} --data-dir {params.t} -i 3 \
          -f "Root; d__{{k}}; p__{{p}}; c__{{c}}; o__{{o}}; f__{{f}}; g__{{g}}" \
          --fill-miss-rank | \
        cut -f 1,2,4 > {output}
        """


