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
        db = "/home/edwa0468/ncbi/taxonomy/taxonomy.sqlite3"
    resources:
        mem_mb=16000
    shell:
        """
        python3 /home/edwa0468/GitHubs/atavide/workflow/scripts/superfocus_to_taxonomy.py -f {input} --tophit -d {params.db} > {output}
        """

rule count_sf_taxonomy:
    input:
        os.path.join(RBADIR, "{sample}", "superfocus", "{sample}_good_out.taxonomy")
    output:
        os.path.join(RBADIR, "{sample}", "superfocus", "{sample}_tax_counts.tsv")
    params:
        s = "{sample}"
    resources:
        cpus=4,
        mem_mb=16000
    shell:
        """
        cut -d$'\t' -f 2- {input} | sed -e 's/\t/|/g' | \
        awk -v F={params.s} 'BEGIN {{print "Superkingdom|Phylum|Class|Order|Family|Genus|Species|Taxid\t"F}} s[$0]++ {{}} END {{ for (i in s) print i"\t"s[i] }}' \
        > {output}
        """

rule join_superfocus_taxonomy:
    input:
        expand(os.path.join(RBADIR, "{smps}", "superfocus", "{smps}_tax_counts.tsv"), smps=SAMPLES)
    output:
        os.path.join(RBADIR, "superfocus_taxonomy.tsv")
    shell:
        """
        perl /home/edwa0468/GitHubs/atavide/workflow/scripts/joinlists.pl -h -z {input} > {output}
        """


rule compress_superfocus_taxonomy:
    input:
        os.path.join(RBADIR, "{sample}", "superfocus", "{sample}_good_out.taxonomy"),
        os.path.join(RBADIR, "{sample}", "superfocus", "{sample}_tax_counts.tsv")
    output:
        os.path.join(RBADIR, "{sample}", "superfocus", "{sample}_good_out.taxonomy.zip"),
        os.path.join(RBADIR, "{sample}", "superfocus", "{sample}_tax_counts.tsv.zip")
    shell:
        """
        for F in {input}; do zip $F.zip $F; done
        """

