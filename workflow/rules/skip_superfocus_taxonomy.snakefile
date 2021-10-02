# this snakefile makes an nearly all taxonomy file that 
# satisfies snakemake but not the user!

temptxt = """
Sorry.

We could not make the superfocus taxonomy. 

You will need to download an sqlite taxonomy database from Rob. 
Check his github at
https://github.com/linsalrob/EdwardsLab/tree/master/taxon
or ask him for a new version!

Good luck!
"""

rule fake_superfocus_taxonomy:
    input:
        m8 = os.path.join(RBADIR, "{sample}", "superfocus", "{sample}.good_out_R1.fastq_alignments.m8"),
    output:
        os.path.join(RBADIR, "superfocus_taxonomy.tsv.zip")
    params:
        text = temptext
    shell:
        """
        echo {params.text} > {output}
        """

