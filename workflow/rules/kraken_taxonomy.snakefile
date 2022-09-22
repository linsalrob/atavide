



rule run_kraken:
    input:
        r1 = os.path.join(PSEQDIR_TWO, "{sample}_good_out_R1.fastq"),
        r2 = os.path.join(PSEQDIR_TWO, "{sample}_good_out_R2.fastq")
    output:
        rt = os.path.join(RBADIR, "{sample}", "kraken", "{sample}.report.tsv"),
        ot = os.path.join(RBADIR, "{sample}", "kraken", "{sample}.output.tsv")
    threads: 8
    resources:
        mem_mb=250000,
        load_kraken=25
    conda:
        "../envs/kraken.yaml"
    shell:
        """
        kraken2 --report {output.rt} \
                --output {output.ot} \
                --threads {threads} \
                {input.r1}
        """



rule kraken_taxonomy:
    input:
        os.path.join(RBADIR, "{sample}", "kraken", "{sample}.output.tsv")
    output:
        os.path.join(RBADIR, "{sample}", "kraken", "{sample}.taxonomy.tsv")
    threads: 4
    resources:
        mem_mb=40000
    params:
        t = TAXON
    conda:
        "../envs/taxonkit.yaml"
    shell:
        """
        grep -v ^U {input} | \
        taxonkit lineage -j {threads} -i 3 --data-dir {params.t} | \
        taxonkit reformat -j {threads} --data-dir {params.t} -i 6 -f \
        "Root; d__{{k}}; p__{{p}}; c__{{c}}; o__{{o}}; f__{{f}}; g__{{g}}" \
        --fill-miss-rank | cut -f 2,3,7 > {output}
        """

rule kraken_phyla:
    # extract the kraken families into a two column table
    input:
        os.path.join(RBADIR, "{sample}", "kraken", "{sample}.report.tsv")
    output:
        os.path.join(RBADIR, "{sample}", "kraken", "{sample}.phyla_fraction.tsv")
    shell:
        """
        echo "Family\t{wildcards.sample}" > {output};
        cat {input} | sed -e 's/Candidatus//' | \
        awk '{{if ($4 == "P") {{printf "%s\\t%f\\n", $6, $1}}; }}' >> {output}
        """


rule kraken_families:
    # extract the kraken families into a two column table
    input:
        os.path.join(RBADIR, "{sample}", "kraken", "{sample}.report.tsv")
    output:
        os.path.join(RBADIR, "{sample}", "kraken", "{sample}.family_fraction.tsv")
    shell:
        """
        echo "Family\t{wildcards.sample}" > {output};
        cat {input} | sed -e 's/Candidatus//' | \
        awk '{{if ($4 == "F") {{printf "%s\\t%f\\n", $6, $1}}; }}' >> {output}
        """


rule kraken_genera:
    # extract the kraken families into a two column table
    input:
        os.path.join(RBADIR, "{sample}", "kraken", "{sample}.report.tsv")
    output:
        os.path.join(RBADIR, "{sample}", "kraken", "{sample}.genera_fraction.tsv")
    shell:
        """
        echo "Family\t{wildcards.sample}" > {output};
        cat {input} | sed -e 's/Candidatus//' | \
        awk '{{if ($4 == "G") {{printf "%s\\t%f\\n", $6, $1}}; }}' >> {output}
        """

rule kraken_species:
    # extract the kraken families into a two column table
    input:
        os.path.join(RBADIR, "{sample}", "kraken", "{sample}.report.tsv")
    output:
        os.path.join(RBADIR, "{sample}", "kraken", "{sample}.species_fraction.tsv")
    shell:
        """
        echo "Species\t{wildcards.sample}" > {output};
        cat {input} | sed -e 's/Candidatus//' | \
        awk '{{if ($4 == "S") {{for (i=6; i<=NF; ++i) {{printf "%s ", $i}} printf "\\t%f\\n", $1}}; }}' \
        >> {output}
        """


rule join_kraken_species:
    input:
        expand(os.path.join(RBADIR, "{sample}", "kraken", "{sample}.species_fraction.tsv"), sample=SAMPLES)
    output:
        os.path.join(STATS, "kraken_species.tsv")
    params:
        sct = os.path.join(ATAVIDE_DIR, "scripts/joinlists.pl")
    shell:
        """
        perl {params.sct} -z -h {input} > {output}
        """


rule join_kraken_phyla:
    input:
        expand(os.path.join(RBADIR, "{sample}", "kraken", "{sample}.phyla_fraction.tsv"), sample=SAMPLES)
    output:
        os.path.join(STATS, "kraken_phyla.tsv")
    params:
        sct = os.path.join(ATAVIDE_DIR, "scripts/joinlists.pl")
    shell:
        """
        perl {params.sct} -z -h {input} > {output}
        """

rule join_kraken_families:
    input:
        expand(os.path.join(RBADIR, "{sample}", "kraken", "{sample}.family_fraction.tsv"), sample=SAMPLES)
    output:
        os.path.join(STATS, "kraken_families.tsv")
    params:
        sct = os.path.join(ATAVIDE_DIR, "scripts/joinlists.pl")
    shell:
        """
        perl {params.sct} -z -h {input} > {output}
        """

rule join_kraken_genera:
    input:
        expand(os.path.join(RBADIR, "{sample}", "kraken", "{sample}.genera_fraction.tsv"), sample=SAMPLES)
    output:
        os.path.join(STATS, "kraken_genera.tsv")
    params:
        sct = os.path.join(ATAVIDE_DIR, "scripts/joinlists.pl")
    shell:
        """
        perl {params.sct} -z -h {input} > {output}
        """
