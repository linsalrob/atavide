#####################################################################
#                                                                   #
# Run rarefactions on the kraken2 species predictions.              #
#                                                                   #
# We subsample the fastq with replacement, and use those samples    #
# to run Kraken2. Then we identify the number of species            #
# at each fraction, and put that in a single table                  #
#                                                                   #
#                                                                   #
#                                                                   #
#####################################################################



# how many rarefactions are we going to run?
FRACTIONS = [i/10 for i in range(1, 10)]




rule subsample_fastq:
    input:
        r1 = os.path.join(PSEQDIR_TWO, "{sample}_good_out_R1.fastq")
    output:
        temporary(os.path.join(RBADIR, "{sample}", "kraken", "{sample}_good_out_R1.{frac}.fastq"))
    threads: 8
    resources:
        mem_mb=25000
    conda:
        "../envs/seqtk.yaml"
    params:
        f = "{frac}"
    shell:
        """
        seqtk sample {input.r1} {params.f} > {output}
        """

rule kraken_subsample:
    input:
        r1 = os.path.join(RBADIR, "{sample}", "kraken", "{sample}_good_out_R1.{frac}.fastq")
    output:
        rt = os.path.join(RBADIR, "{sample}", "kraken", "{sample}.report.{frac}.tsv"),
    threads: 8
    resources:
        mem_mb=256000,
        load_kraken=50
    params:
        directory = os.path.join(RBADIR, "{sample}", "kraken")
    conda:
        "../envs/kraken.yaml"
    shell:
        """
        kraken2 --report {output.rt} \
                --threads {threads} \
                {input.r1}
        """

"""

I think the next three rules will make this work, but need to double check
the rule. Essentially we just mix some bash in here too!

cd RBADIR
WD=$PWD; for SEQ in FAME0000*; do echo $SEQ; cd $SEQ/kraken; \
        echo -e "Fraction\t$SEQ" > $SEQ.kraken_rarefaction.tsv; \
        for FRX in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9; \
        do awk -v F=$FRX '$4 == "S" {c++} END {print F,"\t",c}' $SEQ.report.$FRX.tsv; \
        done >> $SEQ.kraken_rarefaction.tsv; \
        awk '$4 == "S" {c++} END {print "1.0\t",c}' $SEQ.report.tsv >> $SEQ.kraken_rarefaction.tsv; \
        cd $WD; done

joinlists.pl -h */kraken/*kraken_rarefaction.tsv > ../statistics/kraken_rarefaction.tsv
"""


rule kraken_summarize_species_ready:
    # this rule just forces all the outputs to complete before we summarize them!
    input:
        expand(os.path.join(RBADIR, "{sample}", "kraken", "{sample}.report.{frac}.tsv"), sample=SAMPLES, frac=FRACTIONS)
    output:
        temporary(os.path.join(STATS, "ready_to_summarize"))
    shell:
        """
        touch {output}
        """


rule kraken_summarize_species:
    input:
        os.path.join(STATS, "ready_to_summarize")
    output:
        os.path.join(RBADIR, "{sample}", "kraken", "{sample}.kraken_species_rarefaction.tsv")
    params:
        frx = FRACTIONS,
        smp = os.path.join(RBADIR, "{sample}", "kraken", "{sample}"),
        sample = "{sample}"
    shell:
        """
        echo -e "Fraction\t{params.sample}" > {output};
        for FRX in {params.frx}; do awk -v F=$FRX '$4 == "S" {{c++}} END {{print F,"\t",c}}' {params.smp}.report.$FRX.tsv; done >> {output};
        awk '$4 == "S" {{c++}} END {{print "1.0\t",c}}' {params.smp}.report.tsv >> {output}
        """



rule join_kraken_subsamples:
    input:
        expand(os.path.join(RBADIR, "{sample}", "kraken", "{sample}.kraken_species_rarefaction.tsv"), sample=SAMPLES)
    output:
        os.path.join(STATS, "kraken_species_rarefaction.tsv")
    params:
        sct = os.path.join(ATAVIDE_DIR, "scripts/joinlists.pl")
    shell:
        """
        perl {params.sct} -h {input} > {output}
        """


