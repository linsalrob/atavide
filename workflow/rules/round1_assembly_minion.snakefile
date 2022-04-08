"""

This round does a megahit assembly, and then maps the reads back to all
the contigs. We then find unassembled reads and concatenate them
"""




rule megahit_assemble_gpu:
    """
    Note this is how to do the assembly using GPUs, but
    according to megahit they don't really support this
    so I bypass it and use CPUs
    """
    input:
        r1 = os.path.join(PSEQDIR_TWO, FASTQ),
    output:
        os.path.join(ASSDIR, "{sample}_gpu/final.contigs.fa"),
        os.path.join(ASSDIR, "{sample}_gpu/log"),
        temporary(os.path.join(ASSDIR, "{sample}_gpu/checkpoints.txt")),
        temporary(os.path.join(ASSDIR, "{sample}_gpu/done")),
        temporary(directory(os.path.join(ASSDIR, "{sample}_gpu/intermediate_contigs"))),
        temporary(os.path.join(ASSDIR, "{sample}_gpu/options.json"))
    params:
        odir = directory(os.path.join(ASSDIR, '{sample}_gpu'))
    threads: 16
    resources:
        mem_mb=64000,
        time=7200,
        gpu=1,
        partition="gpu"
    conda:
        "../envs/megahit.yaml"
    shell:
        """
        rmdir {params.odir} ; 
        megahit -r {input.r1} -o {params.odir} -t {threads} --use-gpu --mem-flag 2
        """

rule megahit_assemble:
    """
    This version uses CPUs to do the assembly
    """
    input:
        r1 = os.path.join(PSEQDIR_TWO, FASTQ),
    output:
        os.path.join(ASSDIR, "{sample}/final.contigs.fa"),
        os.path.join(ASSDIR, "{sample}/log"),
        temporary(os.path.join(ASSDIR, "{sample}/checkpoints.txt")),
        temporary(os.path.join(ASSDIR, "{sample}/done")),
        temporary(directory(os.path.join(ASSDIR, "{sample}/intermediate_contigs"))),
        temporary(os.path.join(ASSDIR, "{sample}/options.json"))
    params:
        odir = directory(os.path.join(ASSDIR, '{sample}'))
    threads: 32
    resources:
        mem_mb=64000,
        time=7200,
    conda:
        "../envs/megahit.yaml"
    shell:
        """
        rmdir {params.odir} ; 
        megahit -r {input.r1} -o {params.odir} -t {threads} --mem-flag 2
        """

rule combine_contigs:
    """
    Here we take all contigs produced by megahit and combine them into a single (redundant)
    fasta file. This also creates a file called round1_contigs.ids that has the original
    contig ids
    """

    input:
        expand(os.path.join(ASSDIR, "{sample}/final.contigs.fa"), sample=SAMPLES)
    output:
        contigs = os.path.join(ASSDIR, "round1_contigs.fa"),
        ids = os.path.join(ASSDIR, "round1_contigs.ids")
    resources:
        mem_mb=128000,
    script:
        "../scripts/renumber_merge_fasta_smk.py"

rule index_contigs:
    """
    build the bowtie2 index of the contigs
    """
    input:
        os.path.join(ASSDIR, "round1_contigs.fa")
    params:
        baseoutput = os.path.join(CRMDIR, "round1_contigs")
    output:
        idx1 = temporary(os.path.join(CRMDIR, "round1_contigs.1.bt2l")),
        idx2 = temporary(os.path.join(CRMDIR, "round1_contigs.2.bt2l")),
        idx3 = temporary(os.path.join(CRMDIR, "round1_contigs.3.bt2l")),
        idx4 = temporary(os.path.join(CRMDIR, "round1_contigs.4.bt2l")),
        ridx1 = temporary(os.path.join(CRMDIR, "round1_contigs.rev.1.bt2l")),
        ridx2 = temporary(os.path.join(CRMDIR, "round1_contigs.rev.2.bt2l"))
    threads: 8
    resources:
        mem_mb=64000,
    conda:
        "../envs/bowtie.yaml"
    shell:
        # note that we build a large index by default, because
        # at some point we will end up doing that, and so this
        # always makes a large index
        """
        mkdir -p {CRMDIR} && \
        bowtie2-build --threads {threads} --large-index {input} {params.baseoutput}
        """
        

rule map_reads:
    """
    Here we map the original reads back to the contigs
    """
    input:
        idx1 = os.path.join(CRMDIR, "round1_contigs.1.bt2l"),
        idx2 = os.path.join(CRMDIR, "round1_contigs.2.bt2l"),
        idx3 = os.path.join(CRMDIR, "round1_contigs.3.bt2l"),
        idx4 = os.path.join(CRMDIR, "round1_contigs.4.bt2l"),
        ridx1 = os.path.join(CRMDIR, "round1_contigs.rev.1.bt2l"),
        ridx2 = os.path.join(CRMDIR, "round1_contigs.rev.2.bt2l"),
        r1 = os.path.join(PSEQDIR_TWO, FASTQ)
    params:
        contigs = os.path.join(CRMDIR, "round1_contigs")
    output:
        os.path.join(CRMDIR, "{sample}.contigs.bam")
    threads: 8
    resources:
        mem_mb=20000,
    conda:
        "../envs/bowtie.yaml"
    shell:
        """
        bowtie2 --mm -x {params.contigs} -U {input.r1} --threads {threads} | \
        samtools view -bh | samtools sort -o {output} -
        """

"""
Using samtools to extract unmapped reads from the bam files

For the minion reads, we just want the unmapped reads (-f 4)

See: https://edwards.sdsu.edu/research/command-line-deconseq/
"""

rule umapped_reads:
    """
    get left reads
    """
    input:
        os.path.join(CRMDIR, "{sample}.contigs.bam")
    output:
        os.path.join(UNASSM, "{sample}.unassembled.fastq")
    conda:
        "../envs/bowtie.yaml"
    threads: 8
    resources:
        mem_mb=32000,
    shell:
        """
        samtools view -@ {threads} -h {input} | 
                awk 'BEGIN {{FS="\t"; OFS="\t"}} 
                {{if (/^@/ && substr($2, 3, 1)==":") {{print}} 
                else if (and($2, 0x40)) {{print}}}}' \
                | samtools bam2fq -@ {threads} > {output}
        """


"""
We concatanate the unassembled reads so we can assmble them all together.

"""


rule compress_unassembled_reads:
    input:
        os.path.join(UNASSM, "{sample}.unassembled.fastq"),
    output:
        temporary(os.path.join(UNASSM, "{sample}.unassembled.fastq.gz")),
    shell:
        """
        for F in {input}; do gzip $F; done
        """



rule concatenate_unassembled:
    """
    Concat reads
    """
    input:
        expand(os.path.join(UNASSM, "{sample}.unassembled.fastq.gz"), sample=SAMPLES)
    output:
        os.path.join(UNASSM, "unassembled.fastq.gz")
    shell:
        "cat {input} > {output}"

