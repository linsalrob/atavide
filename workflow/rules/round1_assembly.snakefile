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
        r1 = os.path.join(PSEQDIR_TWO, "{sample}_good_out_R1.fastq"),
        r2 = os.path.join(PSEQDIR_TWO, "{sample}_good_out_R2.fastq"),
        s1 = os.path.join(PSEQDIR_TWO, "{sample}_single_out_R1.fastq"),
        s2 = os.path.join(PSEQDIR_TWO, "{sample}_single_out_R2.fastq")
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
        megahit -1 {input.r1} -2 {input.r2} -r {input.s1} -r {input.s2} \
                -o {params.odir} -t {threads} --use-gpu --mem-flag 2
        """

rule megahit_assemble:
    """
    This version uses CPUs to do the assembly
    """
    input:
        r1 = os.path.join(PSEQDIR_TWO, "{sample}_good_out_R1.fastq"),
        r2 = os.path.join(PSEQDIR_TWO, "{sample}_good_out_R2.fastq"),
        s1 = os.path.join(PSEQDIR_TWO, "{sample}_single_out_R1.fastq"),
        s2 = os.path.join(PSEQDIR_TWO, "{sample}_single_out_R2.fastq")
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
        megahit -1 {input.r1} -2 {input.r2} -r {input.s1} -r {input.s2} \
                -o {params.odir} -t {threads} --mem-flag 2
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
        r1 = os.path.join(PSEQDIR_TWO, "{sample}_good_out_R1.fastq"),
        r2 = os.path.join(PSEQDIR_TWO, "{sample}_good_out_R2.fastq"),
        s1 = os.path.join(PSEQDIR_TWO, "{sample}_single_out_R1.fastq"),
        s2 = os.path.join(PSEQDIR_TWO, "{sample}_single_out_R2.fastq")
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
        bowtie2 --mm -x {params.contigs} -1 {input.r1} -2 {input.r2} \
        -U {input.s1} -U {input.s2} --threads {threads} | \
        samtools view -bh | samtools sort -o {output} -
        """

"""
Using samtools to extract unmapped reads from the bam files
A samtools flag of -f 77 means the read is paird, neither the read nor mate is mapped, and the it is the first read in the pair, while a flag of -f 141 means the same except it is the second mate in the pair.
Then we use two flags: -f 4 (the read is unmapped) and -F 1 (the read is not paired) to find the single reads that are not mapped

See: https://edwards.sdsu.edu/research/command-line-deconseq/
"""

rule umapped_left_reads:
    """
    get left reads
    """
    input:
        os.path.join(CRMDIR, "{sample}.contigs.bam")
    output:
        os.path.join(UNASSM, "{sample}.unassembled.R1.fastq")
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
                else if (and($2, 0x1) && and($2, 0x40) && 
                (and($2, 0x4) || and($2, 0x8))) {{print}}}}' \
                | samtools bam2fq -@ {threads} > {output}
        """

rule umapped_right_reads:
    """
    get right reads
    """
    input:
        os.path.join(CRMDIR, "{sample}.contigs.bam")
    output:
        os.path.join(UNASSM, "{sample}.unassembled.R2.fastq")
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
                else if (and($2, 0x1) && and($2, 0x80) && 
                (and($2, 0x4) || and($2, 0x8))) {{print}}}}' \
                | samtools bam2fq -@ {threads} > {output}
        """

rule umapped_single_reads:
    """
    get singletons
    """
    input:
        os.path.join(CRMDIR, "{sample}.contigs.bam")
    output:
        os.path.join(UNASSM, "{sample}.unassembled.singles.fastq")
    conda:
        "../envs/bowtie.yaml"
    threads: 8
    resources:
        mem_mb=32000,
    shell:
        "samtools fastq -@ {threads} -f 4 -F 1 {input} > {output}"

"""
We concatanate the unassembled reads into separate R1/R2/s files so
we can assmble them all together.

We also do this in 3 separate threads to take advantage of parallelization
and to make the command easier
"""


rule compress_unassembled_reads:
    input:
        os.path.join(UNASSM, "{sample}.unassembled.R1.fastq"),
        os.path.join(UNASSM, "{sample}.unassembled.R2.fastq"),
        os.path.join(UNASSM, "{sample}.unassembled.singles.fastq")
    output:
        temporary(os.path.join(UNASSM, "{sample}.unassembled.R1.fastq.gz")),
        temporary(os.path.join(UNASSM, "{sample}.unassembled.R2.fastq.gz")),
        temporary(os.path.join(UNASSM, "{sample}.unassembled.singles.fastq.gz"))
    shell:
        """
        for F in {input}; do gzip $F; done
        """



rule concatenate_R1_unassembled:
    """
    Concat R1 reads
    """
    input:
        expand(os.path.join(UNASSM, "{sample}.unassembled.R1.fastq.gz"), sample=SAMPLES)
    output:
        os.path.join(UNASSM, "R1.unassembled.fastq.gz")
    shell:
        "cat {input} > {output}"

rule concatenate_R2_unassembled:
    """
    Concat R2 reads
    """
    input:
        expand(os.path.join(UNASSM, "{sample}.unassembled.R2.fastq.gz"), sample=SAMPLES)
    output:
        os.path.join(UNASSM, "R2.unassembled.fastq.gz")
    shell:
        "cat {input} > {output}"

rule concatenate_single_unassembled:
    """
    Concate singletons
    """
    input:
        expand(os.path.join(UNASSM, "{sample}.unassembled.singles.fastq.gz"),sample=SAMPLES)
    output:
        os.path.join(UNASSM, "single.unassembled.fastq.gz")
    shell:
        "cat {input} > {output}"


