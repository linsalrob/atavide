"""

Here we take all the unassembled reads and reassemble them again
"""


rule assemble_unassembled_gpu:
    """
    assemble the unassembled reads using the gpu.

    See the note in round 1 assemblies. Basically, they don't really
    support GPUs so we've switched back to CPU.
    """
    input:
        r1 = os.path.join(UNASSM, "R1.unassembled.fastq.gz"),
        r2 = os.path.join(UNASSM, "R2.unassembled.fastq.gz"),
        s0 = os.path.join(UNASSM, "single.unassembled.fastq.gz")
    output:
        os.path.join(REASSM, "gpu", "final.contigs.fa"),
        os.path.join(REASSM, "gpu", "log"),
        temporary(os.path.join(REASSM, "gpu", "checkpoints.txt")),
        temporary(os.path.join(REASSM, "gpu", "done")),
        temporary(directory(os.path.join(REASSM, "gpu", "intermediate_contigs"))),
        temporary(os.path.join(REASSM, "gpu", "options.json"))
    params:
        odir = os.path.join(REASSM, "gpu")
    conda:
        "../envs/megahit.yaml"
    threads: 32
    resources:
        mem_mb=128000,
        time=7200,
        gpu=1,
        partition="gpu"
    shell:
        """
        rmdir {params.odir};
        megahit -1 {input.r1} -2 {input.r2} -r {input.s0} -o {params.odir} -t {threads}  --use-gpu --mem-flag 2
        """

rule assemble_unassembled:
    """
    assemble the unassembled reads using CPU (but more threads)
    """
    input:
        r1 = os.path.join(UNASSM, "R1.unassembled.fastq.gz"),
        r2 = os.path.join(UNASSM, "R2.unassembled.fastq.gz"),
        s0 = os.path.join(UNASSM, "single.unassembled.fastq.gz")
    output:
        os.path.join(REASSM, "final.contigs.fa"),
        os.path.join(REASSM, "log"),
        temporary(os.path.join(REASSM, "checkpoints.txt")),
        temporary(os.path.join(REASSM, "done")),
        temporary(directory(os.path.join(REASSM, "intermediate_contigs"))),
        temporary(os.path.join(REASSM, "options.json"))
    params:
        odir = os.path.join(REASSM)
    conda:
        "../envs/megahit.yaml"
    threads: 32
    resources:
        mem_mb=LARGE_MEM,
        time=7200,
    shell:
        """
        rmdir {params.odir};
        megahit -1 {input.r1} -2 {input.r2} -r {input.s0} -o {params.odir} -t {threads} --mem-flag 2
        """

"""
Combine all the contigs and use metaflye to merge the subassmblies
"""


rule concatenate_all_assemblies:
    """
    Again we take all contigs produced by megahit and combine them into a single (redundant)
    fasta file. This also creates a file called merged_contigs.ids that has the 
    contig ids. 
    """

    input:
        expand(os.path.join(ASSDIR, "{sample}/final.contigs.fa"), sample=SAMPLES),
        os.path.join(REASSM, "final.contigs.fa")
    output:
        contigs = os.path.join(REASSM, "merged_contigs.fa"),
        ids = os.path.join(REASSM, "merged_contigs.ids") 
    resources:
        mem_mb=128000,
    script:
        "../scripts/renumber_merge_fasta_smk.py"

rule merge_assemblies_with_flye:
    """
    Run flye on all the merged contigs fromt the last round to merge everything one more time
    """
    input:
        contigs = os.path.join(REASSM, "merged_contigs.fa")
    output:
        os.path.join(CCMO, "assembly.fasta"),
        os.path.join(CCMO, "assembly_graph.gfa"),
        os.path.join(CCMO, "assembly_graph.gv"),
        os.path.join(CCMO, "assembly_info.txt"),
        os.path.join(CCMO, "flye.log"),
    conda:
        "../envs/flye.yaml"
    threads: 32
    resources:
        mem_mb=128000,
        time=7200,
    shell:
        """
        flye --meta --subassemblies {input.contigs} -o {CCMO} --threads {threads}
        """

rule renumber_rename_contigs:
    """
    make a single assembly file in the main output directory with the contigs
    labeled with the sample_id
    """
    input:
        os.path.join(CCMO, "assembly.fasta")
    output:
        contigs = os.path.join(config['directories']['assemblies'], f"{SAMPLE_ID}_assembly.fasta"),
        ids = os.path.join(config['directories']['assemblies'], f"{SAMPLE_ID}_assembly.original_ids")
    params:
        sample_id = SAMPLE_ID
    resources:
        mem_mb=128000,
    script:
        "../scripts/renumber_merge_fasta_smk.py"

rule index_final_contigs:
    """
    build the bowtie2 index of the final assembly
    """
    input:
        os.path.join(config['directories']['assemblies'], f"{SAMPLE_ID}_assembly.fasta")
    params:
        baseoutput = os.path.join(config['directories']['assemblies'], "bowtie2_index", f"{SAMPLE_ID}_assembly")
    output:
        idx1 = os.path.join(config['directories']['assemblies'], "bowtie2_index", f"{SAMPLE_ID}_assembly.1.bt2l"),
        idx2 = os.path.join(config['directories']['assemblies'], "bowtie2_index", f"{SAMPLE_ID}_assembly.2.bt2l"),
        idx3 = os.path.join(config['directories']['assemblies'], "bowtie2_index", f"{SAMPLE_ID}_assembly.3.bt2l"),
        idx4 = os.path.join(config['directories']['assemblies'], "bowtie2_index", f"{SAMPLE_ID}_assembly.4.bt2l"),
        rid1 = os.path.join(config['directories']['assemblies'], "bowtie2_index", f"{SAMPLE_ID}_assembly.rev.1.bt2l"),
        rid2 = os.path.join(config['directories']['assemblies'], "bowtie2_index", f"{SAMPLE_ID}_assembly.rev.2.bt2l")
    threads: 16
    resources:
        mem_mb=64000,
    conda:
        "../envs/bowtie.yaml"
    shell:
        # note that we build a large index by default, because
        # at some point we will end up doing that, and so this
        # always makes a large index
        """
        mkdir -p {RMRD} && \
        bowtie2-build --threads {threads} --large-index {input} {params.baseoutput}
        """
        

rule map_reads_to_final:
    """
    Here we map the original reads back to the contigs
    """
    input:
        idx1 = os.path.join(config['directories']['assemblies'], "bowtie2_index", f"{SAMPLE_ID}_assembly.1.bt2l"),
        idx2 = os.path.join(config['directories']['assemblies'], "bowtie2_index", f"{SAMPLE_ID}_assembly.2.bt2l"),
        idx3 = os.path.join(config['directories']['assemblies'], "bowtie2_index", f"{SAMPLE_ID}_assembly.3.bt2l"),
        idx4 = os.path.join(config['directories']['assemblies'], "bowtie2_index", f"{SAMPLE_ID}_assembly.4.bt2l"),
        rid1 = os.path.join(config['directories']['assemblies'], "bowtie2_index", f"{SAMPLE_ID}_assembly.rev.1.bt2l"),
        rid2 = os.path.join(config['directories']['assemblies'], "bowtie2_index", f"{SAMPLE_ID}_assembly.rev.2.bt2l"),
        r1 = os.path.join(PSEQDIR_TWO, "{sample}_good_out_R1.fastq"),
        r2 = os.path.join(PSEQDIR_TWO, "{sample}_good_out_R2.fastq"),
        s1 = os.path.join(PSEQDIR_TWO, "{sample}_single_out_R1.fastq"),
        s2 = os.path.join(PSEQDIR_TWO, "{sample}_single_out_R2.fastq")
    params:
        contigs = os.path.join(config['directories']['assemblies'], "bowtie2_index", f"{SAMPLE_ID}_assembly")
    output:
        os.path.join(RMRD, "{sample}." + SAMPLE_ID + ".assembly.bam")
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

rule bai_final_bams:
    input:
        os.path.join(RMRD, "{sample}." + SAMPLE_ID + ".assembly.bam")
    output:
        os.path.join(RMRD, "{sample}." + SAMPLE_ID + ".assembly.bam.bai")
    conda:
        "../envs/bowtie.yaml"
    shell:
        "samtools index {input}"

rule count_mapped_reads:
    input:
        bam = os.path.join(RMRD, "{sample}." + SAMPLE_ID + ".assembly.bam"),
        bai = os.path.join(RMRD, "{sample}." + SAMPLE_ID + ".assembly.bam.bai")
    output:
        os.path.join(RMRD, "{sample}_contig_hits.tsv")
    threads: 8
    resources:
        mem_mb=8000,
    conda:
        "../envs/bowtie.yaml"
    shell:
        """
        samtools idxstats -@ {threads} {input.bam} | cut -f 1,3 > {output}
        """

rule count_rpkm:
    input:
        r1 = os.path.join(PSEQDIR_TWO, "{sample}_good_out_R1.fastq"),
        contig_len = os.path.join(STATS, "sequence_lengths.tsv"),
        hits = os.path.join(RMRD, "{sample}_contig_hits.tsv")
    output:
        rpkm = os.path.join(RMRD, "{sample}_rpkm.tsv")
    resources:
        mem_mb=8000,
    script:
        "../scripts/rpkm.py"



