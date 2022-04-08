############################################################
#                                                          #
# A snakefile to remove host contamination.                #
#                                                          #
# We have several versions of contamination removal        #
# this uses bowtie2.                                       #
#                                                          #
# For more information see this site:                      #
# https://edwards.sdsu.edu/research/command-line-deconseq/ #
#                                                          #
############################################################

import os

# set this to the path to the indexed host genome
# that you want to filter out. This should be indexed
# with `bowtie2-build`, and should be just the base 
# filename of the .1.bt2, .2.bt2, indices or the bt2l for large


host_bt_index = os.path.join(config['directories']['host_dbpath'], config['options']['host_dbname'])

# set this to the location where you would like the results
# written. We will make the directory if it doesn't exist
INTERIMDIR = "host_mapped"


if not os.path.exists(host_bt_index + ".1.bt2") and not os.path.exists(host_bt_index + ".1.bt2l"):
    sys.stderr.write(f"ERROR: Could not find the bowtie2 indexes for {host_bt_index}\n")
    sys.stderr.write(f" - Did you build them with bowtie2-build?\n")
    sys.stderr.write(f" - You may need to edit `host_dbpath` and `host_dbname` in the snakefile config\n")
    sys.exit()


rule btmap:
    input:
        r1 = os.path.join(READDIR, FASTQ),
    output:
        os.path.join(INTERIMDIR, '{sample}.hostmapped.bam')
    params:
        idx = host_bt_index
    conda:
        "../envs/bowtie.yaml"
    threads: 16
    resources:
        mem_mb=64000,
    shell:
        """
		bowtie2 --mm -p {threads} -x {params.idx} -U {input.r1} \
         | samtools view -@ {threads} -bh | samtools sort -o {output} -
        """

rule reads_mapped:
    input:
        os.path.join(INTERIMDIR, '{sample}.hostmapped.bam')
    output:
        os.path.join(INTERIMDIR, "{sample}_R1_mapped_to_host.fastq")
    threads: 8
    resources:
        mem_mb=16000,
    conda:
        "../envs/bowtie.yaml"
    shell:
        "samtools fastq -@ {threads} -G 4 {input} > {output}"

        
rule reads_unmapped:
    input:
        os.path.join(INTERIMDIR, '{sample}.hostmapped.bam')
    output:
        r1 = os.path.join(PSEQDIR_TWO, FASTQ),
    threads: 8
    resources:
        mem_mb=16000,
    conda:
        "../envs/bowtie.yaml"
    shell:
        "samtools fastq -@ {threads} -f 4  {input} > {output}"


rule compress_host_mapped_sequences:
    input:
        r1 = os.path.join(PSEQDIR_TWO, FASTQ),
    output:
        r1 = os.path.join(PSEQDIR_TWO, f"{FASTQ}.gz"),
    shell:
        """
        for F in {input}; do gzip $F; done
        """


