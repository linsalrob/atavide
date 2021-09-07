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
# filename of the .1.bt2, .2.bt2, .3.bt2 indexes etc.
host_bt_index = "/path/to/human"   

# set this to the host name that we will use in the output names
hostname = "human"


# set this to the location where you would like the results
# written. We will make the directory if it doesn't exist
outdir = "host_mapped"


if not os.path.exists(host_bt_index + ".1.bt2") and not os.path.exists(host_bt_index + ".1.bt2l"):
    sys.stderr.write(f"ERROR: Could not find the bowtie2 indexes for {host_bt_index}\n")
    sys.stderr.write(f" - Did you build them with bowtie2-build?\n")
    sys.stderr.write(f" - You may need to edit `host_bt_index` in the snakefile\n")
    sys.exit()


rule btmap:
    input:
        r1 = os.path.join(PSEQDIR, "{sample}_good_out_R1.fastq"),
        r2 = os.path.join(PSEQDIR, "{sample}_good_out_R2.fastq"),
        s1 = os.path.join(PSEQDIR, "{sample}_single_out_R1.fastq"),
        s2 = os.path.join(PSEQDIR, "{sample}_single_out_R2.fastq"),
    output:
        os.path.join(outdir, '{sample}.' + hostname + '.bam')
    params:
        idx = host_bt_index
    conda:
        "envs/bowtie.yaml"
    shell:
        """
		bowtie2 -p {threads} -x {params.idx} -1 {input.r1} -2 {input.r2} | samtools view -bh | samtools sort -o {output} -
        """

rule R1_reads_map_to_ref:
    input:
        os.path.join(outdir, '{sample}.' + hostname + '.bam')
    output:
        os.path.join(outdir, PATTERN_R1 + "_" + hostname + '.mapped.fastq')
    conda:
        "envs/bowtie.yaml"
    shell:
        "samtools fastq -G 12 -f 65 {input} > {output}"

rule R2_reads_map_to_ref:
    input:
        os.path.join(outdir, '{sample}.' + hostname + '.bam')
    output:
        os.path.join(outdir, PATTERN_R2 + "_" + hostname + '.mapped.fastq')
    conda:
        "envs/bowtie.yaml"
    shell:
        "samtools fastq  -G 12 -f 129 {input} > {output}"

rule single_reads_map_to_ref:
    input:
        os.path.join(outdir, '{sample}.' + hostname + '.bam')
    output:
        os.path.join(outdir, '{sample}_singletons_' + hostname + '.mapped.fastq')
    conda:
        "envs/bowtie.yaml"
    shell:
        "samtools fastq  -F 5 {input} > {output}"

rule R1_unmapped:
    input:
        os.path.join(outdir, '{sample}.' + hostname + '.bam')
    output:
        os.path.join(outdir, PATTERN_R1 + "_" + hostname + '.unmapped.fastq')
    conda:
        "envs/bowtie.yaml"
    shell:
        "samtools fastq  -f 77  {input} > {output}"

rule R2_unmapped:
    input:
        os.path.join(outdir, '{sample}.' + hostname + '.bam')
    output:
        os.path.join(outdir, PATTERN_R2 + "_" + hostname + '.unmapped.fastq')
    conda:
        "envs/bowtie.yaml"
    shell:
        "samtools fastq  -f 141 {input} > {output}"

rule single_reads_unmapped:
    input:
        os.path.join(outdir, '{sample}.' + hostname + '.bam')
    output:
        os.path.join(outdir, '{sample}_singletons_' + hostname + '.unmapped.fastq')
    conda:
        "envs/bowtie.yaml"
    shell:
        "samtools fastq  -f 4 -F 1  {input} > {output}"
