#!/usr/bin/env python3

"""
Count the characters in a fasta file. We summarize the longest and shortest reads
and the N50 of the data set.
"""

import os
import sys
from atavide_lib import stream_fastq

__author__ = 'Rob Edwards'


# snakemake.input is the directory of fastq files to count
# snakemake.output.stats is the output file to write

with open(snakemake.output.stats, 'w') as out:
    print("File\tNumber of Sequences\tTotal bp\tShortest length\tLongest length\tN50\tN75\tAuN", file=out)
    for f in os.listdir(snakemake.input):
        if 'fastq' in f or 'fq' in f:
            lens = []
            for (sid, label, seq, qual) in stream_fastq(os.path.join(snakemake.input, f)):
                lens.append(len(seq))
            lens.sort()
            length=sum(lens)

            len_so_far = 0
            n50 = None
            n75 = None
            auN = 0
            for i in lens:
                len_so_far += i
                if not n50 and len_so_far >= length * 0.5:
                    n50 = i
                if not n75 and len_so_far >= length * 0.75:
                    n75 = i
                auN += i**2

            auN /= length

            print(f"{f}\t{len(lens):,}\t{length:,}\t{lens[0]:,}\t" \
                  + f"{lens[-1]:,}\t{n50:,}\t{n75:,}\t{int(auN):,}", file=out)
        else:
            sys.stderr.write(f"Skipped {os.path.join(subdir, f)}. Does not appear to be fastq\n")

