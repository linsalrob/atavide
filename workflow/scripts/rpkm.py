"""
Calculate the RPKM from the contig_hits file. 


We need:
    1. the length of the contigs
    2. the R1 sequences to get the number of sequences
    3. the contig_hits file

"""

import os
import sys
from atavide_lib import stream_fastq

__author__ = 'Rob Edwards'

if not os.path.exists(snakemake.input.r1):
    sys.stderr.write(f"ERROR: snakemake.input.r1 {snakemake.input.r1} not found\n")
if not os.path.exists(snakemake.input.contig_len):
    sys.stderr.write(f"ERROR: snakemake.input.contig_len {snakemake.input.contig_len} not found\n")
if not os.path.exists(snakemake.input.hits):
    sys.stderr.write(f"ERROR: snakemake.input.hits {snakemake.input.hits} not found\n")


numseqs = 0
for seqid, header, seq, qualscores in stream_fastq(snakemake.input.r1):
    numseqs+=1
numseqs /= 1e6

contigs = {}
with open(snakemake.input.contig_len, 'r') as f:
    for l in f:
        r = l.strip()
        if not r:
            continue
        p = r.split("\t")
        if len(p) < 2:
            sys.stderr.write(f"Skipped {r} in {f}\n")
            continue
        contigs[p[0]] = int(p[1])

with open(snakemake.input.hits, 'r') as f, open(snakemake.output.rpkm, 'w') as out:
    for l in f:
        p = l.strip().split("\t")
        if p[0] == '*':
            continue
        if p[0] not in contigs:
            sys.stderr.write(f"ERROR: contig {p[0]} does not have a length.\n")
            contigs[p[0]]=1
        rpkm = int(p[1])/numseqs/contigs[p[0]]
        print(f"{p[0]}\t{rpkm}", file=out)

