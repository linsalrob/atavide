#!/usr/bin/env python3

"""
Renumber multiple fasta files, starting with either number 1 or a number provided on the command line.
At the end, we output the last number written.

Written for snakemake so uses those objects not the command line
"""

from atavide_lib import read_fasta, colours

__author__ = 'Rob Edwards'

# snakemake.input is an array of files to merge
# snakemake.output.contigs is the name of the final contigs file
# snakemake.output.ids is the name of the ids file to write3
# snakemake.params.sample_id is the sample_id that will be prepended to the contig names


idmap = open(snakemake.output.ids, 'w')
out = open(snakemake.output.contigs, 'w')

counter = 0

for f in snakemake.input:
    fa = read_fasta(f)
    for id in fa:
        counter += 1
        if hasattr(snakemake, 'params') and hasattr(snakemake.params, 'sample_id'):
            out.write(f">{snakemake.params.sample_id}_{counter}\n{fa[id]}\n")
            idmap.write(f"{f}\t{id}\t{snakemake.params.sample_id}_{counter}\n")
        else:
            out.write(f">{counter}\n{fa[id]}\n")
            idmap.write(f"{f}\t{id}\t{counter}\n")

idmap.close()
out.close()
