#!/bin/bash/env python

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(description="Extract one or two parts of a fasta sequence and concatanate them into a new fasta (<fasta>_new<.fasta>)")
parser.add_argument("contig", help="contig")
parser.add_argument('start_end', type=int, nargs=2, help="start and end of sequence")
parser.add_argument("--start_end2", type=int, nargs=2, help="start and end of 2nd part of the sequence")
args = parser.parse_args()
print(args)

start, end = args.start_end

contig = SeqIO.read(args.contig, 'fasta')
seq1 = contig.seq[start - 1: end]

if args.start_end2 != None:
    start2, end2 = args.start_end2
    seq2 = contig.seq[start2 - 1: end2]
    newseq = seq1 + seq2
else:
    newseq = seq1

newrec = SeqRecord(newseq, id=contig.id, description='')

tigout = "_new.".join(args.contig.rsplit(".", 1))

SeqIO.write(newrec, tigout, 'fasta')
