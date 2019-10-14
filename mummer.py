#!bin/bash python

import subprocess as sb
import argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser


parser = argparse.ArgumentParser(description="Run nucmer analysis.")
parser.add_argument("-r", help="Reference")
parser.add_argument("-q", help="Query")
parser.add_argument("-p", help="Prefix")
parser.add_argument("-v", action="store_true", help="Retrieve stdout of commands")
args = parser.parse_args()
print(args)

def run_mummer(contig1, contig2, prefix, verbiosity):
    mummer_command = ["mummer", "-maxmatch", "-n", "-b", "-c",
                      contig1, contig2]
    mumout = prefix + ".mums"

    with open(mumout, 'w') as out:
        mummer = sb.Popen(mummer_command, stdout=out, stderr=sb.PIPE)
        mummer_out, mummer_err = mummer.communicate()

    if verbiosity:
        print(mummer_out, mummer_err, flush=True)

    return mumout

def run_mummerplot(mumout, prefix, contig, verbiosity):
    axes_range = "[0," + str(get_fasta_length(contig)) + "]"
    mumplot_command = ["mummerplot", "-x", axes_range, "-y", axes_range,
                       "--png", "-p", prefix, mumout]

    mumplot = sb.check_output(mumplot_command, stderr=sb.STDOUT)

    if verbiosity:
        print(mumplot, flush=True)

def get_fasta_length(fasta):
    with open(fasta, 'r') as infile:
        for name, seq in SimpleFastaParser(infile):
            seqlen = len(seq)
    return seqlen
    

def main(args):
    mumout = run_mummer(args.r, args.q, args.p, args.v)
    run_mummerplot(mumout, args.p, args.q, args.v)

main(args)
