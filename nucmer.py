#!bin/bash python

import subprocess as sb
import argparse


parser = argparse.ArgumentParser(description="Run nucmer analysis.")
parser.add_argument("-r", help="Reference")
parser.add_argument("-q", help="Query")
parser.add_argument("-p", help="Prefix")
parser.add_argument("-v", action="store_true", help="Retrieve stdout of commands")
args = parser.parse_args()
print(args)

def run_nucmer(contig1, contig2, prefix, verbiosity):
    nucmer_command = ["nucmer", "--maxmatch", "--nosimplify",
                         "-p", prefix, contig1, contig2]

    nucmer = sb.Popen(nucmer_command, stdout=sb.PIPE, stderr=sb.PIPE)
    nucmer_out, nucmer_err = nucmer.communicate()

    if verbiosity:
        print(nucmer_out, nucmer_err, flush=True)

def run_showcoord(prefix, verbiosity):
    delta = prefix + ".delta"
    outfile = prefix + ".nucmer"
    showcoord_command = ["show-coords", "-lrcT", delta]

    with open(outfile, 'w') as out:
        showcoord = sb.Popen(showcoord_command, stdout=out, stderr=sb.PIPE)
        showcoord_out, showcoord_err = showcoord.communicate()

    if verbiosity:
        print(showcoord_out, showcoord_err, flush=True)

def main(args):
    run_nucmer(args.r, args.q, args.p, args.v)
    run_showcoord(args.p, args.v)

main(args)
