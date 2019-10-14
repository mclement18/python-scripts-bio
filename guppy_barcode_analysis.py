#!/bin/bash/env python3

import os
import pandas as pd
import argparse
from glob import glob
from Bio.SeqIO.QualityIO import FastqGeneralIterator

parser = argparse.ArgumentParser(description="Get stats from guppy_barcoder binning, merge multiple fastq files outputed by guppy_basecaller and split fastq into bins.")
parser.add_argument("-s", "--summary", help="barcoding_summary.txt from guppy_barcoder.")
parser.add_argument("-o", "--outdir", help="Output directory.")
parser.add_argument("--stats", help="Get stats from the demultiplexing (must give '--barcodes').", action='store_true')
parser.add_argument("--merge", help="Merge fastq files outputed by all guppy_bascaller threads.", action='store_true')
parser.add_argument("--split", help="Split fastq into barcode bins.", action='store_true')
parser.add_argument("-b", "--barcodes", nargs="*", help="Space-separated list of barcodes used for the run (e.g. 'barcode01 barcode02').")
parser.add_argument("-f", "--fastq", help="FASTQ file to split if not using merge.")
parser.add_argument("-v", "--verbose", help="Print stdout from commands.", action='store_true')
args = parser.parse_args()

def load_data(data):
    df = pd.read_table(data, sep="\t")
    print("barcode_summary.txt loaded!", end='\n\n', flush=True)
    return df


def binning(data):
    barcode2read = {}
    stats = {}
    stats_rear = {}
    for read in data.itertuples():
        barcode2read[read.read_id] = read.barcode_arrangement
        if read.barcode_arrangement in stats:
            stats[read.barcode_arrangement] += 1 
        else:
            stats[read.barcode_arrangement] = 1
        if not read.barcode_rear_id == '[none]':
            if read.barcode_arrangement in stats:
                stats_rear[read.barcode_arrangement] += 1
            else:
                stats_rear[read.barcode_arrangement] = 1

    print("Binning terminated!", end='\n\n', flush=True)

    return barcode2read, stats, stats_rear



def make_total(stats, stats_rear, true_barcode):
    tot_nb_reads = 0
    tot_nb_reads_with_rear = 0
    true_nb_reads = 0
    true_nb_reads_with_rear = 0

    for name, value in stats.items():
        print("%s contains %d reads." %(name, value), flush=True)
        tot_nb_reads += value
        ex = ''
        try:
            print("%s contains %d reads with rear barcode." %(name, stats_rear[name]), end='\n\n', flush=True)
            tot_nb_reads_with_rear += stats_rear[name]
        except KeyError as e:
            ex = e
            print("%s does not contain reads with rear barcode." %(name), end='\n\n', flush=True)
        if name in true_barcode:
            print("%s was used in the sequencing run." %(name), end='\n\n', flush=True)
            true_nb_reads += value
            if ex == '':
                true_nb_reads_with_rear += stats_rear[name]

    print("Total number of reads: %d" %(tot_nb_reads), flush=True)
    print("Total number of reads with a rear barcode: %d" %(tot_nb_reads_with_rear), end='\n\n', flush=True)
    print("Total number of reads from barcodes used in sequening run: %d" %(true_nb_reads), flush=True)
    print("Total number of reads from barcodes used in sequening run with a rear barcode: %d" %(true_nb_reads_with_rear), end='\n\n', flush=True)



def merge_fastq(outdir):
    fastq_list = glob(outdir + "*.fastq")
    merged_fastq = fastq_list[0].rsplit("_", 1)[0] + ".fastq"
    
    with open(merged_fastq, 'a') as new:
        for fastq in fastq_list:
            with open(fastq, 'r') as fq:
                for line in fq:
                    new.write(line)
    print("FASTQ files merged!", end='\n\n', flush=True)

    return merged_fastq


def split_fastq_in_bin(outdir, barcodes2read, fastq):
    with open(fastq, 'r') as in_fq:
        for title, seq, qual in FastqGeneralIterator(in_fq):
            read_id = title.split(None, 1)[0]
            bin_name = barcodes2read[read_id] + ".fastq"
            current_binned_fastq = os.path.join(outdir, bin_name)
            with open(current_binned_fastq, 'a') as out_fq:
                out_fq.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))

    print("FASTQ splitted into bins!", end='\n\n', flush=True)


def main(args):
    if args.merge:
        fastq = merge_fastq(args.outdir)
    else:
        fastq = args.fastq
    if args.stats or args.split:
        data = load_data(args.summary)
        barcodes2read, stats, stats_rear = binning(data)
    if args.stats:
        make_total(stats, stats_rear, args.barcodes)
    if args.split:
        split_fastq_in_bin(args.outdir, barcodes2read, fastq)
    print("Done.", flush=True)


main(args)
