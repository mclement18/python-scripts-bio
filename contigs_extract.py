#!bin/bash python

import os
import subprocess
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Extract each contig from an assembly fasta file.")
parser.add_argument("-s", "--strain", help="strain name to be analysed")
parser.add_argument("-f", "--folder", help="folder containing the assembly -- use {Canu, Circlator, Pilon} in combination with '--strain' for precomputed folder path and assembly name")
parser.add_argument("-a", "--assembly", help="Name of the assembly")
args = parser.parse_args()

def contig_extract(assembly, folder):
	print("Folder used: " + folder, flush=True)
	print("Assembly used: " + assembly, flush=True)

	contigs_list = os.path.join(folder, "contigs_list.txt")
	assembly_path = os.path.join(folder, assembly)

	with open(contigs_list, "w") as outfile:
		for contig in SeqIO.parse(assembly_path, "fasta"):
			contig_path = os.path.join(folder, contig.id + ".fasta")
			SeqIO.write(contig, contig_path, "fasta")
			outfile.write(contig.id + "/n")


if args.folder == 'Canu':
	folder_path = os.path.join("~/1.Assembly/Canu/", args.strain)
	assembly_name = args.strain + ".contigs.fasta"
	contig_extract(assembly_name, folder_path)

elif args.folder == 'Circlator':
	folder_path = os.path.join("~/2.Circlator/", args.strain)
	assembly_name = "06.fixstart.fasta"
	contig_extract(assembly_name, folder_path)

elif args.folder == 'Pilon':
	folder_path = os.path.join("~/3.Polishing/", args.strain)
	assembly_name = args.strain + ".final.fasta"
	contig_extract(assembly_name, folder_path)
else:
	contig_extract(args.assembly, args.folder)
