#!bin/bash python

import os
import subprocess
import argparse

parser = argparse.ArgumentParser(description="Extract each contig from an assembly fasta file.")
parser.add_argument("-s", "--strain", help="strain name to be analysed")
parser.add_argument("-f", "--folder", help="folder containing the assembly -- use {Canu, Circlator, Pilon} in combination with '--strain' for precomputed folder path and assembly name")
parser.add_argument("-a", "--assembly", help="Name of the assembly")
args = parser.parse_args()

def contig_extract(assembly, folder):
	print("Folder used: " + folder, flush=True)
	print("Assembly used: " + assembly, flush=True)
	
	contigs_names = subprocess.getoutput("grep '^>' " + folder + assembly)
	contigs_full_list = contigs_names.split("\n")
	contigs_list = []
	
	for contig_long in contigs_full_list:
		contig_short = contig_long.split(" ")[0][1:]
		contigs_list.append(contig_short)
		os.system("samtools faidx " + folder + assembly + " " + contig_short + " > " + folder + contig_short + ".fa") 


if args.folder == 'Canu':
	folder_path = "~/1.Assembly/Canu/" + args.strain + "/"
	assembly_name = args.strain + ".contigs.fasta"
	contig_extract(assembly_name, folder_path)
	
elif args.folder == 'Circlator':
	folder_path = "~/2.Circlator/" + args.strain + "/"
	assembly_name = "06.fixstart.fasta"
	contig_extract(assembly_name, folder_path)
	
elif args.folder == 'Pilon':
	folder_path = "~/3.Polishing/" + args.strain + "/"
	assembly_name = args.strain + ".final.fasta"
	contig_extract(assembly_name, folder_path)
else:
	if not args.folder.endswith("/"):
		args.folder = args.folder + "/"
	contig_extract(args.assembly, args.folder)
