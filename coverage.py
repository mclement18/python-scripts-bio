#!/bin/bash/python

import os
import argparse
import string

parser = argparse.ArgumentParser(description="Plot the coverage of a BAM file.")
parser.add_argument("-i", "--BAM" ,help="Input BAM file", required=True)
parser.add_argument("-o", "--folder", help="Output folder", required=True)
parser.add_argument("-p", "--prefix", help="Output prefix", required=True)
args = parser.parse_args()

if not args.folder.endswith("/"):
	args.folder = args.folder + "/"

os.system("samtools depth " + args.BAM + " > " + args.folder + args.prefix + ".coverage")

def module(command, *arguments):
	commands = os.popen('/software.el7/lmod/lmod/libexec/lmod python %s %s' % (command, " ".join(arguments))).read()
	exec(commands)

module('load', 'R')

os.system("Rscript ~/R/scripts/covplot.R " + args.folder + args.prefix + ".coverage " + args.folder + " " + args.prefix)
