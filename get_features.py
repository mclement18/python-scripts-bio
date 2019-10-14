#!/bin/bash/env python

import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-i", help="input")
parser.add_argument("-t", help="target")
parser.add_argument("-o", help="output")

args = parser.parse_args()

targets = []
with open(args.t, "r") as target:
    for line in target:
        targets.append(line[:-1])

print(targets)

infile = open(args.i, "r")
outfile = open(args.o, "w")

with open(args.o, "w") as outfile:
    with open(args.i, "r") as infile:
        for line in infile:
            name = line.split("\t")[9].split(";")[3].split("=")[1]  
            if name in targets:
                outfile.write(line)

