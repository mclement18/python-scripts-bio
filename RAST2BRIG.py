#!/bin/bash/env python

import os


tsv = "7102_58-plasmid-DHA-annotation.tsv"

infile = open(tsv, 'r')

forward = []
reverse = []
subsystem = []

for line in infile:
    if line.startswith("F"):
        header = line
    else:
        if line.split("\t")[6] == "-":
            reverse.append(line)
        else:
            forward.append(line)
        if line.split("\t")[9] not in subsystem:
            subsystem.append(line.split("\t")[9])

def splitbysubsys(L, sub, sens):
    reformL = []
    if sens == "forward":
        deco = "clockwise-arrow\n"
    elif sens == "reverse":
        deco = "counterclockwise-arrow\n"
    for item in L:
        elements = item.split("\t")
        if int(elements[3]) - int(elements[4]) < 0:
            nelements = [elements[3], elements[4], elements[8], deco, elements[9]]
        else:
            nelements = [elements[4], elements[3], elements[8], deco, elements[9]]
        nitem = "\t".join(nelements)
        reformL.append(nitem)

    classL = {}
    s2c = {}
    n = 0
    c = ["gray", "fuchsia", "teal", "purple", "orange", "olive", "silver", "navy", "lime", "marron", "yellow", "aqua", "green"]
    for s in sub:
        classL[s] = []
        s2c[s] = c[n]
        n += 1
    for item in reformL:
        st = item.split("\t")[4]
        nelements = item.split("\t")[:-1]
        nitem = "\t".join([nelements[0], nelements[1], nelements[2], s2c[st], nelements[3]])
        classL[st].append(nitem)

    return classL

reforward = splitbysubsys(forward, subsystem, "forward")
rereverse = splitbysubsys(reverse, subsystem, "reverse")

def writeout(dic, sens, sub):
    for s in sub:
        outf = s.replace(" ", "") + sens + ".tsv"
        with open(outf, 'w') as outfile:
            for item in dic[s]:
                outfile.write(item)

writeout(reforward, "forward", subsystem)
writeout(rereverse, "reverse", subsystem)
