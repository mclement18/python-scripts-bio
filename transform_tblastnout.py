#!bin/bash/env python

import pandas as pd

inputfile = "tblastn_rpsJ.out"

outfile = "pseudo_annotation_rpsJ.gff"

dic = {"contig":[], "start":[], "end":[], "cds":[], "point":[], "strand":[], "Protein_homo":[], "CDS":[], "zero":[], "note":[]}

with open(inputfile, "r") as infile:
    for line in infile:
        elements = line.split("\t")
        dic["contig"].append(elements[1])
        dic["cds"].append("cds")
        dic["point"].append(".")
        dic["Protein_homo"].append("Protein_homo")
        dic["CDS"].append("CDS")
        dic["zero"].append(0)
        dic["note"].append(elements[0])
        if int(elements[9]) - int(elements[8]) > 0:
            dic["start"].append(int(elements[8]) - 1)
            dic["end"].append(elements[9])
            dic["strand"].append("+")
        else:
            dic["start"].append(int(elements[9]) - 1)
            dic["end"].append(elements[8])
            dic["strand"].append("-")

table = pd.DataFrame(dic)

table.to_csv(outfile, sep="\t", header=False, index=False)
