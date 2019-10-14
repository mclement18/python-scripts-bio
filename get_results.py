#!/bin/bash/env python

import pandas as pd
import os

table = pd.DataFrame()

with open("strains.txt", "r") as strains:
    for strain in strains:
        mutfile = strain[:-1] + "_SNP2AAchange.out"
        mutfilepath = os.path.join(strain[:-1], mutfile)
        dic = {"Strain" : strain[:-1]}
        with open(mutfilepath, "r") as infile:
            for line in infile:
                if line.startswith("Gene"):
                    gene = line.split(" ")[1]
                    dic[gene] = []
                elif line.startswith("AA"):
                    dic[gene].append(line.split(" ")[2][:-1])
        table = table.append(dic, ignore_index=True)

for index, row in table.iterrows():
    for key, value in row.iteritems():
        if value == []:
            table.at[index, key] = "None"
        elif key != "Strain":
            table.at[index, key] = ", ".join(value)

table.set_index("Strain", inplace=True)

table.to_csv("SNP2AAmutations.tsv", sep="\t")
