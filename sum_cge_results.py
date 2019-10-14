#!/bin/bash/env python

import os
import pandas as pd
import json
import argparse
from glob import glob
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Summarise CGE results in a table")
parser.add_argument("-s", "--strains", help="File containing strain list")
parser.add_argument("-f", "--infolder", help="Folder containing strains folders")
parser.add_argument("-o", "--outfolder", help="Outfolder")
args = parser.parse_args()

def loop_strains(strain_list, folder):
    all_plasmids = pd.DataFrame(columns = ["Plasmids", "Strains", "Contigs", "Resgenes"])
    all_resgenes = pd.DataFrame(columns = ["Strains", "Resgenes", "Contigs", "Position", "Plasmid_type/Chromosome"])
    all_STs = pd.DataFrame(columns = ["Strains", "ST", "Alleles", "Notes"])
    all_pointmuts = pd.DataFrame(columns = ["Strains", "Mutations", "NucChange", "Res"])
    with open(strain_list, 'r') as strains:
        for strain in strains:
            print(strain)
            strain_path = os.path.join(folder, strain[:-1])
            chromo = get_chromo(strain[:-1], strain_path)
            plasmids = plasmid(strain[:-1], strain_path)
            resgenes, plasmids_w_res = resgene(strain[:-1], strain_path, plasmids, chromo)
            STs = ST(strain[:-1], strain_path)
            pointmut, present = PMR(strain[:-1], strain_path)
            all_plasmids = all_plasmids.append(plasmids_w_res, ignore_index = True)
            all_resgenes = all_resgenes.append(resgenes, ignore_index = True)
            all_STs = all_STs.append(STs, ignore_index = True)
            if not present == "No":
                all_pointmuts = all_pointmuts.append(pointmut, ignore_index = True)
    tables = {'palsmids_table': all_plasmids, 'resistancegenes_table': all_resgenes, 'ST_table': all_STs, 'pointmutation_table': all_pointmuts}

    return tables


def plasmid(strain, strain_path):
    table = os.path.join(strain_path, "plasmidfinder", "data.json")
    plas = pd.DataFrame(columns = ["Plasmids", "Strains", "Contigs", "Resgenes"])
    with open(table, 'r') as plas_list:
        data = json.load(plas_list)
    plasmid_hits = data["plasmidfinder"]["results"]["Enterobacteriaceae"]["enterobacteriaceae"]
    if not plasmid_hits == "No hit found":
        hit_list = {}
        for hit in plasmid_hits.keys():
            if not plasmid_hits[hit]["contig_name"] in hit_list.keys():
                hit_list[plasmid_hits[hit]["contig_name"]] = []
                hit_list[plasmid_hits[hit]["contig_name"]].append(plasmid_hits[hit]["plasmid"])
            else:
                hit_list[plasmid_hits[hit]["contig_name"]].append(plasmid_hits[hit]["plasmid"])
        for tig in hit_list.keys():
            p_type = "-".join(hit_list[tig])
            row = {"Plasmids": p_type, "Strains": strain, "Contigs": tig, "Resgenes": []}
            plas = plas.append(row, ignore_index = True)

        return plas
    else:
        return pd.DataFrame()


def resgene(strain, strain_path, plasmids, chromo):
    def resinplas(tig, gene, plasmids, chromo):
        ptype = []
        if tig == chromo:
            ptype.append("Chromosomal")
        elif not plasmids.empty:
            for plas in plasmids.itertuples():
                if plas.Contigs == tig:
                    ptype.append(plas.Plasmids)
                    plasmids.at[plas.Index, "Resgenes"].append(gene)
        return ptype, plasmids
    table = os.path.join(strain_path, "resfinder", "results_tab.txt")
    res = pd.DataFrame(columns = ["Strains", "Resgenes", "Contigs", "Position", "Plasmid_type/Chromosome"])
    with open(table,'r') as res_list:
        hit_list = {}
        for hit in res_list:
            if not hit.startswith("Resistance"):
                gene = hit.split("\t")[0]
                tig = hit.split("\t")[5]
                pos = hit.split("\t")[6]
                if tig not in hit_list.keys():
                    hit_list[tig] = []
                    hit_list[tig].append(pos)
                    ptype, plasmids = resinplas(tig, gene, plasmids, chromo)
                    row = {"Strains": strain, "Resgenes": gene, "Contigs": tig, "Position": pos, "Plasmid_type/Chromosome": ptype}
                    res = res.append(row, ignore_index = True)
                elif pos not in hit_list[tig]:
                    hit_list[tig].append(pos)
                    ptype, plasmids = resinplas(tig, gene, plasmids, chromo)
                    row = {"Strains": strain, "Resgenes": gene, "Contigs": tig, "Position": pos, "Plasmid_type/Chromosome": ptype}
                    res = res.append(row, ignore_index = True)

    return res, plasmids


def ST(strain, strain_path):
    table = os.path.join(strain_path, "mlst", "data.json")
    st_table = pd.DataFrame(columns = ["Strains", "ST", "Alleles", "Notes"])
    with open(table, 'r') as data:
        mlst = json.load(data)
    st = mlst['MLST']['results']['sequence_type']
    notes = ''
    if st  == 'Unknown':
        st = mlst['MLST']['results']['nearest_sts']
        notes = mlst['MLST']['results']['notes']
    alleles = []
    for allele, info in mlst['MLST']['results']['allele_profile'].items():
        allele_nb = allele + info['allele_name']
        alleles.append(allele_nb)
    row = {"Strains": strain, "ST": st, "Alleles": alleles, "Notes": notes}
    st_table = st_table.append(row, ignore_index = True)

    return st_table



def PMR(strain, strain_path):
    folder = os.path.join(strain_path, "pointfinder", "*_blastn_results.tsv")
    try:
        table = glob(folder)[0]
    except IndexError:
        print("%s has no pointfinder data!" %(strain), flush=True)
        return None, "No"
    else:
        point = pd.DataFrame(columns = ["Strains", "Mutations", "NucChange", "Res"])
        with open(table, 'r') as point_list:
            for line in point_list:
                mut = line.split("\t")[0]
                NC = line.split("\t")[1]
                res = line.split("\t")[3]
                row = {"Strains": strain, "Mutations": mut, "NucChange": NC, "Res": res}
                point = point.append(row, ignore_index = True)

        return point, "Yes"


def get_chromo(strain, strain_path):
    genome_name = strain + ".final.fasta"
    genome_path = os.path.join(strain_path, genome_name)
    length = 0
    for record in SeqIO.parse(genome_path, "fasta"):
        if length < len(record.seq):
            length = len(record.seq)
            chromo = record.id

    return chromo


def write_out(outfolder, tables):
    for name, table in tables.items():
        table = clean_table(table, name)
        table_name = name + ".tsv"
        outtable = os.path.join(outfolder, table_name)
        table.to_csv(outtable, index = False, sep = "\t")


def clean_table(table, name):
    if name == 'palsmids_table':
        col = "Resgenes"
    elif name == 'resistancegenes_table':
        col = "Plasmid_type/Chromosome"
    elif name == 'ST_table':
        col = "Alleles"
    else:
        col = "no"
    if not col == "no":
        for index, row in table.iterrows():
            if row[col] == []:
                table.at[index, col] = "None"
            else:
                table.at[index, col] = ", ".join(table.at[index, col])

    return table


def main(args):
    tables = loop_strains(args.strains, args.infolder)
    write_out(args.outfolder, tables)


main(args)
