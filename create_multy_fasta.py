#!/bin/bash python

import pandas as pd

def load_data():
    df = pd.read_table('snp_sequence_table.txt')
    return df

def create_fasta(ID, sequence):
    header = "".join(['>', ID])
    fasta = "\n".join([header, sequence])
    return fasta

def loop_through_table(data):
    fasta_list = []
    for row in data.itertuples():
        fasta = create_fasta(row.Strains, row.SNP_Sequence)
        fasta_list.append(fasta)
    multy_fasta = "\n".join(fasta_list)
    return multy_fasta

def save_out(data):
    with open('snp_sequences.fasta', 'w') as file:
        file.write(data)

def main():
   data = load_data()
   fasta = loop_through_table(data)
   save_out(fasta)

main()
