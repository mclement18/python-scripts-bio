#!/bin/bash/python

import json
import pandas as pd
import subprocess as sb

def perform_blast(query, subject):
    blast_command = ['blastn', 
                     '-subject', subject, 
                     '-outfmt', '15']
    fasta = query.encode('utf-8')
    blast_out = sb.run(blast_command, input=fasta, stdout=sb.PIPE)
    blast_stdout = blast_out.stdout.decode('utf-8')
    return blast_stdout


def read_json(blast_out):
    blast_json = json.loads(blast_out)
    if len(blast_json['BlastOutput2'][0]['report']['results']['bl2seq'][0]['hits']) == 0:
        blast_alignment = 'No hits found'
    else:
        blast_alignment = blast_json['BlastOutput2'][0]['report']['results']['bl2seq'][0]['hits'][0]['hsps'][0]
    return blast_alignment


def get_snp(alignment, alternative):
    missing_snp = 0
    if alignment == 'No hits found':
        snp = '-'
        missing_snp = 1
    else:
        snp_position = alignment['midline'].find(' ')
        if snp_position == -1:
            snp = alternative
        else:
            snp = alignment['hseq'][snp_position]  
    return snp, missing_snp


def get_genome_path(strain):
    strain_folder = 'genomes/' + strain
    ls_command = ['ls', strain_folder]
    ls_out = sb.run(ls_command, stdout=sb.PIPE)
    ls = ls_out.stdout.decode('utf-8').split('\n')
    for line in ls:
        if line.endswith('.fna'):
            genome = line
    genome_path = "/".join([strain_folder, genome])
    return genome_path


def create_fasta_query(ID, sequence):
    header = "".join(['>', ID])
    fasta = "\n".join([header, sequence])
    return fasta


def loop_through_strains(strain_list, bait_list):
    sequence_list = []
    for strain in strain_list:
        if "/" in strain:
            strain = strain.replace("/", "_")
        genome_path = get_genome_path(strain)
        sequence, missing_bait, nb_missing_snp = loop_through_baits(genome_path, bait_list)
        sequence_list.append(sequence)
        if nb_missing_snp > 0:
            print('Could not retrieve %d SNPs from strain %s!' % (nb_missing_snp, strain), end='\n\n', flush=True)
            print('A gap was inserted instead of each missing SNP!', end='\n\n', flush=True)
            print('List of baits not found:', flush=True)
            print(missing_bait, sep='\n', flush=True)
            print('', flush=True)
        print('SNPs sequence of strain %s was successfully created!' % (strain), end='\n\n\n', flush=True)
    return sequence_list


def loop_through_baits(genome, bait_list):
    snp_list = []
    missing_bait = []
    nb_missing_snp = 0
    for bait in bait_list.itertuples():
        fasta = create_fasta_query(bait.ID, bait.Sequence)
        blast_output = perform_blast(fasta, genome)
        alignment = read_json(blast_output)
        snp, missing_snp = get_snp(alignment, bait.Ss046_Nucl)
        snp_list.append(snp)
        if missing_snp == 1:
            nb_missing_snp = nb_missing_snp + missing_snp
            missing_bait.append(bait.ID)
    snp_sequence = "".join(snp_list)
    return snp_sequence, missing_bait, nb_missing_snp


def get_data():
    strain_list = []
    with open('genomes/strain_list.txt', 'r') as file:
        for line in file:
            add_to_list = line[:-1]
            strain_list.append(add_to_list)
    bait_list = pd.read_table('table1.txt')
    print('Data uploaded!', end='\n\n\n', flush=True)
    return strain_list, bait_list


def save_output(strain_list, sequence_list):
    out_table = pd.DataFrame({'Strains': strain_list, '97snpSequence': sequence_list})
    out_table.to_csv('snp_sequence_table.txt', sep='\t', index=False)
    print('SNPs sequences table saved!', end='\n\n', flush=True)


def main():
    strain_list, bait_list = get_data()
    sequence_list = loop_through_strains(strain_list, bait_list)
    save_output(strain_list, sequence_list)
    print('Done!', flush=True)


main()
