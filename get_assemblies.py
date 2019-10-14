#!/bin/python

import pandas as pd
import subprocess as sb
import re


def get_ftp_path(strain, df):
    for row in df.itertuples():
        if row.Strain == strain:
            letters, numbers = row.File.split('_')
            number_list = re.findall('...', numbers)
            path_1, path_2, path_3, path_4, path_5, path_6 = row.GenBank_FTP.split("/")
            path_7 = path_6 + "_genomic.fna.gz"
            list_to_join = [path_1, path_2, path_3, path_4, path_5, 'GCA'] + number_list + [path_6, path_7]
            new_ftp = "/".join(list_to_join)
            return new_ftp


def get_assembly(strain, ftp):
    strain_folder = "genomes/" + strain
    wget_command = ['wget', '-P', strain_folder, ftp]
    sb.run(wget_command)


def loop_through_strains(strain_list, df):
    for strain in strain_list:
        if "/" in strain:
            continue
        ftp = get_ftp_path(strain, df)
        get_assembly(strain, ftp)


def get_strains_and_table():
    df = pd.read_table('newtable2.txt')
    strain_list = []
    with open('genomes/strain_list.txt', 'r') as file:
        for line in file:
            current_line = line[:-1]
            strain_list.append(current_line)
    return df, strain_list 


def main():
    df, strain_list = get_strains_and_table()
    loop_through_strains(strain_list, df)

main()
