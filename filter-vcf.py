#!/bin/bash/env python

import argparse
import os
import sys
import subprocess as sb
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML

parser = argparse.ArgumentParser(description="Convert SNP to AA changes using VCF and target genes annotation in BED.") 
parser.add_argument("-i", "--vcf", help="vcf input file")
parser.add_argument("-r", "--ref", help="annotation ref file")
args = parser.parse_args()

def run_bedtools_inersect(vcf, ref):
    with open("temp-ref", "w") as infile:
        infile.write(ref)

    intersect_command = ["bedtools", "intersect",
                         "-a", vcf, "-b", "stdin", "-wa"]
    with open("temp-ref", "r") as infile:
        intersect_out = sb.check_output(intersect_command, stdin=infile)

    gene_sequence = []
    gene_snps = []
    nb_snp = 0
    for line in intersect_out.decode().split("\n")[:-1]:
        gene_sequence.append(line.split("\t")[3])
        if not line.startswith("#") and line.split("\t")[4] != "." and line.split("\t")[6] == "PASS":
            nb_snp += 1
            gene_snps.append(line)
    return gene_sequence, gene_snps, nb_snp

def run_blastp(sequence, snp_sequence, strand, snp_pos):
    if len(sequence) != len(snp_sequence):
        sys.exit("INDEL!\nCheck manually!")
    gene_sequence = Seq(sequence, IUPAC.unambiguous_dna)
    allele_sequence = Seq(snp_sequence, IUPAC.unambiguous_dna)
    if strand == "-":
        gene_prot = gene_sequence.reverse_complement().translate(table=11)
        snp_prot = allele_sequence.reverse_complement().translate(table=11)
        snp_pos = len(sequence) - snp_pos
    else:
        gene_prot = gene_sequence.translate(table=11)
        snp_prot = allele_sequence.translate(table=11)
        snp_pos = snp_pos + 1

    temp_query = "query.fasta"
    temp_subject = "subject.fasta"

    rec_gene = SeqRecord(gene_prot, id="Gene_prot")
    rec_snp = SeqRecord(snp_prot, id="snp_prot")

    SeqIO.write(rec_gene, temp_subject, "fasta")
    SeqIO.write(rec_snp, temp_query, "fasta")

    blastp_cline = NcbiblastpCommandline(query=temp_query, subject=temp_subject, out="-", outfmt=5)
    blastp_out = sb.Popen(str(blastp_cline), stdout=sb.PIPE, stderr=sb.PIPE, universal_newlines=True, shell=(sys.platform!="win32"))
    blast_record = NCBIXML.read(blastp_out.stdout)
    blastp_out.stdout.close()
    blastp_out.stderr.close()
    hsp = blast_record.alignments[0].hsps[0]
    if hsp.match.find(" ") != (-1):
        mut_pos = hsp.match.find(" ")
        mut_type = "non-synonymous"
    elif hsp.match.find("+") != (-1):
        mut_pos = hsp.match.find("+")
        mut_type = "non-synonymous but similar AA properties"
    else:
        mut_type = "synonymous"
        if snp_pos % 3 == 0:
            mut_pos = (snp_pos // 3) - 1
        else:
            mut_pos = snp_pos // 3
    wt_aa = hsp.sbjct[mut_pos]
    mut_aa = hsp.query[mut_pos]

    return mut_pos, wt_aa, mut_aa, mut_type

def loop_through_snps(gene, gene_sequence, snps, nb_snp):
    sequence = "".join(gene_sequence)
    for snp in snps:
        snp_pos = int(snp.split("\t")[1]) - int(gene.split("\t")[1]) - 1
        snp_sequence_list = list(gene_sequence)
        snp_sequence_list[snp_pos] = snp.split("\t")[4] 
        snp_sequence = "".join(snp_sequence_list)
        mut_pos, wt_aa, mut_aa, mut_type = run_blastp(sequence, snp_sequence, gene.split("\t")[5], snp_pos)
        print("NT old: %s NT new: %s Pos_genome: %d Pos_gene: %d" %(gene_sequence[snp_pos], snp_sequence_list[snp_pos], int(snp.split("\t")[1]), (snp_pos + 1)), flush=True)
        print("Mututation type is: %s." %(mut_type), flush=True)
        print("AA change: %s%d%s" %(wt_aa, (mut_pos + 1), mut_aa), flush=True)
        print("\n", flush=True)

def loop_through_genes(genes, vcf):
    with open(genes, "r") as infile:
        for line in infile:
            gene_sequence, snps, nb_snp = run_bedtools_inersect(vcf, line)
            #print("Gene %s contains %d SNPs." %(line.split("\t")[9].split(";")[3].split("=")[1], nb_snp), flush=True)
            print("Gene %s contains %d SNPs." %(line.split("\t")[9].split("|")[2][:-1], nb_snp), flush=True)
            loop_through_snps(line, gene_sequence, snps, nb_snp)

def main(args):
    loop_through_genes(args.ref, args.vcf)

main(args)
