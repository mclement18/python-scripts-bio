#!bin/bash python

import os
import subprocess as sb
import argparse
import sys
from glob import glob
from shutil import copy2

parser = argparse.ArgumentParser(description="Make several round of Pilon polishing. First with '--fix bases' until no change, then with '--fix all' until no change. Do max 15 rounds.")
parser.add_argument("-a", "--assembly", help="Assembly to polish.")
parser.add_argument("-f", "--forward", help="Forward illumina reads.")
parser.add_argument("-r", "--reverse", help="Reverse illumina reads.")
parser.add_argument("-o", "--outdir", help="Output directory.")
parser.add_argument("-p", "--prefix", help="Output prefix.")
parser.add_argument("--fix", choices=['bases', 'all'], default='bases', metavar='fix', help="Define fix type to start with [bases, all]. (default: 'bases')")
parser.add_argument("--limit", type=int, default=15, help="Set the maximal number of polishing rounds. (default: 15)")
parser.add_argument("--restart", help="Restart the polishing from last finished step.", action='store_true')
parser.add_argument("--threads", type=int, default=1, help="Set the number of threads for bowtie2 alignment. (default: 1)")
parser.add_argument("-v", "--verbose", help="Print stdout from commands.", action='store_true')
args = parser.parse_args()

def pilon_polish(assembly, bam, prefix, outdir, fix_type, verbose):
    pilon_command = ['pilon',
                     '--genome', assembly, '--bam', bam,
                     '--output', prefix, '--outdir', outdir,
                     '--fix', fix_type, '--changes']

    print_command(pilon_command)

    pilon = sb.check_output(pilon_command, stderr=sb.STDOUT)

    fasta_name = prefix + ".fasta"
    pilon_fasta = os.path.join(outdir, fasta_name)
    pilon_newfasta = pilon_fasta + ".new"
    with open(pilon_fasta, 'r') as fasta:
        with open(pilon_newfasta, 'w') as newfasta:
            for line in fasta:
                newfasta.write(line.replace("_pilon", ""))

    clean_up(pilon_fasta)
    os.rename(pilon_newfasta, pilon_fasta)

    changes_name = prefix + ".changes"
    pilon_changes = os.path.join(outdir, changes_name)
    changes = 0
    with open(pilon_changes, 'r') as changes_file:
        for change in changes_file:
            changes += 1

    print("Number of pilon changes: %d" %(changes), flush=True)

    if verbose:
        print_out(pilon)

    return pilon_fasta, changes, pilon_changes

def bowtie_align(assembly, one, two, outdir, n, status, threads, verbose):
    def fasta_index(assembly, outdir):
        index = os.path.join(outdir, "index")
        build_index_command = ['bowtie2-build', assembly, index]

        print_command(build_index_command)
        build_index = sb.check_output(build_index_command, stderr=sb.STDOUT)

        genome_index = glob(index + "*")

        return index, genome_index, build_index

    def alignment(index, outdir, one, two, threads, n):
        sam_name = n + "_reads.aligned.sam"
        sam = os.path.join(outdir, sam_name)

        bowtie_command = ["bowtie2", "-x", index,
                          "-1", one, "-2", two, "-p", threads,
                          "--local", "--very-sensitive-local",
                          "-S", sam]
        command2print = bowtie_command
        command2print[8] = str(command2print[8])

        print_command(command2print)
        bowtie = sb.check_output(bowtie_command, stderr=sb.STDOUT)

        return sam, bowtie

    def bamconv(sam, outdir, n):
        bam_name = n + "_reads.aligned.bam"
        bam = os.path.join(outdir, bam_name)

        samtools_view_command = ['samtools', 'view', '-bS', sam]

        print_command(samtools_view_command)
        with open(bam, 'w') as bam_file:
            view = sb.Popen(samtools_view_command, stdout=bam_file, stderr=sb.PIPE)
            view_out, view_err = view.communicate()

        return bam, view_out, view_err

    def bamsort(bam, outdir, n):
        sortbam_name = n + "_reads.aligned.sorted.bam"
        sortbam = os.path.join(outdir, sortbam_name)

        samtoolssort_command = ['samtools', 'sort', bam, '-o', sortbam]

        print_command(samtoolssort_command)
        sort = sb.check_output(samtoolssort_command, stderr=sb.STDOUT)

        return sortbam, sort

    def indexbam(sortbam):
        bam_index = sortbam + ".bai"
        samtoolsindex_command = ['samtools', 'index', sortbam]

        print_command(samtoolsindex_command)
        bamindex = sb.check_output(samtoolsindex_command, stderr=sb.STDOUT)

        return bam_index, bamindex

    if status == None or status == 'fasta':
        index, genome_index, build_index = fasta_index(assembly, outdir)
        sam, bowtie = alignment(index, outdir, one, two, threads, n)
        bam, view_out, view_err = bamconv(sam, outdir, n)
        sortbam, sort = bamsort(bam, outdir, n)
        bam_index, bamindex = indexbam(sortbam)
        new_status = None
    elif status == 'sam':
        build_index = ""
        index = os.path.join(outdir, "index")
        genome_index = glob(index + "*")
        sam_name = n + "_reads.aligned.sam"
        sam = os.path.join(outdir, sam_name)

        bam, view_out, view_err = bamconv(sam, outdir, n)
        sortbam, sort = bamsort(bam, outdir, n)
        bam_index, bamindex = indexbam(sortbam)
        new_status = None
    elif status == 'bam':
        build_index, view_out,  view_err = ["", "", ""]
        index = os.path.join(outdir, "index")
        genome_index = glob(index + "*")
        sam_name = n + "_reads.aligned.sam"
        sam = os.path.join(outdir, sam_name)
        bam_name = n + "_reads.aligned.bam"
        bam = os.path.join(outdir, bam_name)

        sortbam, sort = bamsort(bam, outdir, n)
        bam_index, bamindex = indexbam(sortbam)
        new_status = None
    elif status == 'sortbam':
        build_index, view_out, view_err, sort = ["", "", "", ""]
        index = os.path.join(outdir, "index")
        genome_index = glob(index + "*")
        sam_name = n + "_reads.aligned.sam"
        sam = os.path.join(outdir, sam_name)
        bam_name = n + "_reads.aligned.bam"
        bam = os.path.join(outdir, bam_name)
        sortbam_name = n + "_reads.aligned.sorted.bam"
        sortbam = os.path.join(outdir, sortbam_name)

        bam_index, bamindex = indexbam(sortbam)
        new_status = None

    if verbose:
        print_out(build_index, bowtie, view_out, view_err, sort, bamindex)

    return sam, bam, sortbam, bam_index, genome_index, new_status

def print_out(*variables):
    print(variables, flush=True)

def print_command(command):
    print("Command:", flush=True)
    print(" ".join(command), end='\n\n', flush=True)

def loop(assembly, one, two, prefix, outdir, fix, limit, restart, threads, verbose):
    if restart:
        status, n, new_assembly = check_progress(outdir)
        print("SerialPilon restarted at round %d with status %s, argument '%s' and a round limit of %d." %(n, status, fix, limit), end='\n\n', flush=True)
    else:
        print("SerialPilon started with argument '%s' and a round limit of %d." %(fix, limit), end='\n\n', flush=True)
        status = None
        n = 1
        moved_name = "0_" + prefix + ".fasta"
        moved_assembly = os.path.join(outdir, moved_name)
        copy2(assembly, moved_assembly)

    round_limit = limit + 1
    fix_type = fix
    basedone = False

    while True:
        if n == 1:
            current_assembly = moved_assembly
            current_prefix = str(n) + "_" + prefix
        elif not basedone:
            current_assembly = new_assembly
            current_prefix = str(n) + "_" + prefix
        else:
            current_prefix = str(n) + "_" + prefix
            basedone = False

        print("Round %d" %(n), end='\n\n', flush=True)

        current_sam, current_bam, current_sortbam, bam_index, genome_index, status = bowtie_align(current_assembly, one, two, outdir, str(n), status, threads, verbose)

        new_assembly, changes, current_changes = pilon_polish(current_assembly, current_sortbam, current_prefix, outdir, fix_type, verbose)

        n += 1

        if n >= round_limit:
            print("Reached %d rounds. Polishing stops." %(n - 1), flush=True)
            break

        if changes == 0 and fix_type == "bases":
            fix_type = "all"
            basedone = True
            print("Fix 'bases' done. Change to 'all'.", end='\n\n', flush=True)
            clean_up(current_sam, current_bam, current_sortbam, bam_index, *genome_index, current_changes, new_assembly)
        elif changes == 0 and fix_type == "all":
            print("No more change possible.", end='\n\n', flush=True)
            clean_up(current_sam, current_bam, *genome_index, current_changes)
            final_name = prefix + ".final.fasta"
            final_assembly = os.path.join(outdir, final_name)
            os.rename(new_assembly, final_assembly)
            break
        else:
            print("Continue for a new round of pilon '%s'" %(fix_type), end="\n\n", flush=True)
            clean_up(current_sam, current_bam, current_sortbam, bam_index, *genome_index)

def clean_up(*args):
    for file2delete in args:
        if os.path.exists(file2delete):
            os.remove(file2delete)
        else:
            print("File %s does not exist!" %(file2delete), end='\n\n', flush=True)

def check_paths(path):
    if not os.path.exists(path):
        sys.exit("File %s does not exist!" %(path))


def check_progress(outdir):
    import re
    def sorted_nicely( l ):
        """ Sorts the given iterable in the way that is expected.
 
        Required arguments:
        l -- The iterable to be sorted.
 
        """
        convert = lambda text: int(text) if text.isdigit() else text
        alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
        return sorted(l, key = alphanum_key)

    sam = "*reads.aligned.sam"
    bam = "*reads.aligned.bam"
    sortbam = "*reads.aligned.sorted.bam"
    bai = "*reads.aligned.sorted.bam.bai"
    fasta = "*.fasta"

    sam_path = os.path.join(outdir, sam)
    bam_path = os.path.join(outdir, bam)
    sortbam_path = os.path.join(outdir, sortbam)
    bai_path = os.path.join(outdir, bai)
    fasta_path = os.path.join(outdir, fasta)

    current_sam = glob(sam_path)
    current_bam = glob(bam_path)
    current_sortbam = glob(sortbam_path)
    current_bai = glob(bai_path)
    current_fastas = glob(fasta_path)
    last_fasta = sorted_nicely(current_fastas)[-1]

    if current_bai != []:
        status = "sortbam"
        n = int(float(current_sortbam[0].split("/")[-1].split("_")[0]))
    elif current_sortbam != []:
        status = "bam"
        n = int(float(current_bam[0].split("/")[-1].split("_")[0]))
    elif current_bam != []:
        status = "sam"
        n = int(float(current_sam[0].split("/")[-1].split("_")[0]))
    elif current_sam != [] or current_sam == []:
        status = "fasta"
        n = int(float(last_fasta.split("/")[-1].split("_")[0])) + 1

    if int(float(last_fasta.split("/")[-1].split("_")[0])) == n:
        last_fasta = sorted_nicely(current_fastas)[-2]

    return status, n, last_fasta


def main(args):
    check_paths(args.assembly)
    check_paths(args.forward)
    check_paths(args.reverse)

    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir, mode=0o775)

    loop(args.assembly, args.forward, args.reverse, args.prefix, args.outdir, args.fix, args.limit, args.restart, args.threads, args.verbose)
    print("Done.", flush=True)

main(args)
