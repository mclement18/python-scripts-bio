#!bin/bash python

import os
import subprocess as sb
import argparse


parser = argparse.ArgumentParser(description="Run CGE softwares (ResFinder, PointFinder and PlasmidFinder)")
parser.add_argument("-a", "--assembly", help="Name of the assembly")
parser.add_argument("-o", "--folder", help="Folder containing the assembly")
parser.add_argument("--all", action="store_true", help="Run all programs")
parser.add_argument("--res", action="store_true", help="Run ResFinder")
parser.add_argument("--point", action="store_true", help="Run PointFinder")
parser.add_argument("--plasmid", action="store_true", help="Run PlasmidFinder")
parser.add_argument("--mlst", action="store_true", help="Run MLST")
parser.add_argument("--pmlst", action="store_true", help="Run pMLST")
parser.add_argument("-rc", "--rescoverage", default=0.6, help="Coverage threshold for ResFinder. [0.6]")
parser.add_argument("-ri", "--residentity", default=0.9, help="Identity threshold for ResFinder. [0.9]")
parser.add_argument("-pc", "--pointcoverage", default=0.6, help="Coverage threshold for PointFinder. [0.6]")
parser.add_argument("-pi", "--pointidentity", default=0.9, help="Identity threshold for PointFinder. [0.9]")
parser.add_argument("-mc", "--plasmcoverage", default=0.6, help="Coverage threshold for PlasmidFinder. [0.6]")
parser.add_argument("-mi", "--plasmidentity", default=0.9, help="Identity threshold for PlasmidFinder. [0.95]")
parser.add_argument("--specie", metavar='SPECIE', choices=['escherichia_coli', 'salmonella', 'plasmodium_falciparum', 'neisseria_gonorrhoeae', 'mycobacterium_tuberculosis', 'enterococcus_faecalis', 'enterococcus_faecium'], help="Reference specie to use for PointFinder")
parser.add_argument("--specie_mlst", metavar='SPECIE', choices=['abaumannii', 'abaumannii_2', 'achromobacter', 'aeromonas', 'afumigatus', 'aphagocytophilum', 'arcobacter', 'bbacilliformis', 'bcc', 'bcereus', 'bhampsonii', 'bhenselae', 'bhyodysenteriae', 'bintermedia', 'blicheniformis', 'bordetella', 'borrelia', 'bpilosicoli', 'bpseudomallei', 'brachyspira', 'brucella', 'bsubtilis', 'calbicans', 'campylobacter', 'cbotulinum', 'cconcisus', 'cdifficile', 'cdiphtheriae', 'cfetus', 'cfreundii', 'cglabrata', 'chelveticus', 'chlamydiales', 'chyointestinalis', 'cinsulaenigrae', 'ckrusei', 'clanienae', 'clari', 'cmaltaromaticum', 'cronobacter', 'csepticum', 'csinensis', 'csputorum', 'ctropicalis', 'cupsaliensis', 'dnodosus', 'ecloacae', 'ecoli', 'ecoli_2', 'edwardsiella', 'efaecalis', 'efaecium', 'fpsychrophilum', 'ganatis', 'hcinaedi', 'hinfluenzae', 'hparasuis', 'hpylori', 'hsuis', 'kaerogenes', 'kkingae', 'koxytoca', 'kpneumoniae', 'kseptempunctata', 'leptospira', 'leptospira_2', 'leptospira_3', 'liberibacter', 'llactis', 'lmonocytogenes', 'lsalivarius', 'mabscessus', 'magalactiae', 'mbovis', 'mcanis', 'mcaseolyticus', 'mcatarrhalis', 'mhaemolytica', 'mhyopneumoniae', 'mhyorhinis', 'miowae', 'mmassiliense', 'mplutonius', 'mpneumoniae', 'msynoviae', 'mycobacteria', 'neisseria', 'orhinotracheale', 'otsutsugamushi', 'pacnes', 'paeruginosa', 'pdamselae', 'pfluorescens', 'pgingivalis', 'plarvae', 'pmultocida_multihost', 'pmultocida_rirdc', 'ppentosaceus', 'pputida', 'psalmonis', 'ranatipestifer', 'rhodococcus', 'sagalactiae', 'saureus', 'sbsec', 'scanis', 'sdysgalactiae', 'senterica', 'sepidermidis', 'sgallolyticus', 'shaemolyticus', 'shominis', 'sinorhizobium', 'slugdunensis', 'smaltophilia', 'soralis', 'sparasitica', 'spneumoniae', 'spseudintermedius', 'spyogenes', 'ssuis', 'sthermophilus', 'sthermophilus_2', 'streptomyces', 'suberis', 'szooepidemicus', 'taylorella', 'tenacibaculum', 'tpallidum', 'tvaginalis', 'ureaplasma', 'vcholerae', 'vcholerae2', 'vibrio', 'vparahaemolyticus', 'vtapetis', 'vvulnificus', 'wolbachia', 'xfastidiosa', 'yersinia', 'ypseudotuberculosis', 'yruckeri'], help="Reference specie to use for MLST")
parser.add_argument("--scheme", metavar='SCHEME', choices=['incac', 'incf', 'inchi1', 'inchi2', 'inci1', 'incn'], help="Reference scheme to use for pMLST")
parser.add_argument("--allmut", action="store_true", help="Look for uncharaterized mutation with PointFinder")
parser.add_argument("--extendedout", action="store_true", help="Give extended output of PlasmidFinder, MLST and pMLST")
parser.add_argument("-v", "--verbiosity", action="store_true", help="Retrieve stdout of commands")
args = parser.parse_args()

def run_resfinder(assembly, out_dir, identity, coverage, verbiosity):
    res_db = "Software/resfinder/resfinder_db"
    prog = "resfinder"
    outdir = create_outdir(out_dir, prog)
    resfinder_command = ["resfinder", "-i", assembly, "-o", outdir,
                         "-b", "blastn", "-p", res_db,
                         "-l", str(coverage), "-t", str(identity)]

    resfinder = sb.Popen(resfinder_command, stdout=sb.PIPE, stderr=sb.PIPE)
    resfinder_out, resfinder_err = resfinder.communicate()

    if verbiosity:
        print(resfinder_out, resfinder_err, flush=True)

def run_pointfinder(assembly, out_dir, identity, coverage, specie, allmut, verbiosity):
    point_db = "Software/pointfinder/pointfinder_db"
    prog = "pointfinder"
    blastn_full_path=".conda/envs/biopy/bin/blastn"
    outdir = create_outdir(out_dir, prog)
    pointfinder_command = ["pointfinder", "-i", assembly, "-o", outdir,
                           "-s", specie, "-p", point_db,
                           "-m", "blastn", "-m_p", blastn_full_path, 
                           "-l", str(coverage), "-t", str(identity)]

    if allmut:
        pointfinder_command = pointfinder_command + ["-u"]

    pointfinder = sb.Popen(pointfinder_command, stdout=sb.PIPE, stderr=sb.PIPE)
    pointfinder_out, pointfinder_err = pointfinder.communicate()

    if verbiosity:
        print(pointfinder_out, pointfinder_err, flush=True)

def run_plasmidfinder(assembly, out_dir, identity, coverage, extendedout, verbiosity):
    plasmid_db = "Software/plasmidfinder/plasmidfinder_db"
    prog = "plasmidfinder"
    outdir = create_outdir(out_dir, prog)
    plasmidfinder_command = ["plasmidfinder", "-i", assembly, "-o", outdir,
                             "-p", plasmid_db, "-mp", "blastn", 
                             "-l", str(coverage), "-t", str(identity)]

    if extendedout:
        plasmidfinder_command = plasmidfinder_command + ["-x"]

    plasmidfinder = sb.Popen(plasmidfinder_command, stdout=sb.PIPE, stderr=sb.PIPE)
    plasmidfinder_out, plasmidfinder_err = plasmidfinder.communicate()

    if verbiosity:
        print(plasmidfinder_out, plasmidfinder_err, flush=True)


def run_mlst(assembly, out_dir, specie, extendedout, verbiosity):
    mlst_db = "Software/mlst/mlst_db"
    prog = "mlst"
    outdir = create_outdir(out_dir, prog)
    mlst_command = ["mlst", "-i", assembly, "-o", outdir,
                    "-s", specie, "-p", mlst_db,
                    "-mp", "blastn"]
    if extendedout:
        mlst_command = mlst_command + ["-x"]

    mlst = sb.Popen(mlst_command, stdout=sb.PIPE, stderr=sb.PIPE)
    mlst_out, mlst_err = mlst.communicate()

    if verbiosity:
        print(mlst_out, mlst_err, flush=True)


def run_pmlst(assembly, out_dir, scheme, extendedout, verbiosity):
    pmlst_db = "Software/pmlst/pmlst_db"
    prog = "pmlst"
    outdir = create_outdir(out_dir, prog)
    pmlst_command = ["pmlst", "-i", assembly, "-o", outdir,
                    "-s", scheme, "-p", pmlst_db,
                    "-mp", "blastn"]
    if extendedout:
        pmlst_command = pmlst_command + ["-x"]

    pmlst = sb.Popen(pmlst_command, stdout=sb.PIPE, stderr=sb.PIPE)
    pmlst_out, pmlst_err = pmlst.communicate()

    if verbiosity:
        print(pmlst_out, pmlst_err, flush=True)


def create_outdir(out_dir, prog):
    outdir = os.path.join(out_dir, prog)
    if not os.path.exists(outdir):
        os.mkdir(outdir, 0o755)
    return outdir

def main(args):
    if args.all:
        run_resfinder(args.assembly, args.folder, args.residentity, args.rescoverage, args.verbiosity)
        run_pointfinder(args.assembly, args.folder, args.pointidentity, args.pointcoverage, args.specie, args.allmut, args.verbiosity)
        run_plasmidfinder(args.assembly, args.folder, args.plasmidentity, args.plasmcoverage, args.extendedout, args.verbiosity)
        run_mlst(args.assembly, args.folder, args.specie_mlst, args.extendedout, args.verbiosity)
    else:
        if args.res:
            run_resfinder(args.assembly, args.folder, args.residentity, args.rescoverage, args.verbiosity)
        if args.point:
            run_pointfinder(args.assembly, args.folder, args.pointidentity, args.pointcoverage, args.specie, args.allmut, args.verbiosity)
        if args.plasmid:
            run_plasmidfinder(args.assembly, args.folder, args.plasmidentity, args.plasmcoverage, args.extendedout, args.verbiosity)
        if args.mlst:
            run_mlst(args.assembly, args.folder, args.specie_mlst, args.extendedout, args.verbiosity)
        if args.pmlst:
            run_pmlst(args.assembly, args.folder, args.scheme, args.extendedout, args.verbiosity)


main(args)
