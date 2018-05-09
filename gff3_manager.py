#!bin/usr/python
# -*- coding: utf-8 -*-.
# Adrian Garcia Moreno

import os 
import subprocess
import re as re
from input_controller import *
import argparse

# To Develop this tool
# MUST USE ! ftp://ftp.ncbi.nlm.nih.gov/genomes/README_GFF3.txt 

def gff3_geneExtractor(input_gff3, not_chr):
    print("Processing {}".format(input_gff3))
    output_gff3 = "GeneIDs_{}{}".format(*getFilename(input_gff3))
    output_gff3 = gzip.open(output_gff3, "wt")

    # Need to solve the problem with the NT_ and NW_ 
    # Since I cannot get the chr from those I will only work with the NC_    
    gene_line = re.compile(r"^NC.*\tgene\t.*GeneID:([0-9]+)")
    
    # gene_line should be the regex for any given nomenclature
    # desired to be extracted and it should be an argument to the
    # function. It will be implemented with the nomenclature identifier     
    
    input_gff3 = openFile(input_gff3)
    
    if not_chr:
        for line in input_gff3:
            if gene_line.search(line) is not None:
                chro = chr_from_accession(line.split("\t")[0])
                start = line.split("\t")[3]
                end = line.split("\t")[4]
                geneID = gene_line.search(line).group(1)
                output_gff3.write("\t".join([chro,start,end,geneID])+"\n")
        print("Genes extracted and Transformed\nYour file: {}".format(output_gff3.name))
    else:
        for line in input_gff3:
            if gene_line.search(line) is not None:
                accesion = line.split("\t")[0]
                start = line.split("\t")[3]
                end = line.split("\t")[4]
                geneID = gene_line.search(line).group(1)
                output_gff3.write("\t".join([accesion,start,end,geneID])+"\n")
        print("Genes extracted\nYour file: {}".format(output_gff3.name))
    input_gff3.close()
    output_gff3.close()        

def chr_from_accession(accession):
    # Very important
    # https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly    
    number_of_accession = int(accession.split(".")[0].split("_")[1])    
    if number_of_accession == 23:
        chro = "chrX"
    elif number_of_accession == 24:
        chro = "chrY"
    elif len(str(number_of_accession)) > 2:
        chro = "chrM"
    else: 
        chro = "chr"+str(number_of_accession)
    return chro

def gff3_download(gff3_source):
    # https://www.gencodegenes.org/releases/(19,20).html
    # ftp://ftp.ncbi.nih.gov/genomes/Homo_sapiens/
    if gff3_source == "gcode_hg19":
        gff3_url = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gff3.gz"
    elif gff3_source == "gcode_hg20":
        gff3_url = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_20/gencode.v20.annotation.gff3.gz"
    elif gff3_source == "ncbi_hg19":
        gff3_url = "ftp://ftp.ncbi.nih.gov/genomes/Homo_sapiens/GRCh37.p13_interim_annotation/interim_GRCh37.p13_top_level_2017-01-13.gff3.gz"
    else:# gff3_source == "ncbi_hg20":
        gff3_url = "ftp://ftp.ncbi.nih.gov/genomes/Homo_sapiens/GFF/ref_GRCh38.p12_top_level.gff3.gz"

    print("Downloading {} GFF3 from\n{}".format(gff3_source, gff3_url))
    subprocess.run(["wget", gff3_url], 
    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    gff3_donwloaded = "{}{}".format(*getFilename(gff3_url))
    print("Downloaded {}:\n\t{}".format(gff3_source, gff3_donwloaded))
    return(gff3_donwloaded)

parser = argparse.ArgumentParser(  
    prog = 'GFF3 Manager',
    prefix_chars = '-',
    formatter_class = argparse.ArgumentDefaultsHelpFormatter,
    description = '''To download gff3 files and extract genes (GeneID by 
    now) and if the gff3 comes from ncbi to transform from accesion 
    numbers (only NC_xxxxx.x) to chromosomes location nomenclature
    unless the opposite is indicated with its flag.''')

input_group = parser.add_mutually_exclusive_group(required=True)

input_group.add_argument('-gff3','--input_gff3', 
help='GFF3 file as input', nargs='?')

input_group.add_argument('-down_gff3','--gff3_source', 
help='GFF3 source to download', nargs='?', 
choices=['gcode_hg19','gcode_hg20','ncbi_hg19','ncbi_hg20'])

parser.add_argument('-not_chr','--not_chr', 
help='Flag to not transform RefSeq NC_xxx.x to chromosome', action="store_false")

args = parser.parse_args()

print(args.not_chr)
if args.gff3_source is not None:
    input_gff3 = gff3_download(args.gff3_source)
    gff3_geneExtractor(input_gff3, args.not_chr)

else:
    gff3_geneExtractor(args.input_gff3, args.not_chr)

# Although it is not the case, I must improve this to not extract 
# something from the GFF3
