#!bin/usr/python
# -*- coding: utf-8 -*-.
# Adrian Garcia Moreno

import os 
import subprocess
import re as re
from input_controller import *
import argparse
import sys
import time
import urllib

# To Develop this tool
# MUST USE ! ftp://ftp.ncbi.nlm.nih.gov/genomes/README_GFF3.txt 

def consortium_regex(consortium):
    if consortium == "ncbi":
        regex = re.compile(r"^NC.*\tgene\t.*GeneID:([0-9]+)")
    # Need to solve the problem with the NT_ and NW_ 
    # Since I cannot get the chr from those I will only work with the NC_    

    elif consortium == "gcode":
        regex = re.compile(r"\tgene\t.*gene_id=(ENSG[0-9]+)")
    
    # more regex for any given nomenclature should be added according to
    # to be able and extract genes and it should be an argument to the
    # function. It will be implemented with the nomenclature identifier     
    return regex

def gff3_geneExtractor(input_gff3, to_chr, consortium):
    print("Processing {}".format(input_gff3))
    output_gff3 = "GeneIDs_{}{}".format(*getFilename(input_gff3))
    output_gff3 = gzip.open(output_gff3, "wt")
    gene_line = consortium_regex(consortium)    
    input_gff3 = openFile(input_gff3)
    
    if to_chr:
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
        gff3_url = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gff3.gz"
    elif gff3_source == "ncbi_hg19":
        gff3_url = "ftp://ftp.ncbi.nih.gov/genomes/Homo_sapiens/GRCh37.p13_interim_annotation/interim_GRCh37.p13_top_level_2017-01-13.gff3.gz"
    elif gff3_source == "ncbi_hg20":
        gff3_url = "ftp://ftp.ncbi.nih.gov/genomes/Homo_sapiens/GFF/ref_GRCh38.p12_top_level.gff3.gz"
    else: 
        # Makes little sense this... because I force the options with argparse
        print("Something went wrong")

    print("Downloading {} GFF3 from\n{}".format(gff3_source, gff3_url))
    gff3_donwloaded = "{}{}".format(*getFilename(gff3_url))
    save(gff3_url, gff3_donwloaded)
    
    print("Downloaded {}:\n\t{}".format(gff3_source, gff3_donwloaded))
    return(gff3_donwloaded)

parser = argparse.ArgumentParser(  
    prog = 'GFF3 Manager',
    prefix_chars = '-',
    formatter_class = argparse.ArgumentDefaultsHelpFormatter,
    description = '''To download gff3 files and extract genes (GeneID by 
    now) and if the gff3 comes from ncbi you can transform the accesion 
    numbers (only NC_xxxxx.x) to chromosomes location nomenclature
    specifying with its flag.''')

input_group = parser.add_mutually_exclusive_group(required=True)

input_group.add_argument('-gff3','--input_gff3', 
help='GFF3 file as input', nargs='?')

parser.add_argument('-consortium','--consortium', 
help='If an gff3 file is directly given as input, you must specify the\
consortium which it comes from', choices=["ncbi","gcode"])

input_group.add_argument('-down_gff3','--gff3_source', 
help='GFF3 source url (Warning both from ncbi might need to be \
transformed from accession to chr if that is tour case.\
Genecode here given is already in chr location base \
(Shall I  do an opposite transformation from chr to accesion?)', nargs='?', 
choices=['gcode_hg19','ncbi_hg19', 'gcode_hg20','ncbi_hg20'])

parser.add_argument('-to_chr','--to_chr', 
help='Flag to transform RefSeq NC_xxx.x to chromosome', 
action="store_true")

args = parser.parse_args()

if args.gff3_source is not None:
    consortium = args.gff3_source.split("_")[0]
    input_gff3 = gff3_download(args.gff3_source)
    gff3_geneExtractor(input_gff3, args.to_chr, consortium)

else:
    gff3_geneExtractor(args.input_gff3, args.to_chr, args.consortium)

# Although it is not the case, I must improve this to not extract 
# something from the GFF3
