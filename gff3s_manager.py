#!bin/usr/python
# -*- coding: utf-8 -*-.
# Adrian Garcia Moreno

import io
import re as re
import os 
import subprocess
import sys

def ncbi_formatter(source):
    output = subprocess.check_output("grep ^NC "+source+"_raw.gff3 | cut -f1 | sort -u", shell=True)
    contigs = output.split("\n")
    if contigs[-1] == "":
        contigs = contigs[:-1]

    ncbicontig_chr_dict = {}
    for contig in contigs:
        contig_chr = int(contig.split(".")[0].split("_")[1])
        if contig_chr == 23:
            chro = "chrX"
        elif contig_chr == 24:
            chro = "chrY"
        elif len(str(contig_chr)) > 2:
            chro = "chrM"
        else: 
            chro = "chr"+str(contig_chr)
        ncbicontig_chr_dict[contig] = chro
      
    original_ncbigff3 = open(source+"_raw.gff3", "r")
    new_ncbigff3 = open(source+".gff3", "w")

    for line in original_ncbigff3:
        tabed_elements = line.split("\t")
        if tabed_elements[0] in ncbicontig_chr_dict:
            bed_str = ncbicontig_chr_dict[tabed_elements[0]]+"\t"+"\t".join(tabed_elements[1:10])
            new_ncbigff3.write(bed_str)
        else:
            continue
          
    original_ncbigff3.close()
    new_ncbigff3.close()

def gff3_download(source):
    # https://www.gencodegenes.org/releases/(19,20).html
    # ftp://ftp.ncbi.nih.gov/genomes/Homo_sapiens/
    if source == "gcode_hg19":
        os.system("wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gff3.gz")
        print "Downloaded, time to gunzip"
        os.system("gunzip -c gencode.v19.annotation.gff3.gz > gcode_hg19.gff3")
    if source == "gcode_hg20":
        os.system("wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_20/gencode.v20.annotation.gff3.gz")
        print "Downloaded, time to gunzip"
        os.system("gunzip -c gencode.v20.annotation.gff3.gz > gcode_hg20.gff3")
    if source == "ncbi_hg19":
        os.system("wget ftp://ftp.ncbi.nih.gov/genomes/Homo_sapiens/GRCh37.p13_interim_annotation/interim_GRCh37.p13_top_level_2017-01-13.gff3.gz")    
        print "Downloaded, time to gunzip"
        os.system("gunzip -c interim_GRCh37.p13_top_level_2017-01-13.gff3.gz > ncbi_hg19_raw.gff3")
        ncbi_formatter(source)
    if source == "ncbi_hg20":
        os.system("wget ftp://ftp.ncbi.nih.gov/genomes/Homo_sapiens/GFF/ref_GRCh38.p7_top_level.gff3.gz")
        print "Downloaded, time to gunzip"
        os.system("gunzip -c ref_GRCh38.p7_top_level.gff3.gz > ncbi_hg20_raw.gff3")
        ncbi_formatter(source)
