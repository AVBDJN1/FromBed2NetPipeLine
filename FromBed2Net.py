#!bin/usr/python
# -*- coding: utf-8 -*-.
# Adrian Garcia Moreno

import io
import re as re
import os 
import subprocess
import sys

from input_checker import input_checker
from criteria_divider import score_bed_selector, criteria_divider
from gff3s_manager  import ncbi_formatter, gff3_download
from bed_tsv import annotation, tsver
from gene_extractor import gene_extractor

script_path = os.path.dirname(sys.argv[0])
mode = sys.argv[1]

if mode == "--help":
    try:
        print open("../README", "r").read()
    except:
        print open("README", "r").read()

if mode == "-download_annot":
    annotation_source = sys.argv[2]
    gff3_download(annotation_source.split(":")[1])

        
if mode == "-HP_selection":
    input_file = sys.argv[2]
    col_of_score, threshold, col_of_feature = sys.argv[3].split(",")
    criteria_divider(score_bed_selector(input_file, col_of_score, threshold), col_of_feature)    
    
if mode == "-annotation":
    input_dir_file = sys.argv[2]
    annotation_source = sys.argv[3]
    output_folder = sys.argv[4]
    if annotation_source.split(":")[0] == "-download":
        gff3_download(annotation_source.split(":")[1])
        annotation(input_dir_file, annotation_source.split(":")[1]+".gff3", output_folder)
    else:
        annotation(input_dir_file, annotation_source, output_folder)

if mode == "-tsver":
    input_dir_file = sys.argv[2]
    tsver(input_dir_file)
    
if mode == "-extract_genes":
    input_dir_file = sys.argv[2]
    gene_extractor(input_dir_file)

if mode == "-full":
    input_file = sys.argv[2]
    first_names = os.path.splitext(os.path.basename(input_file))[0]
    col_of_score, threshold, col_of_feature = sys.argv[3].split(",")
    output_folder = sys.argv[5]
    
    if output_folder[-1] != "/":
        output_folder = output_folder+"/"

    try:
        os.mkdir(output_folder)
    except:
        answer = raw_input("Folder already exists, some files could be overwritten\nDo you want to continue? (Y/N)")
        if answer == "N":
            sys.exit()

    input_dir_file = criteria_divider(score_bed_selector(input_file, col_of_score, threshold, output_folder), col_of_feature, output_folder)


    annotation_source = sys.argv[4]
    
    if annotation_source.split(":")[0] == "-download":
        gff3_download(annotation_source.split(":")[1])
        annotation(input_dir_file, "{}.gff3".format(annotation_source.split(":")[1]), output_folder)
        tsver(output_folder)
        gene_extractor(output_folder, col_of_score)        
    else:
        annotation(input_dir_file, annotation_source, output_folder)
        tsver(output_folder)
        gene_extractor(output_folder, col_of_score)        
    #os.system("mv {}* {}".format(first_names, output_folder))
    print "WARNING if your initial file is also in this folder you will find it in the output directory"

    os.system("Rscript {}/TopGOer.r {}/genes_lists/ 9606_geneID2GO.map {}".format(script_path, output_folder, first_names))
    os.system("python {}/sumarizer.py {} 0.003".format(script_path, output_folder))
    
    genes_of_goes_paths = subprocess.check_output("find {}/*GOanalysis/*/genes_of_goes -type d".format(output_folder), shell=True).rstrip("\n")
    for genes_of_goes_path in genes_of_goes_paths.split("\n"):
        os.system("python {}/HTMLizer.py {}/".format(script_path, genes_of_goes_path))
        
    os.system("python {}/HTMLizer.py {}/*summed_up_annot/".format(script_path, output_folder))
