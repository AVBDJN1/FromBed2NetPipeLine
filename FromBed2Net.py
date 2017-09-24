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


mode = sys.argv[1]

if mode == "--help":
    try:
        print open("../README", "r").read()
    except:
        print open("README", "r").read()
        
if mode == "-HP_selection":
    input_file = sys.argv[2]
    col_of_threshold, threshold, ncol_to_divide = sys.argv[3].split("-")
    criteria_divider(score_bed_selector(input_file, col_of_threshold, threshold), ncol_to_divide)    
    
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
    first_names = os.path.splitext(os.path.basename(input_file))[0].split("_")[1]
    col_of_threshold, threshold, ncol_to_divide = sys.argv[3].split("-")
    input_dir_file = criteria_divider(score_bed_selector(input_file, col_of_threshold, threshold), ncol_to_divide)
    annotation_source = sys.argv[4]
    output_folder = sys.argv[5]
    if annotation_source.split(":")[0] == "-download":
        gff3_download(annotation_source.split(":")[1])
        annotation(input_dir_file, annotation_source.split(":")[1]+".gff3", output_folder)
        tsver(output_folder)
        gene_extractor(output_folder)        
    else:
        annotation(input_dir_file, annotation_source, output_folder)
        tsver(output_folder)
        gene_extractor(output_folder)        
    os.system("mv {}* {}".format(first_names, output_folder))
