#!bin/usr/python
# -*- coding: utf-8 -*-.
# Adrian Garcia Moreno

import io
import re as re
import os 
import subprocess
import sys
from input_checker import input_checker

def gene_extractor(input_dir_file, col_of_score = 0):
    counter = 1
    genes_counter = 0
    paths, input_files = input_checker(input_dir_file)
    
    if os.path.isfile(input_files[0]) and len(input_files) == 1:
        path = paths[0] 
    else:
        for path in paths:
            if re.match('.*tsvs/*', path) is not None: 
                input_files = [path+tsv for tsv in os.listdir(path)]
                break
                
    to_process = len(input_files)
    
    try:
        os.mkdir("{}/../genes_lists".format(path))
    except:
        answer = raw_input("Folder already exists, some files could be overwritten\nDo you want to continue? (Y/N)")
        if answer == "N":
            sys.exit()

    geneID_match = re.compile(r"\tgene\t[0-9]+.*GeneID\t([0-9]+|ENSG[0-9]+)\.*[0-9]*\t")
        
    for tsv in input_files:
        print "Extracting genes of {}\t{} of {}".format(tsv, counter, to_process)
        counter += 1
        
        geneID_list_name = "_".join(os.path.splitext(os.path.basename(tsv))[0].split("_")[1:])
        geneID_list = open("{}../genes_lists/GeneIDs_{}.txt".format(path, geneID_list_name), "w")

        tsv_handle = open(tsv,'r')
        for line in tsv_handle:
            if geneID_match.search(line) is not None and col_of_score != 0:                
                geneID_list.write("{} {}\n".format(geneID_match.search(line).group(1), line.split("\t")[int(col_of_score)-1]))
                genes_counter += 1
            elif geneID_match.search(line) is not None and col_of_score == 0:
                geneID_list.write("{}\n".format(geneID_match.search(line).group(1)))
                genes_counter += 1
        geneID_list.close()
        tsv_handle.close()
    
    os.system("ls -1 {}../genes_lists/ | wc -l".format(path))
    os.system("find {}../genes_lists/. -name 'GeneIDs_*' | xargs wc -l".format(path))
    print "Your lists of genes are ready, a total of {} genes were obtained".format(genes_counter)
    return "{}../genes_lists".format(path)
# return "{}../genes_lists".format(path)
    
    
