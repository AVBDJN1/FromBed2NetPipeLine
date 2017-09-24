#!bin/usr/python
# -*- coding: utf-8 -*-.
# Adrian Garcia Moreno

import io
import re as re
import os 
import subprocess
import sys

def gene_extractor(input_dir_file):
    i = 1
    genes_counter = 0
    paths, input_files = input_checker(input_dir_file)
# This forcing can be problematic if there is no a folder with that name... 
# but I will trust that people do not change my nomenclature
    for path in paths:
        if re.match('.*tsvs/', path): 
            input_files = os.listdir(path)
            break
    length = len(input_files)
    geneID_match = re.compile(r"\tgene\t[0-9]+.*ID\t([0-9]+|ENSG[0-9]+)\.*[0-9]*\t")
    os.system("mkdir "+path+"/../genes_lists")
    for tsv in input_files:
        print "Extracting genes of {}\t{} of {}".format(tsv, i, length)
        i += 1
        geneID_list_name = os.path.splitext(os.path.basename(tsv))[0].split("_")[-1]
        geneID_list = open(path+"../genes_lists/GeneIDs_"+geneID_list_name+".txt", "w")
        with io.open(path+tsv,'r',encoding="utf-8") as tsv_handle:
            for line in tsv_handle:
                if geneID_match.search(line) is not None:
                    geneID_list.write(geneID_match.search(line).group(1)+" "+line.split("\t")[4]+"\n")
                    genes_counter += 1
        geneID_list.close()
    os.system("ls -1 "+path+"../genes_lists/ | wc -l")
    os.system("find "+path+"../genes_lists/. -name 'GeneIDs_*' | xargs wc -l")
    print "Your lists of genes are ready, a total of {} genes were obtained".format(genes_counter)
    return path+"/../genes_lists"
    
    
