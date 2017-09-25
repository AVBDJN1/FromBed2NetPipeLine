#!bin/usr/python
# -*- coding: utf-8 -*-.
# Adrian Garcia Moreno

import io
import re as re
import os 
import subprocess
import sys

from input_checker import input_checker

def gene_extractor(input_dir_file):
    counter = 1
    genes_counter = 0
    paths, input_files = input_checker(input_dir_file)
    
    if os.path.isfile(input_files[0]):
       print "Single gene extraction mode with {}".format(os.path.basename(input_files[0]))
       path = paths[0] 
       tsv_handle = open(input_files[0],'r')
    else:
        for path in paths:
            if re.match('.*tsvs/', path): 
                input_files = os.listdir(path)
                break
                
    to_process = len(input_files)
    geneID_match = re.compile(r"\tgene\t[0-9]+.*ID\t([0-9]+|ENSG[0-9]+)\.*[0-9]*\t")
    os.system("mkdir "+path+"/../genes_lists")
    
    if len(input_files) <= 1:
        geneID_list_name = os.path.splitext(os.path.basename(input_files[0]))[0].split("_")[-1]
        print "Extracting genes of {}\t{} of {}".format(input_files[0], counter, to_process)
        geneID_list = open(path+"../genes_lists/GeneIDs_"+geneID_list_name+".txt", "w")
        for line in tsv_handle:
            if geneID_match.search(line) is not None:
                geneID_list.write(geneID_match.search(line).group(1)+" "+line.split("\t")[4]+"\n")
                genes_counter += 1
        geneID_list.close()
        tsv_handle.close()
    else:
        for tsv in input_files:
            geneID_list_name = os.path.splitext(os.path.basename(tsv))[0].split("_")[-1]
            print "Extracting genes of {}\t{} of {}".format(tsv, counter, to_process)
            counter += 1
            geneID_list = open(path+"../genes_lists/GeneIDs_"+geneID_list_name+".txt", "w")
            tsv_handle = open(path+tsv,'r')
            for line in tsv_handle:
                if geneID_match.search(line) is not None:
                    geneID_list.write(geneID_match.search(line).group(1)+" "+line.split("\t")[4]+"\n")
                    genes_counter += 1
            geneID_list.close()
            tsv_handle.close()

    os.system("ls -1 "+path+"../genes_lists/ | wc -l")
    os.system("find "+path+"../genes_lists/. -name 'GeneIDs_*' | xargs wc -l")
    print "Your lists of genes are ready, a total of {} genes were obtained".format(genes_counter)
    return path+"/../genes_lists"
    
    
