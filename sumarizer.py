#!bin/usr/python
# -*- coding: utf-8 -*-.
# Adrian Garcia Moreno

import io
import re
import os 
import subprocess
import sys

def genes_loci_score_writer(path_to_tsv, summary_handle):    
    
    tsv_handle = open(path_to_tsv, "r")
    geneID_match = re.compile(r"\tgene\t[0-9]+.*ID\t([0-9]+|ENSG[0-9]+)\.*[0-9]*\t")
    
    loci_genes = {}
    
    for line in tsv_handle:
        if geneID_match.search(line) is not None:
            locus_score = "\t".join(line.split("\t")[0:3]+[line.split("\t")[4]])
            gene = geneID_match.search(line).group(1)
            if locus_score in loci_genes:
                loci_genes[locus_score].append(gene)   
            else:
                loci_genes[locus_score] = [gene]
            
    tsv_handle.close()
    
    summary_handle.write("{}\n".format(feature_name))
    if loci_genes == {}:
        summary_handle.write("Not a single gene was found so that means that,\n")
    else:
        for locus_score in loci_genes:
            summary_handle.write("{}\n{}\n\n".format(locus_score, "\n".join(loci_genes[locus_score])))

def GO_analysis_writter(go_results_path, threshold, summary_handle):
    errors = []
    feature_name = go_results_path.split("/")[-2]    
    go_result_matches = (go_result for go_result in os.listdir(go_results_path) if re.search(".txt$", go_result) is not None)
    MF, CC, BP = False, False, False # if one of these comes true it means that none have passed the threshold of the summaryzer

    for go_result in go_result_matches:
                
        if go_result[0:2] == "MF":
            go_result_path = go_results_path+go_result
            go_result_read = open(go_result_path, "r")
            summary_handle.write("###MOLECULAR FUNCTION###\n")
            for line in go_result_read:
                try:
                    score = float(line.split("\t")[-1])
                    if score <= threshold:
                        MF = True
                        summary_handle.write(line)
                except:
                    summary_handle.write(line)
            summary_handle.write("\n") # this is \n to mantain the format for the HTMLizer
            
        if go_result[0:2] == "CC":
            go_result_path = go_results_path+go_result
            go_result_read = open(go_result_path, "r")
            summary_handle.write("###CELULAR COMPONENT###\n")
            for line in go_result_read:
                try:
                    score = float(line.split("\t")[-1])
                    if score <= threshold:
                        CC = True
                        summary_handle.write(line)
                except:
                    summary_handle.write(line)
            summary_handle.write("\n")
            
        if go_result[0:2] == "BP":
            go_result_path = go_results_path+go_result
            go_result_read = open(go_result_path, "r")
            summary_handle.write("###BIOLOGICAL PROCESS###\n")
            for line in go_result_read:
                try:
                    score = float(line.split("\t")[-1])
                    if score <= threshold:
                        BP = True
                        summary_handle.write(line)
                except:
                    summary_handle.write(line)
            summary_handle.write("\n")
            
        if CC and BP and MF:
            errors.append(feature_name)
    
    if len(errors) >= 1:
        errors_summary_handle = open("{}_errors_summary.txt".format(threshold), "w")
        errors_summary_handle.write("The following features did not pass the summarizer threshold:\n{}".format("\n".join(errors)))
        errors_summary_handle.close()
            
input_folder = sys.argv[1]
if len(sys.argv) >= 3:
    threshold = float(sys.argv[2])
else:
    threshold = 0.05

annot_folders = os.listdir(input_folder)

tsv_list = ["{}tsvs/{}".format(input_folder, tsv_file) for tsv_file in os.listdir("{}tsvs/".format(input_folder))] 
                                                       # this loop give us directly the full path

go_results_folder = next(x for x in annot_folders if re.search(".*GOanalysis", x) is not None) # to find the name of the folder
go_result_folders = os.listdir("{}{}/".format(input_folder, go_results_folder))

summary_output = "{}/{}_summed_up_annot/".format(input_folder, threshold)

try:
    os.mkdir(summary_output)
except:
    answer = raw_input("Folder already exists, some files could be overwritten\nDo you want to continue? (Y/N)")
    if answer == "N":
        sys.exit()
        
for index in xrange(len(tsv_list)):
    feature_name = "_".join(os.path.splitext(os.path.basename(tsv_list[index]))[0].split("_")[1:])
    
    if feature_name in go_result_folders:
        go_results_path = "{}{}/{}/".format(input_folder, go_results_folder, feature_name)
        summary_handle = open("{}{}_summary.txt".format(summary_output, feature_name), "w")
        genes_loci_score_writer(tsv_list[index], summary_handle)
        GO_analysis_writter(go_results_path, threshold, summary_handle)
        summary_handle.close()

    else:
        summary_handle = open("{}{}_summary.txt".format(summary_output, feature_name), "w")
        genes_loci_score_writer(tsv_list[index], summary_handle)
        summary_handle.write("{} has no genes mapped to GO terms\n".format(feature_name))
        summary_handle.close()

if os.path.isfile("{}_errors_summary.txt".format(threshold)):
    os.system("mv *errors_summary.txt {}".format(summary_output))
