#!bin/usr/python
# -*- coding: utf-8 -*-.
# Adrian Garcia Moreno

import io
import re
import os 
import subprocess
import sys

def genes_loci_score_writer(path_to_tsv, summary_handle):    
    
    tsv_handle = open(path_to_tsv, "r")#input_folder+"tsvs/"+tsv_list[index], "r")
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
    
    summary_handle.write(">>{}\n".format(HPname))
    if loci_genes == {}:
        summary_handle.write("Not a single gene was found so that means that,\n".format(HPname))
    else:
        for locus_score in loci_genes:
            summary_handle.write("{}\n{}\n\n".format(locus_score, "\n".join(loci_genes[locus_score])))

def GO_analysis_writter(go_results_path, summary_handle):
    
    go_result_matches = (go_result for go_result in os.listdir(go_results_path) if re.search("[0-9]$", go_result) is not None)
     
    for go_result in go_result_matches:
                
        if go_result[0:2] == "MF":
            go_result_path = go_results_path+go_result#input_folder+go_results_folders+"/results/"+go_results[0]+"/"+i
            go_result_read = open(go_result_path, "r").read()
            summary_handle.write("###MOLECULAR FUNCTION {}###\n".format(HPname))
            summary_handle.write("{}\n".format(go_result_read))
            
        if go_result[0:2] == "CC":
            go_result_path = go_results_path+go_result#input_folder+go_results_folders+"/results/"+go_results[0]+"/"+i
            go_result_read = open(go_result_path, "r").read()
            summary_handle.write("###CELULAR COMPONENT {}###\n".format(HPname))
            summary_handle.write("{}\n".format(go_result_read))
            
        if go_result[0:2] == "BP":
            go_result_path = go_results_path+go_result#input_folder+go_results_folders+"/results/"+go_results[0]+"/"+i
            go_result_read = open(go_result_path, "r").read()
            summary_handle.write("###BIOLOGICAL PROCESS {}###\n".format(HPname))
            summary_handle.write("{}\n".format(go_result_read))
            
input_folder = sys.argv[1]
#input_folder = "/home/bioinfo/Desktop/RESULTS/hg19_ncbi_run/annotation_results/"

annot_folders = os.listdir(input_folder)

tsv_list = os.listdir(input_folder+"tsvs/")

go_results_folder = next(x for x in annot_folders if re.search(".*GOanalysis", x) is not None)
go_result_folders = os.listdir(input_folder+go_results_folder+"/results")


for index in xrange(len(tsv_list)):
    HPname = os.path.splitext(tsv_list[index])[0][-10:]
    
    path_to_tsv = "{}tsvs/{}".format(input_folder, tsv_list[index])
    if HPname in go_result_folders:
        go_results_path = "{}{}/results/{}/".format(input_folder, go_results_folder, HPname)
        summary_handle = open("{}_summary.txt".format(HPname), "w")
        genes_loci_score_writer(path_to_tsv, summary_handle)
        GO_analysis_writter(go_results_path, summary_handle)
        summary_handle.close()

    else:
        summary_handle = open("{}_summary.txt".format(HPname), "w")
        genes_loci_score_writer(path_to_tsv, summary_handle)
        summary_handle.write("{} has no genes mapped to GO terms\n".format(os.path.splitext(tsv_list[index])[0]))
        summary_handle.close()

os.mkdir("{}../summed_up_annot/".format(input_folder))
os.system("mv *_summary.txt {}../summed_up_annot/".format(input_folder))
