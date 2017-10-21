#!bin/usr/python
# -*- coding: utf-8 -*-.
# Adrian Garcia Moreno

import io
import re as re
import os 
import subprocess
import sys

def score_bed_selector(input_file, col_of_score, threshold):
    
    HPs = open(input_file, "r")
    name = os.path.splitext(os.path.basename(input_file))[0]
    output = "{}_{}-up.bed".format(name, threshold)
    selected_HPs = open(output, "w")
    for HP in HPs:
        if len(HP.split("\t")) > 2 and float(HP.split("\t")[int(col_of_score)-1]) >= float(threshold):
            selected_HPs.write(HP)
    selected_HPs.close()
    return output

def criteria_divider(data_file_or_dir, col_of_feature):
    
    output_folder = os.path.splitext(data_file_or_dir)[0]+"_divided"
    os.system("mkdir "+output_folder)
    with io.open(data_file_or_dir, 'r', encoding="utf-8") as file_handle:
        criteria = ""
        for line in file_handle:
            line_columns = line.split("\t")
            if line_columns[int(col_of_feature)-1] != criteria:
                criteria = line_columns[int(col_of_feature)-1]
                output_files = open(output_folder+"/"+criteria+".bed", "w")
                output_files.write(line)
            else:
                output_files.write(line)
        output_files.close()
    return "{}/".format(output_folder)
