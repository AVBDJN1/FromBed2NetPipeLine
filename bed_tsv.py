#!bin/usr/python
# -*- coding: utf-8 -*-.
# Adrian Garcia Moreno

from input_checker import input_checker
import io
import re as re
import os 
import subprocess
import sys

def annotation(input_dir_file, annotation_file, output_folder):

    try:
        os.makedirs("{}/annotated".format(output_folder))
    except:
        answer = raw_input("Folder already exists, some files could be overwritten\nDo you want to continue? (Y/N)")
        if answer == "N":
            sys.exit()
        
    counter = 0

    paths, input_files = input_checker(input_dir_file)
    
    if os.path.isfile(input_files[0]) and len(input_files) == 1:
        path = paths[0] 
    else:
        input_files = [paths[0]+bed_file for bed_file in input_files]

    to_process = len(input_files)
    for bed_file in input_files:
        counter += 1
        name_bed_file = os.path.basename(bed_file)
        print "Annotating {}\t{} of {}".format(name_bed_file, counter, to_process)
        try:
            os.system("bedtools intersect -a {} -b {} -wa -wb > {}/annotated/annotated_{}".format(bed_file, annotation_file, output_folder, name_bed_file))
        except: 
            print "bedtools is not installed, please use:\n    sudo apt-get install bedtools"

    os.system("ls -1 {}/annotated/ | wc -l".format(output_folder))
    os.system("find {}/annotated/ -name 'annotated_*' | xargs wc -l".format(output_folder))
    print "---Done annotation---"
    return "{}/annotated/".format(output_folder)
# return "{}/annotated/".format(output_folder)

def tsver(input_dir_file):
    #This is functional only for gff3 files bynow
    counter = 1
    paths, input_files = input_checker(input_dir_file)
    
    if os.path.isfile(input_files[0]) and len(input_files) == 1:
        path = paths[0] 
    else:
        for path in paths:
            if re.match('.*annotated/*', path) is not None:
                input_files = [path+annotated for annotated in os.listdir(path)]
                break
                
    to_process = len(input_files)
    try:
        os.mkdir("{}/../tsvs".format(path))
    except:
        answer = raw_input("Folder already exists, some files could be overwritten\nDo you want to continue? (Y/N)")
        if answer == "N":
            sys.exit()
                    
    for annotated in input_files:
        print "Transforming {} to tsv\t{} of {}".format(annotated, counter, to_process)
        counter += 1
        
        tsv_file_name = os.path.splitext(os.path.basename(annotated))[0]    
        tsv_file = open("{}../tsvs/{}.tsv".format(path, tsv_file_name), "w")
        
        annotation_handle = open(annotated,'r')
        for line in annotation_handle:
            bytabs = line.split("\t")
            if len(bytabs) >= 3:
                annot_str = re.sub('\=|\;|\:|\,', '\t', bytabs[-1])
                string = "\t".join(bytabs[0:-1])+"\t"+annot_str
                tsv_file.write(string)
        tsv_file.close()
        annotation_handle.close()
        
    os.system("ls -1 {}../tsvs/ | wc -l".format(path))
    os.system("find {}../tsvs/ -name '*.tsv' | xargs wc -l".format(path))
    
    print "---Done gene tsv outlay---\n\
    tsvs saved at: {}".format(os.path.realpath("{}/../tsvs".format(path)))
    return "{}/../tsvs".format(path)
# return "{}/../tsvs".format(path)
