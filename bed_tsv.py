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
    contador = 1
    if os.path.isdir(input_dir_file):
        files = os.listdir(input_dir_file)
        length = len(files)
        for bed_file in files:
            print "Annotating {}\t{} of {}".format(bed_file, contador, length)
            contador += 1
            try:
                os.system("bedtools intersect -a "+input_dir_file+"/"+bed_file+" -b "+annotation_file+" -wa -wb > "+"annoted_"+bed_file[len(bed_file)-14:len(bed_file)-4]+".txt")
            except: 
                print "bedtools is not installed, please use: sudo apt-get install bedtools"
    else:
        os.system("bedtools intersect -a "+input_dir_file+" -b "+annotation_file+" -wa -wb > "+"annoted_"+input_dir_file[len(input_dir_file)-14:len(input_dir_file)-4]+".txt")
    os.system("mkdir "+output_folder)
    os.system("mkdir "+output_folder+"/annotated")
    os.system("mv annoted_* "+output_folder+"/annotated")
    os.system("ls -1 "+output_folder+"/annotated/ | wc -l")
    os.system("find "+output_folder+"/annotated/ -name 'annoted_*' | xargs wc -l")
    print "---Done annotation---"
    return output_folder+"/annotated"

def tsver(input_dir_file):
    contador = 1
    paths, input_files = input_checker(input_dir_file)
    
    if os.path.isfile(input_files[0]):
       print "Single annotation mode with {}".format(os.path.basename(input_files[0]))
       path = paths[0] 
       annotation_handle = open(input_files[0],'r')
    else:
        for path in paths:
            if re.match('.*annotated/*', path) is not None:
                input_files = os.listdir(path)
                break
                
    to_process = len(input_files)
    os.system("mkdir "+path+"/../tsvs")
    if len(input_files) <= 1:
        tsv_file_name = os.path.splitext(os.path.basename(input_files[0]))[0]
        print "Transforming {} to tsv\t{} of {}".format(input_files[0], contador, to_process)
        tsv_file = open(path+"../tsvs/"+tsv_file_name+".tsv", "w")
        for line in file_handle:
            bytabs = line.split("\t")
            if len(bytabs) >= 14:
                annot_str = re.sub('\=|\;|\:|\,', '\t', bytabs[-1])
                string = "\t".join(bytabs[0:-1])+"\t"+annot_str
                tsv_file.write(string)
        tsv_file.close()
        annotation_handle.close()
    else:
        for annotated in input_files:
            print "Transforming {} to tsv\t{} of {}".format(annotated, contador, to_process)
            contador += 1
            tsv_file_name = os.path.splitext(os.path.basename(annotated))[0]
            tsv_file = open(path+"../tsvs/"+tsv_file_name+".tsv", "w")
            annotation_handle = open(path+annotated,'r')
            for line in annotation_handle:
                bytabs = line.split("\t")
                if len(bytabs) >= 14:
                    annot_str = re.sub('\=|\;|\:|\,', '\t', bytabs[-1])
                    string = "\t".join(bytabs[0:-1])+"\t"+annot_str
                    tsv_file.write(string)
        tsv_file.close()
        annotation_handle.close()
    os.system("ls -1 "+path+"/../tsvs/ | wc -l")
    os.system("find "+path+"/../tsvs/ -name '*.tsv' | xargs wc -l")
    print "---Done gene tsv outlay---\n\
    tsvs saved at: {}".format(os.path.realpath(path+"/../tsvs/"))
    return path+"/../tsvs"
