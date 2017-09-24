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
    # in case that bedtools is not installed use: apt-get install bedtools
    contador = 1
    if os.path.isdir(input_dir_file):
        files = os.listdir(input_dir_file)
        length = len(files)
        for bed_file in files:
            print "Annotating {}\t{} of {}".format(bed_file, contador, length)
            contador += 1
            os.system("bedtools intersect -a "+input_dir_file+"/"+bed_file+" -b "+annotation_file+" -wa -wb > "+"annoted_"+bed_file[len(bed_file)-14:len(bed_file)-4]+".txt")
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
# This forcing can be problematic if there is no a folder with that name... 
# but I will trust that people do not change my nomenclature
    if not os.path.isfile(input_files[0]):
        for path in paths:
            if re.match('.*annotated/', path) is not None:
                input_files = os.listdir(path)
                break
    to_process = len(input_files)
    os.system("mkdir "+path+"/../tsvs")
    for annotated in input_files:
        print "Transforming {} to tsv\t{} of {}".format(annotated, contador, to_process)
        contador += 1
        tsv_file_name = os.path.splitext(os.path.basename(annotated))[0]
        tsv_file = open(path+"../tsvs/"+tsv_file_name+".tsv", "w")
        with io.open(path+annotated,'r',encoding="utf-8") as file_handle:
            for line in file_handle:
                bytabs = line.split("\t")
                if len(bytabs) >= 14:
                    annot_str = re.sub('\=|\;|\:|\,', '\t', bytabs[-1])
                    string = "\t".join(bytabs[0:-1])+"\t"+annot_str
                    tsv_file.write(string)
        tsv_file.close()
    os.system("ls -1 "+path+"/../tsvs/ | wc -l")
    os.system("find "+path+"/../tsvs/ -name '*.tsv' | xargs wc -l")
    print "---Done gene tsv outlay---\n\
    tsvs saved at: {}".format(os.path.realpath(path+"/../tsvs/"))
    return path+"/../tsvs"
