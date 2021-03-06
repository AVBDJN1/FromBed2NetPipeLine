#!bin/usr/python
# -*- coding: utf-8 -*-.
# Adrian Garcia Moreno

import io
import re as re
import os 
import subprocess
import sys

def input_checker(input_dir_file):
    paths = []
    input_files = []
    if os.path.isdir(input_dir_file):
        dirs_and_files = os.listdir(input_dir_file)
        for dir_or_file in dirs_and_files:
            if os.path.isdir(input_dir_file+dir_or_file):  
                paths.append(input_dir_file+dir_or_file+"/")
            else:
                if input_dir_file not in paths:
                    paths.append(input_dir_file)
            input_files.append(dir_or_file)
    elif os.path.isfile(input_dir_file):
        input_files = [input_dir_file]
        paths = ["{}/".format(os.path.dirname(input_dir_file))]
    else:
        print "No folder nor file found with that name"
        sys.exit()
    return paths, input_files 
    
