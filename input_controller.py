#!bin/usr/python3.5
# -*- coding: utf-8 -*-.
# Adrian Garcia Moreno

import sys
import os
import gzip
import time
import urllib.request
import json

def getFilename(inputFile):
    if type(inputFile) != str:
        inputFile = inputFile.name
    fileName = os.path.splitext(os.path.basename(inputFile))[0]
    extension = os.path.splitext(os.path.basename(inputFile))[1]
    return fileName, extension

def is_directory(path_str):
    if os.path.isdir(path_str):
        return path_str
    else:
        msg = "{} is not a directory".format(path_str)
        raise argparse.ArgumentTypeError(msg)    

def input_checker(infile, directory):
    if infile == None:
        return [directory+File for File in os.listdir(directory)]
    else:
        return [infile]

def input_format_reader(inputfile):
    name, extension = getFilename(inputfile)
    if extension == ".gz":
        file_handle = gzip.open(inputfile, "rt")
    else:
        file_handle = open(inputfile, "r")
    return file_handle

def openFile(file_handle, skipLines = 0):

    if type(file_handle) == str:
        file_handle = input_format_reader(file_handle)

    while skipLines > 0: 
        skipLines -= 1
        file_handle.readline()
    # Maybe it is interesting to check directly where is the first line
    # that offers the csv-like format and return the file_handle there
    return file_handle

def createOutputFile(outputFolder, inputFile, outname = "", prefix = "", 
                   subfolders = False, writingMode = "w"):
    
    filename, extension = getFilename(inputFile)
    
    if subfolders:
        os.makedirs(outputFolder+filename, exist_ok=True)
        output = "{}{}/{}{}{}".format(outputFolder,filename,prefix,outname,extension)
        out_hndl = open(output, writingMode)
        
    else:
        if inputFile == outname:
            outname = filename
        os.makedirs(outputFolder, exist_ok=True)
        output = "{}{}{}{}".format(outputFolder,prefix,outname,extension)
        out_hndl = open(output, writingMode)
        
    return out_hndl

def verbositier_n_timer(infile, numofFiles, counter, total_elapsed_time):
    start_time = time.time()
    single_elapsed_time = time.time() - start_time
    total_elapsed_time += single_elapsed_time
    
    single_elapsed_time = time.strftime("%H:%M:%S", time.gmtime(single_elapsed_time))
    total_elapsed_timestr = time.strftime("%H:%M:%S", time.gmtime(total_elapsed_time))
    
    counter += 1
    
    print("Processing {}\t{} of {} Taken {} of so far {}".format(
    infile, counter, numofFiles, single_elapsed_time, total_elapsed_timestr))
    
    
    
    return counter, total_elapsed_time

def verbositier(infile, numofFiles, counter):
    
    counter += 1
    print("Processing {}\t{} of {}".format(
    infile, counter, numofFiles))
    
    return counter

def reporthook(count, block_size, total_size):
    global start_time
    if count == 0:
        start_time = time.time()
        return
    duration = time.time() - start_time
    progress_size = int(count * block_size)
    speed = int(progress_size / (1024 * duration))
    percent = int(count * block_size * 100 / total_size)
    sys.stdout.write("\r...%d%%, %d MB, %d KB/s, %d seconds passed" %
                    (percent, progress_size / (1024 * 1024), speed, duration))
    sys.stdout.flush()

def save(url, filename):
    urllib.request.urlretrieve(url, filename, reporthook)
    sys.stdout.write("\nDone!\n")
