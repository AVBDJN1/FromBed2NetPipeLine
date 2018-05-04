#!bin/usr/python3.5
# -*- coding: utf-8 -*-.
# Adrian Garcia Moreno

import sys
import os

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

def openFile(file_handle, skipLines = 0):

    if type(file_handle) == str:
        file_handle = open(file_handle, "r")

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

def getFilename(inputFile):
    if type(inputFile) != str:
        inputFile = inputFile.name
    fileName = os.path.splitext(os.path.basename(inputFile))[0]
    extension = os.path.splitext(os.path.basename(inputFile))[1]
    return fileName, extension







