#!bin/usr/python3
# -*- coding: utf-8 -*-.
# Adrian Garcia Moreno

import os 
import argparse
from input_controller import *

def score_selector(line, delimiter, score_clm, cutoff):
    if float(line.rstrip().split(delimiter)[score_clm]) >= cutoff:            
        return True

def feature_divider(line, delimiter, feature_clm, feature, outputFolder, subfolders, out_hndl, inputFile):
    if line.rstrip().split(delimiter)[feature_clm] != feature:
        feature = line.rstrip().split(delimiter)[feature_clm]
        out_hndl = createOutputFile(outputFolder, inputFile, feature, "", 
                   subfolders, writingMode = "a")
                   
    out_hndl.write(line)
    
    return feature, out_hndl

def criteria_divider(infile, directory, feature_clm, score_clm, cutoff, 
                     delimiter, skiplines, outputFolder):
                         
    if directory is not None and feature_clm is not None:
        subfolders = True
    else:
        subfolders = False

    inputs = input_checker(infile, directory)
    outputs = set()
    if score_clm != None and feature_clm != None:
        outputFolder = "{}{}-cutted-off/Cdivided/".format(outputFolder, cutoff)
        for inpt in inputs:
            input_handle = openFile(inpt, skiplines)
            feature = ""
            out_hndl = ""
            for line in input_handle:
                if line == "":
                    continue
                if score_selector(line, delimiter, score_clm, cutoff):
                    feature, out_hndl = feature_divider(line, delimiter, 
                                        feature_clm, feature, outputFolder, 
                                        subfolders, out_hndl, inpt)
                    outputs.update([out_hndl.name])
            if type(out_hndl) != str:         
                out_hndl.close()
            input_handle.close()
        return outputs
    
    elif args.feature_clm is not None:    
        for inpt in inputs:
            input_handle = openFile(inpt, skiplines)
            feature = ""
            out_hndl = ""
            for line in input_handle:
                if line == "":
                    continue
                feature, out_hndl = feature_divider(line, delimiter, 
                                    feature_clm, feature, outputFolder, 
                                    subfolders, out_hndl, inpt)
                outputs.update([out_hndl.name])
            out_hndl.close()
            input_handle.close()
        return outputs
        
    else:
        outputFolder = "{}({}-cutted-off)/Cdivided/".format(outputFolder, cutoff)
        for inpt in inputs:
            input_handle = openFile(inpt, skiplines)
            out_hndl = createOutputFile(outputFolder, inpt, inpt, "", subfolders)
            for line in input_handle:
                if line == "":
                    continue
                if score_selector(line, delimiter, score_clm, cutoff):
                    out_hndl.write(line)
            outputs.update([out_hndl.name])                    
            out_hndl.close()
            input_handle.close()
        return outputs

def CDiv(infile, directory, feature_clm, score_clm, cutoff, delimiter, skiplines, outputFolder):
    outs = set()
    for cuts in cutoff:
        outputs = criteria_divider(infile, directory, feature_clm, score_clm, 
                         cuts, delimiter, skiplines, outputFolder)
        outs.update(outputs)
    return outs

# ~ def argparser():
    # ~ parser = argparse.ArgumentParser(  
        # ~ prog = 'Score and Criteria Divider',
        # ~ prefix_chars = '-',
        # ~ formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        # ~ description = '''To divide in different files the lines that have 
        # ~ a specific string in one of their column and that surpasss the 
        # ~ cutoff''')
    
    # ~ input_group = parser.add_mutually_exclusive_group(required=True)
    
    # ~ input_group.add_argument('-file','--infile', 
    # ~ help='Csv-like file/s as input', nargs='?')
    
    # ~ input_group.add_argument('-dir','--directory', 
    # ~ help='Directory/ies with all and only your csv-like files', nargs='?',
    # ~ type=is_directory)
    
    # ~ parser.add_argument('-skip', '--skiplines', type=int, default = 0,
    # ~ help='Number of initial lines (header) in the infile to skip')
    
    # ~ parser.add_argument('-d', '--delimiter', default='\t',
    # ~ help='Delimiter of your columns the default is \\t',
    # ~ nargs = 1)
    
    # ~ parser.add_argument('-f_clm','--feature_clm', nargs = "?",
    # ~ help='The column where the features are placed', 
    # ~ type=int)
    
    # ~ parser.add_argument('-s_clm','--score_clm', nargs = "?",
    # ~ help='The column where the score is placed', 
    # ~ type=int)
    
    # ~ parser.add_argument('-cuts','--cutoffs', nargs='*',
    # ~ help='Determine the minimun value or values of with which to filter the data\
    # ~ acording to the score selected with --score_clm', 
    # ~ type=float, default = 0.0)
    
    # ~ parser.add_argument('-o','--outputfolder', 
    # ~ help='The name of the folder where your results will be dumped', nargs='?',
    # ~ default = "./")
    
    # ~ args = parser.parse_args()
    # ~ return args

# ~ args = argparser()
# ~ infile = args.infile
# ~ directory = args.directory
# ~ feature_clm = args.feature_clm - 1 if args.feature_clm != None else None
# ~ score_clm = args.score_clm - 1 if args.score_clm != None else None
# ~ delimiter = args.delimiter
# ~ directory = args.directory
# ~ cutoffs = args.cutoffs
# ~ skiplines = args.skiplines
# ~ outputFolder = args.outputfolder+"/" if args.outputfolder[-1] != "/" else args.outputfolder

# ~ CDiv(infile, directory, feature_clm, score_clm, cutoffs, delimiter, skiplines, outputFolder)
