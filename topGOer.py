#!bin/usr/python3
# -*- coding: utf-8 -*-.
# Adrian Garcia Moreno

import argparse
import subprocess
import sys
import os

def topGOer(inp, mapfile, genCol, outputfolder, gen_col_cutoff, taxid, 
            mode, GOcutoff, rscript_path, multiprocess = False):
    if gen_col_cutoff is None:
        gen_col_cutoff = ""
    else:
        gen_col_cutoff = "-score {}".format(gen_col_cutoff)
    if taxid is None:
        taxid = ""
    else:
        taxid = "-taxid {}".format(taxid)
    if mode is None:
        mode = ""
    else:
        mode = "-mode {}".format(mode)
    if GOcutoff is None:
        GOcutoff = ""
    else:
        GOcutoff = "-pval_thres {}".format(GOcutoff)
     
    os.makedirs(outputfolder, exist_ok=True)        
    RscriptCommand = "Rscript {}TopGOer.r {} {} {} {} {} {} {} {}".format(
    rscript_path, inp, mapfile, genCol, outputfolder, gen_col_cutoff, 
    taxid, mode, GOcutoff)
    print (RscriptCommand)
    
    subprocess.run(RscriptCommand, shell=True, 
    stdout=subprocess.PIPE, stderr=subprocess.PIPE)

# ~ def argparser():
    # ~ parser = argparse.ArgumentParser(  
        # ~ prog = 'TopGOer Wrapper',
        # ~ prefix_chars = '-',
        # ~ formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        # ~ description = '''This is just an intermediary of TopGOer.r 
        # ~ as a wrapper to ease the workflow. You will need to have installed,
        # ~ R, Rscript and the package topGO and its dependencies''')
    
    # ~ parser.add_argument('-in','--input', 
    # ~ help='The file or folder with the list of genes', nargs='?')
    
    # ~ parser.add_argument('-map','--mapfile', 
    # ~ help='Mapping file/s with the desired ontology. Be sure that genes are \
    # ~ annotated with the same nomenclature than your genes lists. You could \
    # ~ download, the EntreID-GO mapping file writting "gene2go:download" \
    # ~ or if you already hava that dowloaded write "gene2go" and in both last two \
    # ~ cases determine the taxid to extract those annoatations with -taxid',
     # ~ nargs='?')
    
    # ~ parser.add_argument('-genCol','--genCol', 
    # ~ help='The column where the genes are', nargs='?', default = "1")
    
    # ~ parser.add_argument('-o','--outputfolder', 
    # ~ help='The name of the folder where your summary files will be dumped', nargs='?')
    
    # ~ # Optionals
    # ~ parser.add_argument('-genSCO','--score_col_cutoff', nargs='?',
    # ~ help='If given the column and minimun valor of the gene separated by \
    # ~ , or - (no space), ie 5,2 genes are in column 5 and It will select only those with\
    # ~ score above 2')
    
    # ~ parser.add_argument('-taxid','--taxid', nargs= '?',
    # ~ help='Taxid to retrieve the desired GO annotation data human\
    # ~ ,9606, is the default')
    
    # ~ parser.add_argument('-mode','--mode', 
    # ~ help='The GO categories to annalyse: BP(biological process)\
    # ~ CC(cellular component) and MF(mollecular function), separated \
    # ~ by ","',nargs='?', default = 'BP,CC,MF',choices=['CC,BP','BP,CC', 
    # ~ 'MF,BP','BP,MF','CC,MF','MF,CC','BP','CC','MF'])
    
    # ~ parser.add_argument('-GOcut','--GOcutoff', nargs='?',
    # ~ help='Determine the minimun valor of the GO p-value to filter the\
     # ~ results, default is 0.05')
    
    # ~ args = parser.parse_args()
    # ~ return args

# ~ args = argparser()
# ~ inp = args.input
# ~ mapfile = args.mapfile
# ~ genCol = args.genCol
# ~ outputfolder = args.outputfolder
# ~ gen_col_cutoff = args.gen_col_cutoff
# ~ taxid = args.taxid
# ~ mode = args.mode
# ~ GOcutoff = args.GOcutoff
# ~ topGOer(inp, mapfile, genCol, outputfolder, gen_col_cutoff, taxid, mode, GOcutoff)
    
