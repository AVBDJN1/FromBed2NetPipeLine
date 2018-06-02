#!bin/usr/python
# -*- coding: utf-8 -*-.
# Adrian Garcia Moreno

import io
import re as re
import argparse
import os
import subprocess
import sys

from input_controller import *
from criteria_divider import *
from gff3_manager  import *
from bedtooling import *
from topGOer import *
from summary_n_htmlizer import *
        
def outputter(subfolders, infolder, name):
    if subfolders:
        outfolder = infolder.split("/")
        outfolder[-3] = name
        outfolder = "/".join(outfolder)
    else:
        outfolder = infolder.split("/")
        outfolder[-2] = name
        outfolder = "/".join(outfolder)
    return outfolder

parser = argparse.ArgumentParser(  
    prog='FromBed2[Net]',
    prefix_chars='-',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='''A functional annalysis workflow that
    stars from a bed file and ends with an HTML report
    of the Gene Ontology results''',
    epilog='''Under improvement''',
    )

input_group = parser.add_mutually_exclusive_group(required=True)
input_group.add_argument('-file','--infile', 
help='Bed file/s as input', nargs='?')
input_group.add_argument('-dir','--directory', 
help='Directory/ies with all and only your query bed files', nargs='?',
type=is_directory)

parser.add_argument('-o','--outputfolder', 
help='The name of the folder where your results will be dumped', nargs='?')

### Cdiv
parser.add_argument('-f_clm','--feature_clm', nargs = "?",
help='The column where the features are placed', 
type=int)

parser.add_argument('-s_clm','--score_clm', nargs = "?",
help='The column where the score is placed', 
type=int)

parser.add_argument('-LOCUScut','--LOCUScutoff', nargs='*',
help='Determine the minimun value of with which to filter the \
locus of the initial data acording to the score selected with --score_clm', 
type=float, default = 0.0)

parser.add_argument('-skip', '--skiplines', type=int, default = 0,
help='Number of initial lines (header) in the infile to skip')

parser.add_argument('-d', '--delimiter', default='\t',
help='Delimiter of your columns the default is \\t',
nargs = 1)

### GFF3
gff3_group = parser.add_mutually_exclusive_group(required=True)
gff3_group.add_argument('-annot','--annotationfile', 
help='Annotation file [i.e. gff3] to use it directly with bedtools', 
nargs='?')
gff3_group.add_argument('-down_gff3','--down_gff3', 
help='GFF3 source url (Warning ncbi source might need to be \
transformed from accession to chr location if that is your case (use \
-to_chr). Genecode, here given, is already in chr location base \
from it genes will be extracted to use later with bedtools', 
nargs='?', choices=['gcode_hg19','ncbi_hg19', 'gcode_hg20','ncbi_hg20'])
gff3_group.add_argument('-gff3','--input_gff3', 
help='Gff3 file to extract the genes and use later with with bedtools \
it requieres to specify -consortium', nargs='?')

parser.add_argument('-consortium','--consortium', 
help='If a gff3 file is directly given as input, you must specify the \
consortium which it comes from to extract the genes', choices=["ncbi","gcode"],
default = "ncbi")

parser.add_argument('-to_chr','--to_chr', 
help='Flag to transform RefSeq NC_xxx.x to chromosome \
(Shall I  do an opposite transformation from chr to accesion?)', 
action="store_true")

### topGOer

parser.add_argument('-map','--mapfile', 
help='Mapping file with the genes to GO ontology. Be sure that genes are \
annotated with the same nomenclature than your genes lists', nargs='?')

parser.add_argument('-genCol','--genCol', 
help='The column where the genes are', nargs='?', default = "1")

parser.add_argument('-genSCO','--score_col_cutoff', nargs='?',
help='If given the column and minimun valor of the gene separated by \
"," or "-" (no space), ie 5,2 genes are in column 5 and It will select only \
those with score above 2')

parser.add_argument('-taxid','--taxid', nargs= '?',
help='Taxid to retrieve the desired GO annotation data human\
,9606, is the default')

parser.add_argument('-mode','--mode', 
help='The GO categories to annalyse: BP(biological process)\
CC(cellular component) and MF(mollecular function), separated \
by ","',nargs='?', default = 'BP,CC,MF',choices=['CC,BP', 'BP,CC', 
'MF,BP', 'BP,MF', 'CC,MF', 'MF,CC', 'BP', 'CC', 'MF'])

parser.add_argument('-GOcut','--GOcutoff', nargs='?',
help='Determine the minimun valor of the GO p-value to filter the \
results', default = '0.05')

### Summ n HTML
parser.add_argument('-SGOcut','--SGOcutoff', nargs='?',
help='Determine the minimun valor of the p-value to filter the \
results of the GO enrichment annalysis to be shown in the summary report \
((values: between) - default: %(default)s)', default=0.05, type=float, 
choices=map(lambda x: x/100.0, range(0,105, 5)))

parser.add_argument('-SLOCcut','--SLOCcut', nargs='?',
help='Determine the minimun valor of the chromosomic regions to filter \
the results ((values: between) - default: %(default)s)', type=float, default=0.0)

parser.add_argument('-Fnames','--feature_names', 
help='JSON file with the data organised as "Feature_ID":"Feature_name"\
if none is provided the feature_name will be not included', nargs='?')

parser.add_argument('-to_txt','--to_plaintxt', 
help='Flag to left the summary in plain text instead of HTML', 
action="store_true")

parser.add_argument('-v', '--version', 
action='version', version='%(prog)s 0.9.0')

args = parser.parse_args()

infile = args.infile
directory = args.directory
outputfolder = args.outputfolder

# Cdiv
feature_clm = args.feature_clm - 1
score_clm = args.score_clm - 1
LOCUScutoffs = args.LOCUScutoff
skiplines = args.skiplines
delimiter = args.delimiter
# GFF3
gff3_source  = args.down_gff3
to_chr = args.to_chr
input_gff3 = args.input_gff3
consortium = args.consortium
# BedTooling
annotationfile = args.annotationfile
# topGOer
mapfile = args.mapfile
genCol = args.genCol
score_col_cutoff = args.score_col_cutoff
taxid = args.taxid
mode = args.mode
GOcutoff = args.GOcutoff
# Summ n HTML
onts = mode.split(",")
SGOcut = args.SGOcutoff
SLOCcut = args.SLOCcut
feature_names = args.feature_names
to_plaintxt = args.to_plaintxt
        
rscript_path = os.path.dirname(sys.argv[0])+"/"
rscript_path = "" if rscript_path == "/" else rscript_path

if directory is not None and feature_clm is not None:
    subfolders = True
else:
    subfolders = False

outs = CDiv(infile, directory, feature_clm, score_clm, LOCUScutoffs, 
            delimiter, skiplines, outputfolder)
if annotationfile == None: 
    annotationfile = gff3_manager(down_gff3, to_chr, input_gff3, consortium)

outsfolders = getuniquefolders(outs, 1)
bedouts = set()
for folder in outsfolders:
    outputfolder = outputter(subfolders, folder, "bedtooled")
    bedouts.update([outputfolder])
    bedtooling(outputfolder, None, folder, annotationfile)
GOouts = set()
for folder in bedouts:
    outputfolder = outputter(subfolders, folder, "GOresults")
    GOouts.update([outputfolder])
    topGOer(folder, mapfile, genCol, outputfolder, score_col_cutoff, taxid, mode, GOcutoff, "")
for folder in GOouts:
    GOResFold = outputter(subfolders, folder, "GOresults")
    annotFold = outputter(subfolders, folder, "bedtooled")
    outputFold = outputter(subfolders, folder, "HTMLs")
    summary_n_htmlizer(GOResFold, onts, SGOcut, annotFold, 
                SLOCcut, outputFold, feature_names, to_plaintxt)



