#!bin/usr/python3
# -*- coding: utf-8 -*-.
# Adrian Garcia Moreno

from input_controller import *
import subprocess
import argparse

def annotation(bed_file, annotation_file, output_folder):
    
    name_bed_file = "annoted_{}{}".format(*getFilename(bed_file))
    outfile = "{}{}".format(output_folder, name_bed_file)
    
    bedtoolscommand = "bedtools intersect -a {} -b {} -wa -wb > {}\
    ".format(bed_file, annotation_file, outfile)
    
    subprocess.run(bedtoolscommand, shell=True,
    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    return outfile

def bedtooling(outputfolder, infile, directory, annotationfile):
 
    output_folder = "{}/".format(outputfolder)
    os.makedirs(output_folder, exist_ok=True)
    bed_files = input_checker(infile, directory)

    counter = 0
    total_elapsed_time = 0
    numofFiles = len(bed_files)
    outputs = set()
    for bed_file in bed_files:
        counter, total_elapsed_time = verbositier_n_timer(bed_file, 
                                numofFiles, counter, total_elapsed_time)
        outfile = annotation(bed_file, annotationfile, output_folder)
        outputs.update([outfile])
    return outputs

# ~ def argparser():
    # ~ parser = argparse.ArgumentParser(  
        # ~ prog='Bedtooling',
        # ~ prefix_chars='-',
        # ~ formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        # ~ description='''A wrapper of bedtools intersect conserving
        # ~ all the information of both files. Please be sure that you
        # ~ hve bedtools installed on your computer.''',
        # ~ )
    
    # ~ input_group = parser.add_mutually_exclusive_group(required=True)
    
    # ~ input_group.add_argument('-f','--infile', 
    # ~ help='Bed file/s as input', nargs='?', type=open)
    
    # ~ input_group.add_argument('-d','--directory', 
    # ~ help='Directory/ies with all and only your query bed files', nargs='?',
    # ~ type=is_directory)
    
    # ~ parser.add_argument('-annot','--annotationfile', 
    # ~ help='Annotation file/s [i.e. gff3] from where to obtain the extract\
    # ~ the desired information', nargs='?')
    
    # ~ parser.add_argument('-o','--outputfolder', 
    # ~ help='The name of the folder where your results will be dumped', nargs='?')
    
    # ~ args = parser.parse_args()
    # ~ return args

# ~ args = argparser()
# ~ outputfolder = args.outputfolder 
# ~ infile = args.infile 
# ~ directory = args.directory 
# ~ annotationfile = args.annotationfile    
# ~ bedtooling(outputfolder, infile, directory, annotationfile)
