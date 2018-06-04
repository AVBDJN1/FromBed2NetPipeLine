#!bin/usr/python3
# -*- coding: utf-8 -*-.
# Adrian Garcia Moreno

import io
import re
import os 
import subprocess
import json
from decimal import *
from input_controller import *
import argparse

def getname(filein):
    name = os.path.splitext(filein)[0].split("_")[-1]
    return name

def get_dic(path_to_json):
    with open(path_to_json, "r") as jsonhndl:
        dic = json.load(jsonhndl)
    return dic

def get_value(key, dictionary, errors):
    if dictionary == None:
        return ""
    else:
        try:    
            value = dictionary[key]
        except:
            errors.append(key)
            value = ""
        return value, errors

def loci_gens_reader(annotpath, name, cutoff, score_col = 4, gen_col = 9):
    cut_off = Decimal(cutoff).quantize(Decimal('.01'))
    loci_gens = {}
    annotated = "{}annoted_{}.txt".format(annotpath, name)
    annotated = open(annotated, "r")
    for line in annotated:
        info = line.rstrip().split("\t")
        score = Decimal(info[score_col]).quantize(Decimal('.01'))
        if score < cutoff:
            continue
        locus = ":".join(info[0:3])
        gen = info[gen_col]
        if locus in loci_gens:
            loci_gens[locus].add(gen)
        else:
            loci_gens[locus] = set([gen])
    annotated.close()
    return(loci_gens) #{locus:{gens}, ...}
    
def loci_gens_writer(loci_gens):
    string = ""
    for locus in loci_gens:
        gens = ", ".join(loci_gens[locus])
        string += "{}\n{}\n\n".format(locus, gens)
    return(string)
    
def urlizer(text):
    lines = text.split("\n")
    
    geneID_urlbase = '<a href= "https://www.ncbi.nlm.nih.gov/gene/{}" style="color:#5CD7FF" target=new>{}</a>'
    regions_urlbase = '<a href= "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position={}%3A{}-{}" target=new>{}</a>'
    #GO_urlbase = '<a href= "http://amigo.geneontology.org/amigo/term/{}" target=new>{}</a>'
    #HP_urlbase = '<h1><a href= "http://www.human-phenotype-ontology.org/hpoweb?id={}" target=new>{}</a></h1>'
    
    GeneID_match = re.compile(r"^(\d+)")
    regions_match = re.compile(r"^(chr\S*)")
    #GOID_match = re.compile(r"^(GO:\d+)")
    #HP_match = re.compile(r"(feature:[0-9]+)")
    
    for numline in range(len(lines)):
        if regions_match.search(lines[numline]) is not None:
            sor_info = lines[numline].split("\t")
            chro, start, end = sor_info[0].split(":")
            sor_info[0] = regions_urlbase.format(chro,start,end, sor_info[0])
            lines[numline] = "\t".join(sor_info)
        
        if GeneID_match.search(lines[numline]) is not None:
            lines[numline] = ", ".join([geneID_urlbase.format(gene, gene) for gene in lines[numline].split(", ")])        

        # if GOID_match.search(lines[numline]) is not None:
            # info = lines[numline].split("\t")
            # GOID = info[0]
            # genes = info[-1]
            # info[0] = GO_urlbase.format(GOID, GOID)
            # info[-1] = ":".join([geneID_urlbase.format(gene, gene) for gene in genes.split(":")])
            # lines[numline] = "\t".join(info)
            
        #if HP_match.search(lines[numline]) is not None:
        #    info = lines[numline].split("\t")
        #    info[0] = HP_urlbase.format(info[0], info[0])
        #    lines[numline] = "\t".join(info[0])
            
    return("\n".join(lines))

def GO_writer(GOanalysis_folder, onts, name, GO_cutoff=0.01, score_col = 5):
    GOanalysis_str = ""
    for ont in onts:
        GOanalysis = "{}{}_{}.txt".format(GOanalysis_folder, ont, name)
        GOanalysis_txt = open(GOanalysis, "r").read().rstrip("\n").split("\n")
        GOanalysis_txt[0] = "Ontology == {}\n{}".format(ont, GOanalysis_txt[0])
        for line in range(len(GOanalysis_txt[1:])):
            score = GOanalysis_txt[line+1].split("\t")[score_col].split(" ")[-1]
            if float(score) > GO_cutoff:
                break
        GOanalysis_str += "{}\n\n".format("\n".join(GOanalysis_txt[:line+1]))
    return (GOanalysis_str.rstrip("\n\n"))

def GO_writer_html(GOanalysis_folder, onts, name, GO_cutoff=0.01, score_col = 5):
    GO_urlbase = get_html_bases("GO")
    geneID_urlbase = get_html_bases("EntrezID")
    GOanalysis_str = ""
    for ont in onts:
        GOanalysis = "{}{}_{}.txt".format(GOanalysis_folder, ont, name)
        GOanalysis_txt = open(GOanalysis, "r").read().rstrip("\n").split("\n")
        GOanalysis_txt[0] = "Ontology == {}\n{}".format(ont, GOanalysis_txt[0])
        for line in range(len(GOanalysis_txt[1:])):
            info = GOanalysis_txt[line+1].split("\t")
            if float(info[score_col].split(" ")[-1]) > GO_cutoff:
                break
            info[0] = GO_urlbase.format(info[0], info[0])
            info[-1] = ":".join([geneID_urlbase.format(gene, gene) for gene in info[-1].split(":")])
            GOanalysis_txt[line+1] = "\t".join(info)
        GOanalysis_str += "{}\n\n".format(html_tabloider("\n".join(GOanalysis_txt[:line+1])))
    return (GOanalysis_str.rstrip("\n\n"))
    
def html_tabloider(go_result_table):
    lines = go_result_table.split("\n")
    lines[0] =  '<table>\n<caption style="text-align:left"><h3>{}</h3></caption>'.format(lines[0])
    lines[1] = "<tr><th>{}</th></tr>".format("</th><th>".join(lines[1].split("\t")))
    for line in range(len(lines[2:])):
        lines[line+2] = "<tr><td>{}</td></tr>".format("</td><td>".join(lines[line+2].split("\t")))
    lines.append("</table>")
    return "\n".join(lines)

def get_html_bases(nomenclature, nomen_urls_dict = None):
    # This dict should be extracted from its proper json, 
    nomen_urls_dict = {'HPO':'<h1><a href= "http://www.human-phenotype-ontology.org/hpoweb?id={}" target=new>{} {}</a></h1>\n',
                       'GO':'<a href= "http://amigo.geneontology.org/amigo/term/{}" target=new>{}</a>',
                       'EntrezID':'<a href= "https://www.ncbi.nlm.nih.gov/gene/{}" style="color:#5CD7FF" target=new>{}</a>'}
    # Many DB offer an easy way to explore their info with base URLs
    # that only differ in the biological element to be explored
    # it is interesting to select which url base use according to the 
    # data under exploration, a function will be created and for
    # example the urls base will be saved in a JSON to extract them
    # This function could also work with the nomenclature identifier
    return nomen_urls_dict[nomenclature]

def summary_n_htmlizer(GOanalysis_folder, onts, GOcutoff, annotpath,
            locusCutoff, outputFolder, Feature_names, to_plaintxt):
                
    features = set([getname(feature) for feature in os.listdir(GOanalysis_folder)])
    Feature_names = get_dic(Feature_names)
    os.makedirs(outputFolder, exist_ok=True)
    
    to_process = len(features) 
    count = 1
    feature_name_errors = []
    
    if to_plaintxt:
        for feature in features:
            summhdnl = open("{}/{}.txt".format(outputFolder, feature), "w")
            feature_name, feature_name_errors = get_value(feature, Feature_names, feature_name_errors)
            summhdnl.write("{} {}\n".format(feature, feature_name))
            loci_gens = loci_gens_reader(annotpath, feature, locusCutoff)
            summhdnl.write("{}".format(loci_gens_writer(loci_gens)))
            summhdnl.write("{}".format(GO_writer(GOanalysis_folder, onts, feature, GOcutoff)))
            summhdnl.close()
    else:    
        HTML_heading = "<!DOCTYPE html>\n<html>\n<head>\n<title>{}</title>\
                        \n<style>\ntable, th, td, capth {{border:1px solid \
                        black; border-collapse: collapse; padding: 5px; \
                        vertical-align: top; text-align: left;}}\n</style>\
                        </head>\n<body>\n<pre>\n"
        HPO_urlbase = get_html_bases("HPO") ## This function can be expanded
        HTML_footing = "</pre>\n</body>\n</html>"
    
        for feature in features:
            summhdnl = open("{}/{}.html".format(outputFolder, feature), "w")
            summhdnl.write(HTML_heading.format(feature))
            feature_name, feature_name_errors = get_value(feature, Feature_names, feature_name_errors)
            summhdnl.write(HPO_urlbase.format(feature, feature, feature_name))
            loci_gens = loci_gens_reader(annotpath, feature, locusCutoff)
            loci_gens = loci_gens_writer(loci_gens)
            loci_gens = urlizer(loci_gens)
            GOanalysis_str = GO_writer_html(GOanalysis_folder, onts, feature, GOcutoff)
            summhdnl.write("{}{}\n{}".format(loci_gens, GOanalysis_str, HTML_footing))
            summhdnl.close()
    
    print("Features without name {}:\n{}".format(len(feature_name_errors),
    "\n".join(feature_name_errors)))

def argparser():         
    parser = argparse.ArgumentParser(  
        prog = 'Summary and HTMLizer',
        prefix_chars = '-',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = '''This summs up, by default in a html file, unless
        plain text is desired, the results of the pipeline''')
    
    parser.add_argument('-annot','--annotFold', 
    help='The folder with the annotated bed files', nargs='?')
    
    parser.add_argument('-GO','--GOResFold', 
    help='The folder with the GO annalysis resulting files', nargs='?')
    
    parser.add_argument('-LOCcut','--locusCutoff', nargs='?',
    help='Determine the minimun valor of the chromosomic regions to filter \
    the results ((values: between) - default: %(default)s)', type=float, default=0.0)
    
    parser.add_argument('-GOcut','--GOcutoff', nargs='?',
    help='Determine the minimun valor of the GO p-value to filter the\
     results ((values: between) - default: %(default)s)', default=0.05, type=float)
    
    parser.add_argument('-onts','--ontologies', 
    help='Ontologies that you want to summup separated by "," (comma)\
    [BP,CC,MF]', nargs='*', default = ['BP', 'CC', 'MF'], choices=['CC,BP', 
    'BP,CC','MF,BP','BP,MF','CC,MF','MF,CC','BP','CC', 'MF'])
    
    parser.add_argument('-Fnames','--Feature_names', 
    help='JSON file with the data organised as "Feature_ID":"Feature_name"\
    if none is provided the feature_name will be not included', nargs='?')
    
    parser.add_argument('-o','--outputfolder', 
    help='The name of the folder where your summary files will be dumped', nargs='?')
    
    parser.add_argument('-to_txt','--to_plaintxt', 
    help='Flag to left the summary in plain text instead of HTML', 
    action="store_true")
    
    args = parser.parse_args()
    return args

args = argparser()

GOResFold = args.GOResFold
onts = args.ontologies
GOcutoff = args.GOcutoff
annotFold = args.annotFold
locusCutoff = args.locusCutoff
outputFolder = args.outputfolder
Feature_names = args.Feature_names
to_plaintxt = args.to_plaintxt

summary_n_htmlizer(GOResFold, onts, GOcutoff, annotFold, locusCutoff, outputFolder, Feature_names, to_plaintxt)
