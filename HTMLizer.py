#!bin/usr/python
# -*- coding: utf-8 -*-.
# Adrian Garcia Moreno

# HTMLizer

import io
import re
import os 
import subprocess
import sys
from input_checker import input_checker
                
def tabloider(go_result_table, summary_warning = ""):
    handle = go_result_table.split("\n")
    go_result_table_htmlizer = ""
    for line in handle:
        if line[0] != "#":
            bytabs = line.split("\t")
            bytabs[0] = "<tr><td>{}".format(bytabs[0])
            bytabs[-1] = "{}</td></tr>\n".format(bytabs[-1])
            go_result_table_htmlizer += "</td><td>".join(bytabs)
        else:
            go_result_table_htmlizer += line+"\n"
    if summary_warning != "":
        summary_warning = "<tr><td colspan=6>{}</td></tr>\n".format(summary_warning)
    go_result_table_htmlizer = "{}{}</table>".format(go_result_table_htmlizer, summary_warning)
    return go_result_table_htmlizer
    
def urlizer(full_text, input_folder_or_file, GOus_dict):
    regions_urlbase = '<a href= "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position={}%3A{}-{}" target=new>{}</a>'#.format(chro,start,end, region)
    HP_urlbase = '<h1><a href= "http://www.human-phenotype-ontology.org/hpoweb?id={}" target=new>{}</a></h1>'#.format(HP, HP)
    GO_urlbase = '<a href= "http://amigo.geneontology.org/amigo/term/{}" target=new>{}</a>'#.format(GOID, GOID)
    geneID_urlbase = '<a href= "https://www.ncbi.nlm.nih.gov/gene/{}" style="color:#7a1919;background-color:#ffffa0" target=new>{}</a>'#.format(GeneID, GeneID)
    full_graph_base = '<a href= "http://amigo.geneontology.org/visualize?format=png&term_data_type=string&mode=amigo&term_data={}" target=new>AmiGO GRAPH</a>'#.format(GOus_dict[lines[line][3]])
    caption_base = '<table>\n<caption>{} {} <a href=..{} target=new>TogGO GRAPH</a></caption>\n'#.format(full_graph, lines[line], full_path_name)
    regions_match = re.compile(r"^chr.+")
    GeneID_match = re.compile(r"^(\d+)")
    caption_match = re.compile(r"^###+.*")
    HP_match = re.compile(r"(HP:[0-9]+)")
    GOID_match = re.compile(r"(GO:[0-9]+).*")
    
    lines = full_text.split("\n")
    
    for line in xrange(len(lines)):
        
        if regions_match.search(lines[line]) is not None:
            chro, start, end, score = lines[line].split("\t")
            region = "\t".join([chro, start, end])
            lines[line] = regions_urlbase.format(chro,start,end, region)+"\t"+score
            
        elif caption_match.search(lines[line]) is not None:
            full_path_name = subprocess.check_output("find {}../*GOanalysis/{}/{}*.pdf".format(input_folder_or_file, feature_name, lines[line][3]), shell=True).rstrip("\n")
            full_path_name = full_path_name.split("..")[-1] # this is to normalise the path and later i will add the .. in the substitution
            full_graph = full_graph_base.format(GOus_dict[lines[line][3]])
            lines[line] = caption_base.format(full_graph, lines[line], full_path_name)

        elif GeneID_match.match(lines[line]) is not None:
            GeneID = GeneID_match.match(lines[line]).group(1)
            lines[line] = geneID_urlbase.format(GeneID,GeneID)
            
        elif HP_match.search(lines[line]) is not None:
            HP = HP_match.search(lines[line]).group(1)
            lines[line] = re.sub(HP, HP_urlbase.format(HP, HP), lines[line])
            
        elif GOID_match.search(lines[line]) is not None:
            GOID = GOID_match.search(lines[line]).group(1)
            lines[line] = re.sub(GOID, GO_urlbase.format(GOID, GOID), lines[line])
        
    full_text = "\n".join(lines)
    return full_text
    
input_folder_or_file = sys.argv[1]

try:
    os.mkdir("{}../htmls".format(input_folder_or_file))
except:
    answer = raw_input("Folder already exists, some files could be overwritten\nDo you want to continue? (Y/N)")
    if answer == "N":
        sys.exit()

paths, input_files = input_checker(input_folder_or_file)
    
if os.path.isfile(input_files[0]) and len(input_files) == 1:
    path = paths[0] 
else:
    path = paths[0] 
    input_files = [paths[0]+bed_file for bed_file in input_files]

files_toprocess = len(input_files)
counter = 1

for input_summedup in input_files:
    print "HTMLizing {}\t{} of {}".format(input_summedup, counter, files_toprocess)
    counter += 1
    output_html = os.path.splitext(os.path.basename(input_summedup))[0]
    results_handle = open(input_summedup, "r")
    feature_name = results_handle.readline().rstrip()
    results_handle.seek(0)
    full_text = results_handle.read()
    parts = full_text.split("\n\n")
    
    GOus_dict = {}
    GOID_match = re.compile(r"(GO:[0-9]+).*")
    if len(parts) > 1 and parts[-2][0] == "#":
        
        GOus = GOID_match.findall(parts[-4])
        if len(GOus) == 0:
            summary_warning = "No GO.ID passed the summary threshold, you can still check the results of the TopGO GRAPH"
        else:
            summary_warning = ""
        GO_category = parts[-4][3]
        GOus_dict[GO_category] = re.sub(":", "%3A", " ".join(GOus))
        parts[-4] = tabloider(parts[-4], summary_warning)

        GOus = GOID_match.findall(parts[-3])
        if len(GOus) == 0:
            summary_warning = "No GO.ID passed the summary threshold, you can still check the results of the TopGO GRAPH"
        else:
            summary_warning = ""
        GO_category = parts[-3][3]
        GOus_dict[GO_category] = re.sub(":", "%3A", " ".join(GOus))
        parts[-3] = tabloider(parts[-3], summary_warning)
        
        GOus = GOID_match.findall(parts[-2])
        if len(GOus) == 0:
            summary_warning = "No GO.ID passed the summary threshold, you can still check the results of the TopGO GRAPH"
        else:
            summary_warning = ""
        GO_category = parts[-2][3]
        GOus_dict[GO_category] = re.sub(":", "%3A", " ".join(GOus))   
        parts[-2] = tabloider(parts[-2], summary_warning)

    full_text = "\n\n".join(parts)
    #print GOus_dict
    full_text = urlizer(full_text, path, GOus_dict)
    
    HTML_heading = "<!DOCTYPE html>\n<html>\n<head>\n<title>{}</title>\
                    \n<style>\ntable, th, td {{border:1px solid black; border-collapse:\
                    collapse;}}\ncapth, td {{padding: 5px;}}\n</style></head>\n<body>\n<pre>".format(feature_name)
    HTML_footing = "</pre>\n</body>\n</html>"

    htmlized = open("{}../htmls/{}.html".format(path, output_html), "w")
    htmlized.write(HTML_heading)
    htmlized.write(full_text)
    htmlized.write(HTML_footing)
    htmlized.close()

# file:///home/bioinfo/Desktop/TFM/Study/hyI/Study/hyI/hg19_ncbi_run/annotation_results/hg19_enrichment_GOanalysis/results/HP:0100704/MFgraph_HP:0100704_classic_10_all.pdf
# ../Study/hyI/hg19_ncbi_run/summed_up_annotation/


#os.system("mv *.html htmls/")
