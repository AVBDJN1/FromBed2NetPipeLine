#!bin/usr/python
# -*- coding: utf-8 -*-.
# Adrian Garcia Moreno

# HTMLizer

import io
import re
import os 
import subprocess
import sys
                
def tabloider(go_result_table):
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
    go_result_table_htmlizer = go_result_table_htmlizer + "</table>"
    return go_result_table_htmlizer
    
def urlizer(full_text, input_folder_or_file, GOus_dict):
    regions_urlbase = '<a href= "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position={}%3A{}-{}" target=new>{}</a>'#.format(chro,start,end, region)
    HP_urlbase = '<h1><a href= "http://www.human-phenotype-ontology.org/hpoweb?id={}" target=new>{}</a></h1>'#.format(HP, HP)
    GO_urlbase = '<a href= "http://amigo.geneontology.org/amigo/term/{}" target=new>{}</a>'#.format(GOID, GOID)
    geneID_urlbase = '<a href= "https://www.ncbi.nlm.nih.gov/gene/{}" style="color:#7a1919;background-color:#ffffa0" target=new>{}</a>'#.format(GeneID, GeneID)
    regions_match = re.compile(r"^chr.+")#(r"(chr[0-9]+|chr[XYM])\t([0-9]+)\t([0-9]+)")
    GeneID_match = re.compile(r"^(\d+)")
    caption_match = re.compile(r"^###+.*")
    HP_match = re.compile(r">+>+(HP:[0-9]+)")
    GOID_match = re.compile(r"(GO:[0-9]+).*")
    lines = full_text.split("\n")
    
    for line in xrange(len(lines)):
        if regions_match.search(lines[line]) is not None:
            chro, start, end, score = lines[line].split("\t")
            region = "\t".join([chro, start, end])
            lines[line] = regions_urlbase.format(chro,start,end, region)+"\t"+score

        if GeneID_match.match(lines[line]) is not None:
            GeneID = GeneID_match.match(lines[line]).group(1)
            lines[line] = geneID_urlbase.format(GeneID,GeneID)

        if caption_match.search(lines[line]) is not None:
            full_path_name = subprocess.check_output("find {}../annotation_results/*GOanalysis/results/{}/{}*.pdf".format(input_folder_or_file, HP, lines[line][3]), shell=True).rstrip("\n")
            full_path_name = full_path_name.split("..")[-1]
            full_graph = '<a href= "http://amigo.geneontology.org/visualize?format=png&term_data_type=string&mode=amigo&term_data={}" target=new>Full GO GRAPH</a>'.format(GOus_dict[lines[line][3]])
            lines[line] = '<table>\n<caption>{} {} <a href=../{} target=new>Minimum GO GRAPH</a></caption>\n'.format(full_graph, lines[line], full_path_name)
            
        if HP_match.search(lines[line]) is not None:
            HP = HP_match.search(lines[line]).group(1)
            lines[line] = re.sub(HP, HP_urlbase.format(HP, HP), lines[line])
            HTML_heading = "<!DOCTYPE html>\n<html>\n<head>\n<title>{}</title>\
            \n<style>\ntable, th, td {{border:1px solid black; border-collapse:\
            collapse;}}\ncapth, td {{padding: 5px;}}\n</style></head>\n<body>\n<pre>".format(HP)

        if GOID_match.search(lines[line]) is not None:
            GOID = GOID_match.search(lines[line]).group(1)
            lines[line] = re.sub(GOID, GO_urlbase.format(GOID, GOID), lines[line])
    full_text = "\n".join(lines)
    return HTML_heading, full_text
    
input_folder_or_file = sys.argv[1]
os.mkdir("{}../htmls".format(input_folder_or_file))

if os.path.isfile(input_folder_or_file):
    summed_up_results = [sys.argv[1]]
else:
    summed_up_results = [sys.argv[1]+summed_up_result for  summed_up_result in os.listdir(sys.argv[1])] # FromBed2NetPipeLine/mygoresults.txt

files_toprocess = len(summed_up_results)
counter = 1

for input_summedup in summed_up_results:
    print "HTMLizing {}\t{} of {}".format(input_summedup, counter, files_toprocess)
    counter += 1
    output_html = os.path.splitext(os.path.basename(input_summedup))[0]
    results_handle = open(input_summedup, "r")
    full_text = results_handle.read()
    parts = full_text.split("\n\n")
    
    GOus_dict = {}
    GOID_match = re.compile(r"(GO:[0-9]+).*")
    if len(parts) > 1 and parts[-2][0] == "#" and len(parts):

        GOus = GOID_match.findall(parts[-4])
        GO_category = parts[-4][3]
        GOus_dict[GO_category] = re.sub(":", "%3A", " ".join(GOus))
        parts[-4] = tabloider(parts[-4])

        GOus = GOID_match.findall(parts[-3])
        GO_category = parts[-3][3]
        GOus_dict[GO_category] = re.sub(":", "%3A", " ".join(GOus))
        parts[-3] = tabloider(parts[-3])
        
        GOus = GOID_match.findall(parts[-2])
        GO_category = parts[-2][3]
        GOus_dict[GO_category] = re.sub(":", "%3A", " ".join(GOus))   
        parts[-2] = tabloider(parts[-2])

    full_text = "\n\n".join(parts)

    HTML_heading, full_text = urlizer(full_text, input_folder_or_file, GOus_dict)
    HTML_footing = "</pre>\n</body>\n</html>"

    htmlized = open("{}../htmls/{}.html".format(input_folder_or_file, output_html), "w")
    htmlized.write(HTML_heading)
    htmlized.write(full_text)
    htmlized.write(HTML_footing)
    htmlized.close()

# file:///home/bioinfo/Desktop/TFM/Study/hyI/Study/hyI/hg19_ncbi_run/annotation_results/hg19_enrichment_GOanalysis/results/HP:0100704/MFgraph_HP:0100704_classic_10_all.pdf
# ../Study/hyI/hg19_ncbi_run/summed_up_annotation/


#os.system("mv *.html htmls/")
