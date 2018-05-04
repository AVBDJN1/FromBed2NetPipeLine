#!bin/usr/python
# -*- coding: utf-8 -*-.
# Adrian Garcia Moreno

# HTMLizer

import io
import re
import os 
import subprocess
import sys
import json
from decimal import *


def getname(filein):
    name = os.path.splitext(filein)[0].split("_")[-1]
    return name

def tidy_reader(annotpath, name, cut_off):
    cut_off = Decimal(cut_off).quantize(Decimal('.01'))
    sor_gens = {}
    annotated = "{}annotated_{}.bed".format(annotpath, name)
    annotated = open(annotated, "r")
    for info in annotated:
        info = info.rstrip().split("\t")
        score = Decimal(info[4]).quantize(Decimal('.01'))
        if score < cut_off:
            continue
        sor = ":".join(info[0:3])
        mutation = info[5].split(".")
        mutation = mutation[2] if len(mutation) > 2 else mutation[0]
        mutation_score = "{}={}".format(mutation, score)
        gen = info[9]
        if sor in sor_gens:
            sor_gens[sor]["Genes"].add(gen)
            sor_gens[sor]["Muts"].add(mutation_score)
        else:
            sor_gens[sor] = {"Genes":set([gen]),"Muts":set([mutation_score])}
    annotated.close()
    return(sor_gens) #{SOR:{MUTs:[MUT=SCORs],GENs:[gens]}}
    
def sor_gens_writer(sor_gens):
    string = ""
    for sor in sor_gens:
        muts = " ".join(sor_gens[sor]["Muts"])
        gens = ", ".join(sor_gens[sor]["Genes"])
        string += "{}\t{}\n{}\n\n".format(sor, muts, gens)
    return(string)
    
def GO_writer(GOanalysis_folder, onts, name):
    GOanalysis_str = ""
    for ont in onts:
        GOanalysis = "{}{}_{}.txt".format(GOanalysis_folder, ont, name)
        GOanalysis_str += "Ontology == {}\n{}\n\n".format(ont,open(GOanalysis, "r").read().rstrip("\n"))
    return (GOanalysis_str.rstrip("\n\n"))

def get_dic(path_to_json):
    with open(path_to_json, "r") as jsonhndl:
        dic = json.load(jsonhndl)
    return dic
    
def urlizer(text):
    lines = text.split("\n")
    geneID_urlbase = '<a href= "https://www.ncbi.nlm.nih.gov/gene/{}" style="color:#5CD7FF" target=new>{}</a>'
    regions_urlbase = '<a href= "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position={}%3A{}-{}" target=new>{}</a>'
    #GO_urlbase = '<a href= "http://amigo.geneontology.org/amigo/term/{}" target=new>{}</a>'
    #HP_urlbase = '<h1><a href= "http://www.human-phenotype-ontology.org/hpoweb?id={}" target=new>{}</a></h1>'
    
    GeneID_match = re.compile(r"^(\d+)")
    regions_match = re.compile(r"^(chr\S*)")
    #GOID_match = re.compile(r"^(GO:\d+)")
    #HP_match = re.compile(r"(HP:[0-9]+)")
    
    for numline in xrange(len(lines)):
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
    
def tabloider(go_result_table):
    handle = go_result_table.split("\n")
    handle[0] =  '<table>\n<caption style="text-align:left"><h3>{}</h3></caption>'.format(handle[0])
    handle[1] = "<tr><th>{}</th></tr>".format("</th><th>".join(handle[1].split("\t")))
    for line in xrange(len(handle[2:])):
        handle[line+2] = "<tr><td>{}</td></tr>".format("</td><td>".join(handle[line+2].split("\t")))
    handle.append("</table>")
    return "\n".join(handle)
    
def get_genename(gene, GI_NAME_dic):
    if gene != "":
        try:
            genename = GI_NAME_dic[gene]["GeneSymbol"] if GI_NAME_dic[gene]["GeneSymbol"] != "" else gene
            return genename
        except:
            return gene
        
def GO_writer_html(GOanalysis_folder, onts, name, GI_NAME_dic, cutoff=0.01):
    GO_urlbase = '<a href= "http://amigo.geneontology.org/amigo/term/{}" target=new>{}</a>'
    geneID_urlbase = '<a href= "https://www.ncbi.nlm.nih.gov/gene/{}" style="color:#5CD7FF" target=new>{}</a>'
    GOanalysis_str = ""
    for ont in onts:
        GOanalysis = "{}{}_{}.txt".format(GOanalysis_folder, ont, name)
        GOanalysis_txt = open(GOanalysis, "r").read().rstrip("\n").split("\n")
        GOanalysis_txt[0] = "Ontology == {}\n{}".format(ont, GOanalysis_txt[0])
        for line in xrange(len(GOanalysis_txt[1:])):
            info = GOanalysis_txt[line+1].split("\t")
            if float(info[5].split(" ")[-1]) > cutoff:
                break
            info[0] = GO_urlbase.format(info[0], info[0])
            info[-1] = ":".join([geneID_urlbase.format(gene, get_genename(gene, GI_NAME_dic)) for gene in info[-1].split(":")])
            GOanalysis_txt[line+1] = "\t".join(info)
        GOanalysis_str += "{}\n\n".format(tabloider("\n".join(GOanalysis_txt[:line+1])))
    return (GOanalysis_str.rstrip("\n\n"))

## Arguments Needed
annotpath = sys.argv[1] 
GOanalysis_folder = sys.argv[2]
cut_off = sys.argv[3]
onts = sys.argv[4].split(",")
outputfolder = sys.argv[5]
HP_NAME_dic = sys.argv[6]
GI_NAME_dic = sys.argv[7]
#cutoff = sys.argv[8]
## Arguments Needed
# annotpath = "Rarebiosis/data/annotated/"
# GOanalysis_folder = "Rarebiosis/GO/"
# cut_off = "2"
# onts = "CC,BP".split(",")
# outputfolder = "Summaries"

HPs = set([getname(hp) for hp in os.listdir(GOanalysis_folder)])
#existingHPS = set([getname(hp) for hp in os.listdir(outputfolder)])
#HPs = HPs.difference(existingHPS)

HTML_heading = "<!DOCTYPE html>\n<html>\n<head>\n<title>{}</title>\
                \n<style>\ntable, th, td, capth {{border:1px solid black; border-collapse:\
                collapse; padding: 5px; vertical-align: top; text-align: left;}}\n</style></head>\n<body>\n<pre>"
HP_urlbase = '<h1><a href= "http://www.human-phenotype-ontology.org/hpoweb?id={}" target=new>{} {}</a></h1>\n'
HTML_footing = "</pre>\n</body>\n</html>"

HP_NAME_dic = get_dic(HP_NAME_dic)
#HP_NAME_dic["HP:0006008"] = "Oboselete? Unilateral brachydactyly"
#HP_NAME_dic["HP:0004488"] = "Oboselete? Macrocephaly at birth"
#HP_NAME_dic["HP:0010282"] = "Oboselete? Thin lower lip vermilion"
#HP_NAME_dic["HP:0010465"] = "Oboselete? Precocious puberty in females"
#HP_NAME_dic["HP:0004686"] = "Oboselete? Short third metatarsal"
#HP_NAME_dic["HP:0009885"] = "Oboselete? Prenatal short stature"
#HP_NAME_dic["HP:0010672"] = "Oboselete? Abnormality of the third metatarsal bone"
#HP_NAME_dic["HP:0002564"] = "Oboselete? Malformation of the heart and great vessels"
#HP_NAME_dic["HP:0011427"] = "Oboselete? Enlarged fetal cisterna magna"
#HP_NAME_dic["HP:0010823"] = "Oboselete? Ridged cranial sutures"
#HP_NAME_dic["HP:0010800"] = "Oboselete? Absent cupid's bow"
#HP_NAME_dic["HP:0001587"] = "Oboselete? Primary ovarian failure"
GI_NAME_dic = get_dic(GI_NAME_dic)

try:
    os.mkdir(outputfolder)
except:
    ""
    
HPwithoutname = []
    
total = len(HPs) 
count = 1
for HP in HPs:
    print "Processing {}, {} of {}".format(HP, count, total)
    count += 1
    with open("{}/{}.html".format(outputfolder, HP), "w") as summhdnl:
        summhdnl.write(HTML_heading.format(HP))
        try:
            HP_name = HP_NAME_dic[HP]
        except:
            HPwithoutname.append(HP)
            HP_name = ""
        summhdnl.write(HP_urlbase.format(HP, HP, HP_name))
        sor_gens = tidy_reader(annotpath, HP, cut_off)
        sor_gens_str = sor_gens_writer(sor_gens)
        sor_gens_str = urlizer(sor_gens_str)
        GOanalysis_str = GO_writer_html(GOanalysis_folder, onts, HP, GI_NAME_dic)
        summhdnl.write("{}{}\n{}".format(sor_gens_str, GOanalysis_str, HTML_footing))

print "HPs without known name {}".format(len(HPwithoutname))
print ",".join(HPwithoutname)

