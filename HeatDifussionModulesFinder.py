#!bin/usr/python
# -*- coding: utf-8 -*-.
# Adrian Garcia Moreno

ndex = __import__("ndex")
import itertools
import os
import io
import pandas as pd
import requests
import json
import base64
import networkx as nx
from ndex.networkn import NdexGraph
import matplotlib.pyplot as plt
import pygraphviz
import sys
import operator as op

def get_cx_from_filename(filename):
    if os.path.isfile(filename):
        with io.open(filename, 'rU') as file_cx:
            return json.load(file_cx)

def insert_heat(G, query_genes):
    query_nodes_ids_name_dict = {}
    for nodeid in G:
        if G.node[nodeid]['name'] in query_genes:
            G.node[nodeid]['diffusion_input'] = 1.0
            query_nodes_ids_name_dict[nodeid] = G.node[nodeid]['name']
        else:
            G.node[nodeid]['diffusion_input'] = 0.0
    query_genes_not_in_network = set(query_genes).difference(query_nodes_ids_name_dict.values())
    return G, query_nodes_ids_name_dict, query_genes_not_in_network
        
def heated_nodes_stats(graph_results, mode='all'):
    G = ndex.networkn.NdexGraph(cx=graph_results.json()['data'])
    results_df = []
    for nodeid in G:
        if mode == 'all':
            if G.node[nodeid]['diffusion_output_heat'] > 0.0:
                results_df.append([nodeid, G.node[nodeid]['name'], G.node[nodeid]['diffusion_output_heat'], G.node[nodeid]['diffusion_output_rank']]) 
        if isinstance(mode, int):
            if G.node[nodeid]['diffusion_output_rank'] < mode:
                results_df.append([nodeid, G.node[nodeid]['name'], G.node[nodeid]['diffusion_output_heat'], G.node[nodeid]['diffusion_output_rank']])
    results_df = sorted(results_df, key=op.itemgetter(3))
    results_df = pd.DataFrame(results_df)
    results_df.columns = ['NodeID', 'Gene', 'Heat', 'Rank']    
    return results_df
    
def modules_checker(query_node_names, heated_genes_ids, my_graph):
    edges = []
    for nodeid in list(itertools.combinations(heated_genes_ids, 2)):
        if my_graph.has_edge(nodeid[0], nodeid[1]) or my_graph.has_edge(nodeid[1], nodeid[0]):
            edges.append((my_graph.node[nodeid[0]]['name'], my_graph.node[nodeid[1]]['name']))
    modules_graph = nx.Graph()
    modules_graph.add_edges_from(edges)
    modules = []
    for module in list(nx.connected_components(modules_graph)):
        if len(set(query_node_names).intersection(map(type(query_node_names[0]), module))) >= 1:
            modules.append(map(type(query_node_names[0]), module))
    return modules, modules_graph
    
def retrieve_query_genes(gene_list_path):
    with io.open(gene_list_path, "rU") as gene_list_handle:
        df = pd.read_csv(gene_list_handle, delimiter = " ", header = None, names=["genes", "scores"])
        gene_list = map(str, list(set(df["genes"].tolist())))
    return gene_list
    
def draw_hierarquical_graph(graph):
    nx.nx_agraph.write_dot(graph,'test.dot')
    pos = nx.nx_agraph.graphviz_layout(graph, prog='dot')
    nx.draw(graph, pos, with_labels=False, arrows=False, node_size=15)
    i = 1
    for p in pos:  # raise text positions
        l = list(pos[p])
        l[1] += 10 * i
        pos[p] = tuple(l)
        i = i * -1
    nx.draw_networkx_labels(graph, pos, font_size=9, font_color="blue")
    plt.show()
    
def full_workflow(gene_list, output_folder, number_of_heated_nodes, graph, url):
    output_name = "".join(os.path.splitext(os.path.basename(gene_list))[0].split("_")[1:])
    gene_list = retrieve_query_genes(gene_list)
    graph_heated, query_nodes_ids_name_dict, query_genes_not_in_network = insert_heat(graph, gene_list)
    query_node_names = query_nodes_ids_name_dict.values()
    # print "{} = Genes-{}\tHeated-{}\tNo Heated-{}\n".format(output_name, len(gene_list), len(query_node_names), len(query_genes_not_in_network), " ".join(sorted(query_node_names)))
    graph_results = requests.post(url, data=json.dumps(graph_heated.to_cx()), headers={'content-type': 'application/json'})
    results_df = heated_nodes_stats(graph_results, number_of_heated_nodes)
    results_df.to_csv("{}HeatData_{}.tsv".format(output_folder, output_name), sep='\t', encoding='utf-8')
    heated_genes_ids = map(int, results_df['NodeID'])
    modules, modules_graph = modules_checker(query_node_names, heated_genes_ids, graph)
    # draw_hierarquical_graph(modules_graph)
    module_counter = 0
    with io.open("{}Modules_{}.txt".format(output_folder, output_name), "w") as modules_handle:
        for module in modules:
            module_counter += 1
            module = u"#{}\n{}\n".format(module_counter, " ".join(module))
            modules_handle.write(module)

graph_path = sys.argv[1]
graph = get_cx_from_filename(graph_path)
graph = ndex.networkn.NdexGraph(cx=graph)

gene_list_path = sys.argv[2]

output_folder = "{}/../HeatDiffusion_results/".format(os.path.dirname(gene_list_path))
try:
    os.mkdir(output_folder)
except:
    print "Folder already exists data will be overwritten"

if os.path.isdir(gene_list_path):
    gene_list_path = ["{}{}".format(gene_list_path, gene_list) for gene_list in os.listdir(gene_list_path)]
elif os.path.isfile(gene_list_path):
    gene_list_path = [gene_list_path]
else:
    print "No folder nor file found with that name"
    sys.exit()

try:
    number_of_heated_nodes = int(sys.argv[3])
except:
    print "All heated nodes will be reported"
    number_of_heated_nodes = "all"

host = sys.argv[4]

if host == "local":
    url = 'http://localhost/'
else:
    url = 'http://v3.heat-diffusion.cytoscape.io'

for gene_list in gene_list_path:
    full_workflow(gene_list, output_folder, number_of_heated_nodes, graph, url)
    
