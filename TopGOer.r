#!/usr/bin/env Rscript
# Adrian Garcia Moreno

library(topGO)

gene2go_downloader <- function(){
  system("wget ftp://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz")
  system("gunzip gene2go.gz")
  return("gene2go")
}

gene2go_tax_extrator <- function(taxid, map){
  system(paste("grep ^", taxid, " ", map," > ", taxid, "_gene2go.tsv", sep = ""))
  return()
}

mapper <- function(taxid){
  gene2goncbi <- read.table(paste(taxid, "_gene2go.tsv", sep = ""), header = F, comment.char="", sep="\t", quote="")
  gene2goncbi <- gene2goncbi[,c(2,3)]
  gene_universe <- unique(gene2goncbi[,1]) 
  dictionar <- c()
  for (gene in gene_universe){
    gos_aviable <- gene2goncbi[gene2goncbi[,1] == gene,]
    uniquegoes <- unique(gos_aviable[,2])
    goes_str <- paste(uniquegoes, collapse = ", ")
    aception <- paste(gene, goes_str, sep = "\t")
    dictionar <- rbind(dictionar, aception)
  }
  write.table(dictionar, file = paste(taxid, "_geneID2GO.map", sep = ""), quote = FALSE, row.names = F, col.names = F)
  return(paste(taxid, "_geneID2GO.map", sep = ""))
}

input_checker <- function (input){
  if(dir.exists(input)){
    output_path <- paste(dirname(input), "/", sep = "")
    input_files <- paste(input, list.files(input), sep = "")
    return(input_files)
  } else if(file.exists(input)){
    output_path <- paste(dirname(input), "/../", sep = "")
    input_files <- c(input)
    return(input_files)
  } else{writeLines(sprintf(
    "%s Don't exist on your files, please check your input", input))}
}

scores_gene_list_extractor <- function(input_file, gene_col, gene_score, gene_score_col){
  bed_file <- read.table(input_file, header = F, sep = "\t", quote = "")
  genes_list <- bed_file[bed_file[,gene_score_col] > gene_score, gene_col]
  return(genes_list)
}

noscores_gene_list_extractor <- function(input_file, gene_col){
  bed_file <- read.table(input_file, header = F, sep = "\t", quote = "")
  genes_list <- bed_file[, gene_col]
  return(genes_list)
}

results_with_sig_genes <- function(TopGOobject, classic, pval_thres, genes_query, results_filename){
  resultsclassic <- getSigGroups(TopGOobject, classic)
  tops <- length(score(resultsclassic)[score(resultsclassic) < pval_thres])
  results_table <- GenTable(TopGOobject, classic = resultsclassic, orderBy = "classic", ranksOf = "classic", topNodes = tops)
  genes_in_goes <- genesInTerm(TopGOobject, results_table$GO.ID)
  colnames(results_table)[6] <- "Fisher Test"
  intersected <- lapply(genes_in_goes, function(x) intersect(x, genes_query))
  genes_of_goes <- unlist(lapply(intersected, function(x) paste(x, collapse = ":")))
  results_table <- cbind(results_table, genes_of_goes)
  colnames(results_table)[7] <- "Significant Genes"
  write.table(results_table, file = results_filename, sep = "\t", quote = FALSE, row.names = F, col.names = T)
}

FullBasicTopGOAnalysis <- function(genes_list, name, map, mode = c("MF", "CC", "BP"), output_folder, classic, pval_thres){
  genes_query <- names(genes_list)[genes_list == 1]
  for (ont in mode){
    print(sprintf("Generating %s Results", ont))
    results_filename <- paste(output_folder, "/", ont, "_", name, ".txt", sep = "")
    TopGOobject <- new("topGOdata", ontology = ont, allGenes = genes_list, annot = annFUN.gene2GO, gene2GO = map)
    results_with_sig_genes(TopGOobject, classic, pval_thres, genes_query, results_filename)
  }
}

orders <- commandArgs(trailingOnly = TRUE)

#orders <- c(orders, "NA")                          # Because I have an uneven number of arguments and...
args_matrix <- matrix(orders, ncol = 2, byrow = T) # to create something similar to a dictionary
                                                   # I create this matrix of 2 columns
# These args are mandatory
input_file <- args_matrix[1,1]
map_file <- args_matrix[1,2]
gene_col <- as.integer(args_matrix[2,1])
output_folder <- args_matrix[2,2]
dir.create(path = output_folder)

# These args are optional
score <- args_matrix[match("-score", args_matrix),2]
if (!is.na(score)){
  scores <- unlist(strsplit(score, ",|-"))
  gene_score_col <- as.integer(scores[1])
  gene_score <- as.numeric(scores[2])
}

taxid <- args_matrix[match("-taxid", args_matrix),2]
if (is.na(taxid)){
  taxid <- "9606"
}

mode <- args_matrix[match("-mode", args_matrix),2]
if (is.na(mode)){
  mode <- c("MF", "CC", "BP")
}else{
  mode <- unlist(strsplit(mode, ",|-"))
  }

pval_thres <- args_matrix[match("-pval_thres", args_matrix),2]
if (is.na(pval_thres)){
  pval_thres <- 0.05
}else{
  pval_thres <- as.numeric(pval_thres)
}

if (substr(map_file, nchar(map_file)-3, nchar(map_file)) == ".map"){
}else{
  if (map_file == "gene2go:download"){
    print("Downloading Gene2Go table")
    map_file <- gene2go_downloader()
    print("Generating map file 1/2")
    gene2go_tax_extrator(taxid, map_file)
    print("Generating map file 2/2")
    map_file <- mapper(taxid)
    print("Map file generated")
  }
  if (map_file == "gene2go"){
    print("Generating map file 1/2")
    gene2go_tax_extrator(taxid, map_file)
    print("Generating map file 2/2")
    map_file <- mapper(taxid)
    print("Map file generated")
  }}

# Rscript Enrichment/U-TopGOFullBasic.r gene_lists/ 
# gene2go:download/gene2go/.*.map <taxid> 
# <mode(MF-CC-BP)> <pval_thres=0.05> output

# USER DEFINED CONSTANTS
# EXAMPLE
# input_folder <- "annotated/" || input_file <- "Desktop/TeoreticalHPGenes/HP:0000089.txt"
# gene_col <- 10
# gene_score_col <- 5
# gene_score <- 2
# pval_thres <- 0.01
#map_file <- "Desktop/Rarebiosis/annotation_files/9606_geneID2GO.map"
# output_folder <- "./"
# mode <- c("MF", "CC", "BP")

errors <- ""
input_files <- input_checker(input_file)
total_files <- length(input_files)
counter <- 1

map <- readMappings(map_file)
gene_universe <- read.table(file = map_file, header = F, comment.char="", sep="\t", quote="")[,1]

classic <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher_Test")

for (input_file in input_files){
  
  name <- unlist(strsplit(basename(input_file), split = "[_|.]"))
  name <- tail(name, n = 2)[1]
  print(sprintf("Reading %s    %i of %i", name, counter, total_files))
  counter <- counter + 1
  
  # 1º Possible Error == Empty File
  if (file.info(input_file)$size == 0){
    errors <- paste(errors, name, "\tNO GENES\n", sep = "")
    next}
  
  if (is.na(score)){
    genes_list <- noscores_gene_list_extractor(input_file, gene_col)    
  } else{
    genes_list <- scores_gene_list_extractor(input_file, gene_col, gene_score, gene_score_col)
  }
  
  genes_list <- factor(as.integer(gene_universe %in% genes_list))
  names(genes_list) <- gene_universe
  
  # 2º Possible Error == No gene-score above cut-off
  if (length(genes_list) == 0){
    errors <- paste(errors, name, "\tNO GENE-SCORE ABOVE CUT-OFF\n", sep = "")
    next}
  
  # 3º Possible Error == No mapping genes
  if(!any(genes_list %in% gene_universe)){
    errors <- paste(errors, name, "\tNO GENE MAP TO ONTOLOGY\n", sep = "")
    next}  
  
  FullBasicTopGOAnalysis(genes_list, name, map, mode, output_folder, classic, pval_thres)
}

if (errors == ""){
  print("Nice! No Warnings!")
  }else{
  error_file <- paste(output_folder, "ERRORS.txt", sep = "")
  write(errors, error_file)
  print("Some HP files couldn't be processed, take a look at ERRORS.txt in your output folder")
}
