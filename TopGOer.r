#!/usr/bin/env Rscript
# Adrian Garcia Moreno

library(topGO)

human_ensmbl2go <- function(){
  library(biomaRt)
  
  mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice",  dataset="hsapiens_gene_ensembl") # HumanGenome19
  attr <- listAttributes(mart)
  EG2GO <- getBM(mart=mart, attributes=c('ensembl_gene_id','go_id','goslim_goa_accession'))
  EG2GO <- EG2GO[EG2GO$go_id != '',]
  all_ensmbl_genes <- unique(EG2GO[,1])
  
  dictionar <- c()
  for (gene in all_ensmbl_genes){
    gos_aviable <- EG2GO[EG2GO[,1] == gene,]
    unique_goes <- unique(gos_aviable[,2])
    unique_slimgoes <- unique(gos_aviable[,3])
    goes_str <- paste(unique_goes, collapse = ", ")
    slim_goes_str <- paste(unique_slimgoes, collapse = ", ")
    all_gous_str <- paste(goes_str, slim_goes_str, sep = ", ")
    aception <- paste(gene, all_gous_str, sep = "\t")
    dictionar <- rbind(dictionar, aception)
  }
  write.table(dictionar, file = "ENSG2GOmap.map", quote = FALSE, row.names = F, col.names = F)  
  return("ENSG2GOmap.map")
}

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

FullBasicTopGOAnalysis <- function(input, map, mode = c("MF", "CC", "BP"), output_folder, pval_thres){
  
  geneID2GOmap <- readMappings(map)
  writeLines("Map Readed")
  geneID2GOmap_table <- read.table(file = map, header = F, comment.char="", sep="\t", quote="")
  gene_universe <- geneID2GOmap_table[,1]
  writeLines("Gene Universe Acknowledged")  
  classic <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher_Test")
  writeLines("Classic Test ready\nCreating output folder\nReading Input")

  #input <- "./example_data/3col/GeneIDs/"
  
  if(dir.exists(input)){
    output_path <- paste(dirname(input), "/", sep = "")
    input_files <- paste(input, list.files(input), sep = "")
  }else{
    output_path <- paste(dirname(input), "/../", sep = "")
    input_files <- c(input)
  }

  output_folder <- paste(output_path, output_folder, "_GOanalysis/", sep = "")
  system(paste("mkdir ", output_folder, sep = ""))
  
  
  total_files <- length(input_files)
  counter <- 1
  for (file in input_files){
    name <- strsplit(basename(file), "[.]")[[1]][1]
    print(sprintf("Reading %s    %i of %i", basename(file), counter, total_files))
    counter <- counter + 1
    problematic_HP <- 1
    
    linescheck <- readLines(paste(file))
    if(identical(linescheck,character(0))){
      write(paste(name, " none genes found", as.character(problematic_HP), sep = ""), file = "warnings_enrichment.txt", append = TRUE, sep = "\n")
      problematic_HP <- problematic_HP + 1
    next}
    query_score_table <- read.table(text = gsub("\t", " ", readLines(file)), header = F, comment.char="", sep=" ", quote="", fill = T)
    query <- query_score_table[,1]  

    genes_list <- factor(as.integer(gene_universe %in% query))
    names(genes_list) <- gene_universe
    
    if (length(levels(genes_list)) == 1){
      write(paste(name, " genes not mapped", as.character(problematic_HP), sep = ""), file = "warnings_enrichment.txt", append = TRUE, sep = "\n")
      problematic_HP <- problematic_HP + 1
    next}

    if ("MF" %in% mode[[1]]){
    print("Generating Mollecular Function Results")
    MF_GOobject <- new("topGOdata", ontology = "MF", allGenes = genes_list, annot = annFUN.gene2GO, gene2GO = geneID2GOmap)
    MF_resultclassic <- getSigGroups(MF_GOobject, classic)
    MF_tops <- length(score(MF_resultclassic)[score(MF_resultclassic) < pval_thres])
    MF_results_table <- GenTable(MF_GOobject, classic = MF_resultclassic, orderBy = "classic", ranksOf = "classic", topNodes = MF_tops)
    write.table(MF_results_table, file = paste("MF", name, sep = "_"), sep = "\t", quote = FALSE, row.names = F, col.names = T)
    printGraph(MF_GOobject, MF_resultclassic, firstSigNodes = 10, fn.prefix = paste("MFgraph", name, sep = "_"), useInfo = "all", pdfSW = T)
    }

    if ("BP" %in% mode[[1]]){
    print("Generating Biological Process Results")
    BP_GOobject <- new("topGOdata", ontology = "BP", allGenes = genes_list, annot = annFUN.gene2GO, gene2GO = geneID2GOmap)
    BP_resultclassic <- getSigGroups(BP_GOobject, classic)
    BP_tops <- length(score(BP_resultclassic)[score(BP_resultclassic) < pval_thres])
    BP_results_table <- GenTable(BP_GOobject, classic = BP_resultclassic, orderBy = "classic", ranksOf = "classic", topNodes = BP_tops)
    write.table(BP_results_table, file = paste("BP", name, sep = "_"), sep = "\t", quote = FALSE, row.names = F, col.names = T)
    printGraph(BP_GOobject, BP_resultclassic, firstSigNodes = 10, fn.prefix = paste("BPgraph", name, sep = "_"), useInfo = "all", pdfSW = T)
    }
    
    if ("CC" %in% mode[[1]]){
    print("Generating Cellular Component Results")
    CC_GOobject <- new("topGOdata", ontology = "CC", allGenes = genes_list, annot = annFUN.gene2GO, gene2GO = geneID2GOmap)
    CC_resultclassic <- getSigGroups(CC_GOobject, classic)
    CC_tops <- length(score(CC_resultclassic)[score(CC_resultclassic) < pval_thres])
    CC_results_table <- GenTable(CC_GOobject, classic = CC_resultclassic, orderBy = "classic", ranksOf = "classic", topNodes = CC_tops)
    write.table(CC_results_table, file = paste("CC", name, sep = "_"), sep = "\t", quote = FALSE, row.names = F, col.names = T)
    printGraph(CC_GOobject, CC_resultclassic, firstSigNodes = 10, fn.prefix = paste("CCgraph", name, sep = "_"), useInfo = "all", pdfSW = T)
    }
    
  individual_results_folder <- paste(output_folder, name, "/", sep = "")
  system(paste("mkdir ", individual_results_folder, sep = ""))
  
  if (length(mode[[1]]) <= 2){
    if (length(mode[[1]]) == 1){
    system(sprintf(paste("mv %s* ", individual_results_folder, sep=""), as.character(mode[[1]][1])))
    }else{system(sprintf(paste("mv %s* %s* ", individual_results_folder, sep=""), as.character(mode[[1]][1]), as.character(mode[[1]][2])))}
  }else{system(paste("mv MF* BP* CC* ", individual_results_folder, sep=""))}
  
  } # This is the one that closes the for loop and the next is an small checking
  
  if ("warnings_enrichment.txt" %in% list.files(".")){
    print("Some HP files couldn't be processed, take a look at warnings_enrichment.txt in your output folder")
    system(paste("mv warnings_enrichment.txt ", output_folder, sep=""))
    }else {print("Nice! No Warnings!")}
}

orders <- commandArgs(trailingOnly = TRUE)

orders <- c(orders, "NA")                          # Because I have an uneven number of arguments and...
args_matrix <- matrix(orders, ncol = 2, byrow = T) # to create something similar to a dic
                                                   # I create this matrix of 2 columns
# These args are mandatory
input <- args_matrix[1,1]
map <- args_matrix[1,2]
output_folder <- args_matrix[length(args_matrix[,1]),1]

# These args are optional
taxid <- args_matrix[match("-taxid", args_matrix),2]
if (is.na(taxid)){
  taxid <- "9606"
}

mode <- args_matrix[match("-mode", args_matrix),2]
if (is.na(mode)){
  mode <- list(c("MF", "CC", "BP"))
}else{
  mode <- strsplit(mode, ", |,|-| ")
  }

pval_thres <- args_matrix[match("-pval_thres", args_matrix),2]
if (is.na(pval_thres)){
  pval_thres <- 0.05
}else{
  pval_thres <- as.numeric(pval_thres)
}

if (substr(map, nchar(map)-3, nchar(map)) == ".map"){
  map <- map
 }else{
  if (map == "gene2go:download"){
    print("Downloading Gene2Go table")
    map <- gene2go_downloader()
    print("Generating map file 1/2")
    gene2go_tax_extrator(taxid, map)
    print("Generating map file 2/2")
    map <- mapper(taxid)
    print("Map file generated")
  }
  if (map == "gene2go"){
    print("Generating map file 1/2")
    gene2go_tax_extrator(taxid, map)
    print("Generating map file 2/2")
    map <- mapper(taxid)
    print("Map file generated")
  }
}

# Rscript Enrichment/U-TopGOFullBasic.r gene_lists/ 
# gene2go:download/gene2go/.*.map <taxid> 
# <mode(MF-CC-BP)> <pval_thres=0.05> output
FullBasicTopGOAnalysis(input, map, mode, output_folder, pval_thres)

