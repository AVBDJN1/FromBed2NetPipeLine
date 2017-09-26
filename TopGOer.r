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

gene2go_tax_extrator <- function(taxid, gene2go_map){
  system(paste("grep ^", taxid, " ", gene2go_map," > ", taxid, "_gene2go.tsv", sep = ""))
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

FullBasicTopGOAnalysis <- function(input, map, mode = c("MF", "CC", "BP"), output){
  if (dir.exists(input)){
    files <- list.files(path = input)
  }else{files <- c(input)}
  
  system(paste("mkdir ", input, "../", output, "_GOanalysis/", sep = ""))
  system(paste("mkdir ", input, "../", output, "_GOanalysis/results/", sep = ""))
  #system(paste("mkdir ", input, "../", output, "_GOanalysis/results/", sep = ""))

  geneID2GOmap <- readMappings(map)
  print("Map Readed")
  geneID2GOmap_table <- read.table(file = map, header = F, comment.char="", sep="\t", quote="")
  gene_universe <- geneID2GOmap_table[,1]
  print("Gene Universe Acknowledged")  
  classic <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher_Test")
  print("Classic Test ready")
  total_files <- length(files)
  counter <- 1
  for (file in files){
    name <- substr(file, nchar(file)-13, nchar(file)-4)
    print(sprintf("Reading %s    %i of %i", file, counter, total_files))
    counter <- counter + 1
    problematic_HP <- 1
    
    if (dir.exists(input)){
      files <- list.files(path = input)
      linescheck <- readLines(paste(input,"/",file, sep = ""))
      if(identical(linescheck,character(0))){
        write(paste(name, " none genes found", as.character(problematic_HP), sep = ""), file = "warnings_enrichment.txt", append = TRUE, sep = "\n")
        problematic_HP <- problematic_HP + 1
        next}
      query_score_table <- read.table(file = paste(input,"/",file, sep = ""), header = F, comment.char="", sep=" ", quote="")
      query <- query_score_table[,1]  
    }else{files <- c(input)
    query_score_table <- read.table(file = input, header = F, comment.char="", sep=" ", quote="")
    query <- query_score_table[,1]
    linescheck <- readLines(input)
    if(identical(linescheck,character(0))){
      write(paste(name, " none genes found", as.character(problematic_HP), sep = ""), file = "warnings_enrichment.txt", append = TRUE, sep = "\n")
      problematic_HP <- problematic_HP + 1
      next}
    }

    genes_list <- factor(as.integer(gene_universe %in% query))
    names(genes_list) <- gene_universe
    skkiped <- ""
    if (length(levels(genes_list)) == 1){
      write(paste(name, " genes not mapped", as.character(problematic_HP), sep = ""), file = "warnings_enrichment.txt", append = TRUE, sep = "\n")
      problematic_HP <- problematic_HP + 1
      next
    }

    if ("MF" %in% mode[[1]]){
    print("Generating Mollecular Function Results")
    MF_GOobject <- new("topGOdata", ontology = "MF", allGenes = genes_list, annot = annFUN.gene2GO, gene2GO = geneID2GOmap)
    MF_resultclassic <- getSigGroups(MF_GOobject, classic)
    MF_tops <- length(score(MF_resultclassic)[score(MF_resultclassic) < 0.05])
    MF_results_table <- GenTable(MF_GOobject, classic = MF_resultclassic, orderBy = "classic", ranksOf = "classic", topNodes = MF_tops)
    write.table(MF_results_table, file = paste("MF", name, sep = "_"), sep = "\t", quote = FALSE, row.names = F, col.names = T)
    printGraph(MF_GOobject, MF_resultclassic, firstSigNodes = 10, fn.prefix = paste("MFgraph", name, sep = "_"), useInfo = "all", pdfSW = T)
    }

    if ("BP" %in% mode[[1]]){
    print("Generating Biological Process Results")
    BP_GOobject <- new("topGOdata", ontology = "BP", allGenes = genes_list, annot = annFUN.gene2GO, gene2GO = geneID2GOmap)
    BP_resultclassic <- getSigGroups(BP_GOobject, classic)
    BP_tops <- length(score(BP_resultclassic)[score(BP_resultclassic) < 0.05])
    BP_results_table <- GenTable(BP_GOobject, classic = BP_resultclassic, orderBy = "classic", ranksOf = "classic", topNodes = BP_tops)
    write.table(BP_results_table, file = paste("BP", name, sep = "_"), sep = "\t", quote = FALSE, row.names = F, col.names = T)
    printGraph(BP_GOobject, BP_resultclassic, firstSigNodes = 10, fn.prefix = paste("BPgraph", name, sep = "_"), useInfo = "all", pdfSW = T)
    }
    
    if ("CC" %in% mode[[1]]){
    print("Generating Cellular Component Results")
    CC_GOobject <- new("topGOdata", ontology = "CC", allGenes = genes_list, annot = annFUN.gene2GO, gene2GO = geneID2GOmap)
    CC_resultclassic <- getSigGroups(CC_GOobject, classic)
    CC_tops <- length(score(CC_resultclassic)[score(CC_resultclassic) < 0.05])
    CC_results_table <- GenTable(CC_GOobject, classic = CC_resultclassic, orderBy = "classic", ranksOf = "classic", topNodes = CC_tops)
    write.table(CC_results_table, file = paste("CC", name, sep = "_"), sep = "\t", quote = FALSE, row.names = F, col.names = T)
    printGraph(CC_GOobject, CC_resultclassic, firstSigNodes = 10, fn.prefix = paste("CCgraph", name, sep = "_"), useInfo = "all", pdfSW = T)
    }
    
  system(paste("mkdir ", input, "../", output, "_GOanalysis/results/", name, "/", sep = ""))
  if (length(mode[[1]]) <= 2){
    if (length(mode[[1]]) == 1){
      system(sprintf(paste("mv %s* ", input, "../", output, "_GOanalysis/results/", name, "/", sep=""), as.character(mode[[1]][1])))
    }else{system(sprintf(paste("mv %s* %s* ", input, "../", output, "_GOanalysis/results/", name, "/", sep=""), as.character(mode[[1]][1]), as.character(mode[[1]][2])))}
  }else{system(paste("mv MF* BP* CC* ", input, "../", output, "_GOanalysis/results/", name, "/", sep=""))}
  }
  if ("warnings_enrichment.txt" %in% list.files(".")){
    print("Some HP files couldn't be processed, take a look at warnings_enrichment.txt in your output folder")
    system(paste("mv warnings_enrichment.txt ", input, "../", output, "_GOanalysis/", sep=""))
    }else {print("Nice! No Warnings!")}
}

orders <- commandArgs(trailingOnly = TRUE)

input <- orders[1]
gene2go_map <- orders[2]

if (length(orders) == 5){
  taxid <- orders[3]
  mode <- orders[4]
  output <- orders[5]
  print("2")
} 
if (length(orders) == 4){
  if (is.na(as.integer(orders[3]))){
    mode <- strsplit(orders[3], ", |,|-| ")
    output <- orders[4]
    print("3")
  }else{
    taxid <- orders[3]
    mode <- list(c("MF", "CC", "BP"))
    output <- orders[4]
    print("4")
  }  
}
if (length(orders) == 3){
  output <- orders[3]
  mode <- list(c("MF", "CC", "BP"))
  cat("5")
}
if (substr(gene2go_map, nchar(gene2go_map)-3, nchar(gene2go_map)) == ".map"){
  map <- gene2go_map
}
if (gene2go_map == "gene2go:download"){
  print("Downloading Gene2Go table")
  gene2go_map <- gene2go_downloader()
  print("Generating map file 1/2")
  gene2go_tax_extrator(taxid, gene2go_map)
  print("Generating map file 2/2")
  map <- mapper(taxid)
}
if (gene2go_map == "gene2go"){
  print("Generating map file 1/2")
  gene2go_tax_extrator(taxid, gene2go_map)
  print("Generating map file 2/2")
  map <- mapper(taxid)
}

# Rscript Enrichment/U-TopGOFullBasic.r ncbi20_run/genes_lists/ analisis_data/Annotation_sources/TopGo_stuff/9606_geneID2GO.map enrich
# Rscript Enrichment/U-TopGOFullBasic.r ncbi20_run/genes_lists/ gene2go:download 9606 enrich

FullBasicTopGOAnalysis(input, map, mode, output)
