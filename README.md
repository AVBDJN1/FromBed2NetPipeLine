# Project Features-Regions-Explorer - by EidrianGM Adrian Garcia Moreno
(all the code given as an example can be run by the folder example_data)
Imagine that we have a bed file where a genomic regions is associated to a certain feature or whatever. This asociation can have an score.
We might like to do the analysis for each feature, or even, for each feautre whose genomic region has a score above a certain threshold.
Here I give you a group of scripts that it might interest you to solve this issue.
I recommend you to use directly the **FULL ANNOTATION MODE** and later if everything goes well.. perform the enrichment.

```
#usr/bin/bash

python ../FromBed2Net.py -HP_selection initial_file_4testing.bed 5-3-4
```

This first mode `-HP_selection` just retrieve and divide the genomic regions by its feature. The input file is `initial_file_4testing.bed`
5-3-4 means that the threshold is in the fifth column and the score should be above or equal to three, and the criteria is in the fourth column.
we should obtain a folder with the division and a .txt with all the features together above the threshold.

## ANNOTATION

In case that you already have a bed file o a directory with bed files you can start annotating using `-annotation` mode later you must give the input and download the annotation or give one, lastly you must set the name of the output folder where results will be stored.

### Download annotation mode:
```
python ../FromBed2Net.py -annotation initial_file_4testing_3-up_divided/ -download:ncbi_hg19 ncbi_hg19_run
```

Other available annotation sources to download: (any further suggestion will be thanked)

* `python ../FromBed2Net.py -annotation initial_file_4testing_3-up_divided/ -download:ncbi_hg20 ncbi_hg20_run`
* `python ../FromBed2Net.py -annotation initial_file_4testing_3-up_divided/ -download:ncbi_hg20 ncbi_hg20_run`
* `python ../FromBed2Net.py -annotation initial_file_4testing_3-up_divided/ -download:gcode_hg20 gcode_hg20_run`

As you might notice, you can download the annotation of the previous or nowadays assembly for the annotation of NCBI(EntrezID) or GENECODE(Ensemble)

A way to run this also is giving it any gff3 file with a format compatible with our initial file (should have the same chromosome annotation)
```
python ../FromBed2Net.py -annotation initial_file_4testing_3-up_divided/ ncbi_hg19.gff3 ncbi_hg19_run
```

works also with a given query file downloading a gff3 or providing it
```
python ../FromBed2Net.py -annotation initial_file_4testing_3-up_divided/Feature_N1.bed ncbi_hg19.gff3 ncbi_hg19_run
```

## TSVER
Middle transformation for an easier annotation extraction handling

works with original output folder or any folder with a folder inside called ".*annotated/"

```python ../FromBed2Net.py -tsver ncbi_hg19_run/```

works with first forced output folder or any folder
```
python ../FromBed2Net.py -tsver ncbi_hg19_run/annotated/
```
works with any annotated file
```
python ../FromBed2Net.py -tsver ncbi_hg19_run/annotated/annoted_Feature_N1.bed
```

## GENE EXTRACTION
takes the input like TSVER but this script extracts the genes and the score related to a feature

```
python ../FromBed2Net.py -extract_genes ncbi_hg19_run/
python ../FromBed2Net.py -extract_genes ncbi_hg19_run/tsvs/
python ../FromBed2Net.py -extract_genes ncbi_hg19_run/tsvs/annoted_Feature_N1.tsv
```

## **FULL ANNOTATION MODE**

All the previous steps can be done with a single command:
```
python ../FromBed2Net.py -full initial_file_4testing.bed 5-3-4 ncbi_hg19.gff3 example/
```
## ENRICHMENT

This following script that I am going to explain is based on TopGO, and performs an enrichment analysis.
You can perform this by having a mapping file gene-to-goID, and genes lists as query. In case that you do not have a mapping file, the script can download for you geneID2GO map of ncbi (<ftp://ftp.ncbi.nih.gov/gene/DATA/>).
This file contains the mapping of several taxons, by default it is going to retrieve the mappings of the tax 9606, homo sapiens, however you can retrieve the tax id that you want. We run this analisys using Rscript (what is between <> is optional):

```
Rscript ../TopGOer.r annotation_results/genes_lists/ gene2go:download/gene2go/.*.map <-taxid="9060" taxid> <-mode="MF-CC-BP" "BP-CC"> <-pval_thres=0.05 "0.03"> output
```

Now the optional arguments must have a previus flag to be identified, as I hope that you understands, each of one has a default value (human tax, all GO categories, pval_threshold 0.05) which will be their value in case that they are not explicitly given. 

### Download mode
```
Rscript ../TopGOer.r annotation_results/genes_lists/ gene2go:download <-taxid="9060" taxid> <-mode="MF-CC-BP" "BP-CC"> <-pval_thres=0.05 "0.03"> output
```
### Extract a new tax of your gene2go file previously downloaded and only two categories
```
Rscript ../TopGOer.r annotation_results/genes_lists/ gene2go -taxid 9596 -mode CC-BP -pval_thres 0.03 enrich
```
### Give a map yourself (it would make no sense determine the taxid here) and only two categories and 0.04 threshold
```
Rscript ../TopGOer.r annotation_results/genes_lists/ 9606_geneID2GO.map -mode CC-MF -pval_thres 0.04 enrich
```

## HOTNET2

>(Working on it)





