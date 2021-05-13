## Given a list of data metrics, this script will go through each data metric, match the gene symbol to official gene symbol

library(dplyr)
source("./R/function.R")
options(stringsAsFactors = FALSE)

DATALIST <- "./DATA/fList_example.txt"
pdir = "./"
setwd("./") # working director

geneAliasDict <- list("mouse" = "./REFERENCE/Mus_musculus.gene_info.gz",
                      "human" = "./REFERENCE/Homo_sapiens.gene_info.gz")

taxid <- list("mouse"= "10090", 
              "human"="9606")

## Calculate evidence factor for each data metrics
PATH <- read.table(DATALIST, header=F) %>%
  magrittr::set_colnames(c("Source", "ID","Input","Models","DataType")) %>%
  dplyr::mutate(rID=gsub("-",".",ID)) %>%
  dplyr::select(-ID) %>%
  dplyr::rename("ID"="rID")
  
## read gene used as labels    
inputsleepGene <- read.table("./DATA/inputSleepGeneList.txt", header=FALSE, sep="\t") %>%
  magrittr::set_colnames(c("GeneSymbol", "Tier"))
inputsleepGeneList <- inputsleepGene$GeneSymbol

## Read the list of gene for evaluation
geneListBlank <- read.table("./DATA/GeneList.txt",sep="\t", quote="", header=T)
geneListBlank <- labelData(geneListBlank, sleepGeneList = inputsleepGeneList)

## Create empty data based on the set of all known gene
rawValuetable <- geneListBlank

## Add new data metrics to rawValuetable
rawValuetable <- addRawValue(PATH, rawValuetable, inputsleepGeneList)
write.csv(rawValuetable, "./OUT/fList_rawValuetable.csv", row.names=FALSE)
