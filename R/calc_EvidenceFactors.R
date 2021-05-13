## Script to calculate evidence factors given a list of data metrics

library(ggplot2)
library(ggridges)
library(dplyr)
library(tidyr)
library(tidyverse)
library(magrittr)
library(ggrepel)
library(ggpubr)
library(cowplot)
source("./R/function.R")
options(stringsAsFactors = FALSE)

pdir = "./"
setwd("./") # working director
curPath <- "./OUT/" # directory to write output file
curDate <- '20210503' ## working date, use to create unique identity for output

geneAliasDict <- list("mouse" = "./REFERENCE/Mus_musculus.gene_info.gz",
                      "drosophila" = "./REFERENCE/Drosophila_melanogaster.gene_info.gz",
                      "human" = "./REFERENCE/Homo_sapiens.gene_info.gz")

taxid <- list("mouse"= "10090", 
              "drosophila" ="7227",
              "human"="9606")
    
## read gene used as labels    
inputsleepGene <- read.table("./DATA/inputSleepGeneList.txt", header=FALSE, sep="\t") %>%
  magrittr::set_colnames(c("GeneSymbol", "Tier"))
inputsleepGeneList <- inputsleepGene$GeneSymbol
minInputSleepGeneCuttoff <- length(inputsleepGeneList)*.25

## Read the list of gene for evaluation
geneListBlank <- read.table("./DATA/GeneList.txt",sep="\t", quote="", header=T)
geneListBlank <- labelData(geneListBlank, sleepGeneList = inputsleepGeneList)

## Create empty data based on the set of all known gene
rawValuetable <- geneListBlank
EFtable <- geneListBlank

## Calculate evidence factor for each data metrics
PATH <- read.table("combine_fList_updateGeneList.txt",header=F) %>%
  magrittr::set_colnames(c("Source", "ID","Input","Models","DataType")) %>%
  dplyr::mutate(rID=gsub("-",".",ID)) %>%
  dplyr::select(-ID) %>%
  dplyr::rename("ID"="rID") %>%
  #dplyr::mutate(maxEF=NA, upper99.EF=NA)
  dplyr::mutate(EFscore=NA)

## read existing table with data metrics raw value
rawValuetable <- read.csv(paste0("./fList_rawValuetable.csv"), header=T) %>%
  dplyr::select(-label)
rawValuetable <- labelData(rawValuetable, sleepGeneList=inputsleepGeneList)

## Add new data metrics to rawValuetable
rawValuetable <- addRawValue(PATH, rawValuetable, inputsleepGeneList)
write.csv(rawValuetable, "./fList_rawValuetable.csv", row.names=FALSE)

pdf(paste0(curPath,"EFbin_",curDate,".pdf"),8,3)
p <- list()
for(i in 3:ncol(rawValuetable)){
  ID <- colnames(rawValuetable)[i]
  if(i%%100==0){print(c("processing column: ", i, ID))}
  
  dataType <- PATH[PATH$ID==ID,"DataType"]
  dat <- rawValuetable[,c(1,2,i)] %>%
    drop_na() %>%
    magrittr::set_colnames(c("GeneSymbol", "label", "value"))
  
  if(nrow(subset(dat, dat$label=="SRG"))<minInputSleepGeneCuttoff){
    print(paste(c("skip as less than 25% sleeep genes found in the data:", ID)))
    next
  }
  
  ef <- calEF_bin(dat,dataType)
  if(nrow(ef) < 1){
    print(paste(c("skip as equal data value:", ID)))
    next
  }
  
  p[[1]] <- drawDensityPlot(dat, nCol=3)
  p[[2]] <- drawEFPlot(ef, name = ID)
  do.call("grid.arrange", c(p, ncol=2))
  
  datEF <- matchEF_bin(dat,ef)

  EFtable <- EFtable %>%
    left_join(datEF, by=c("GeneSymbol")) %>%
    dplyr::rename(!!ID := evidenceFactors)

  datEF <- labelData(datEF, sleepGeneList=inputsleepGeneList)

  if(length(table(ef$evidenceFactors))>1){
    a <- subset(datEF, datEF$label=="SRG")
    #PATH[PATH$ID==ID, "EFscore"] = mean(a[a$evidenceFactors>=quantile(a$evidenceFactors, .9),"evidenceFactors"])
    PATH[PATH$ID==ID, "EFscore"] = max(a$evidenceFactors)
  }
}
dev.off()

write.table(PATH, paste0(curPath, "fList_PATH_",curDate,".txt"),sep="\t", row.names=FALSE)
write.csv(EFtable, paste0(curPath, "fList_EFtable_",curDate,".csv"), row.names=FALSE)

## select feature with high evidence
zz <- data.frame(ID=character(),
                 maxEF=numeric(),
                 stringsAsFactors=FALSE)
for(i in 3:ncol(EFtable)){
  ID <- colnames(EFtable)[i]
  maxEF <- max(EFtable[,i],na.rm=TRUE)
  zz <- rbind(zz, data.frame(ID, maxEF))
}

yy <- zz %>%
  dplyr::arrange(desc(maxEF)) %>%
  dplyr::filter(maxEF>=3)
head(yy, 20)

p1 <- ggplot(zz, aes(x=maxEF)) + 
  geom_histogram(bins=500, fill="darkblue", alpha=0.8) +
  #geom_density() +
  xlab("maxEF") +
  geom_vline(xintercept = 3, linetype="dashed", color="#777777") + 
  scale_x_continuous(trans='log2') +
  theme_bw()
ggsave(p1, filename = paste0(curPath,"EF_distribution",curDate,".pdf"),width = 2.5, height = 2)

cleanRawValuetable <- rawValuetable %>%
  dplyr::select(c("GeneSymbol","label",c(names(yy)))) %>%
  dplyr::mutate(naCount = apply(.[-c(1,2)],1,function(x){sum(is.na(x))})) %>%
  dplyr::filter(naCount<length(yy)* .5) %>%
  dplyr::select(-naCount)

cleanEFtable <- EFtable %>%
  dplyr::select(c("GeneSymbol","label",names(yy))) %>%
  dplyr::mutate(naCount = apply(.[-c(1,2)],1,function(x){sum(is.na(x))})) %>%
  dplyr::filter(naCount<length(yy)* .5) %>%
  dplyr::select(-naCount)

## calculate pairwise correlation for selected samples wiht high evidence factors
corTable <- as.data.frame(cor(cleanRawValuetable[,3:ncol(cleanRawValuetable)], use="pairwise.complete.obs"))

corTable$ID <- 1:(length(yy))
dataName <- colnames(cleanRawValuetable)[-c(1,2)]
vv <- melt(corTable, id.vars = "ID") %>%
  dplyr::mutate(variable=as.numeric(variable)) %>%
  dplyr::filter(ID!=variable, ID<variable, value>=0.8) %>%
  dplyr::mutate(sp1=dataName[.$ID], sp2=dataName[.$variable])

cleanRawValuetable <- rawValuetable %>%
  dplyr::select(c("GeneSymbol","label",c(names(yy)))) %>%
  dplyr::mutate(naCount = apply(.[-c(1,2)],1,function(x){sum(is.na(x))})) %>%
  dplyr::filter(naCount<length(yy)* .5) %>%
  dplyr::select(-vv$sp2) %>%
  dplyr::select(-naCount)

cleanEFtable <- EFtable %>%
  dplyr::select(c("GeneSymbol","label",names(yy))) %>%
  #dplyr::select(-(vv$sp2[vv$sp2%in%names(yy)])) %>%
  dplyr::mutate(naCount = apply(.[-c(1,2)],1,function(x){sum(is.na(x))})) %>%
  dplyr::filter(naCount<length(yy)* .5) %>%
  dplyr::select(-vv$sp2) %>%
  dplyr::select(-naCount)

write.csv(cleanRawValuetable, paste0("fList_cleanRawValue_",curDate,".csv"), row.names=FALSE)
write.csv(cleanEFtable, paste0("fList_cleanEFtable_",curDate,".csv"), row.names=FALSE)
