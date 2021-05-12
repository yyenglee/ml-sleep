## Nov 28, 2020
## Version 2 - bin the metric, use within or out of the metric for evidence factors

pdir = "~/workDIR/sleep/"
setwd("~/workDIR/sleep")

library(ggplot2)
library(ggridges)
library(dplyr)
library(tidyr)
library(tidyverse)
library(magrittr)
library(ggrepel)
library(ggpubr)
library(cowplot)
source(paste0(pdir,"function.R"))
source(paste0(pdir,"configure.R"))
options(stringsAsFactors = FALSE)

inputsleepGene <- read.table("inputSleepGeneList.txt", header=FALSE, sep="\t") %>%
  magrittr::set_colnames(c("GeneSymbol"))
inputsleepGeneList <- inputsleepGene$GeneSymbol
minInputSleepGeneCuttoff <- length(inputsleepGeneList)*.25

## Read the list of gene for evaluation
geneListBlank <- read.table(paste0(pdir,"GeneList.txt"),sep="\t", quote="", header=T)
geneListBlank <- labelData(geneListBlank, sleepGeneList = inputsleepGeneList)

## Create empty data based on the set of all known gene
rawValuetable <- geneListBlank
EFtable <- geneListBlank

## Calculate evidence factor for each data metrics
curDate <- '20210328' ## working date, use to create unique identity for output
curPath <- "./"

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
newDate <- '20210126'
rawValuetable <- addRawValue(PATH, rawValuetable, inputsleepGeneList)
write.csv(rawValuetable, paste0("./fList_rawValuetable_",newDate,".csv"), row.names=FALSE)

pdf(paste0(curPath,"EFbin_20210328.pdf"),8,3)
p <- list()
for(i in 3:ncol(rawValuetable)){
  ID <- colnames(rawValuetable)[i]
  if(i%%100==0){print(c("processing column: ", i, ID))}
  
  dataType <- PATH[PATH$ID==ID,"DataType"]
  dat <- rawValuetable[,c(1,2,i)] %>%
    drop_na() %>%
    magrittr::set_colnames(c("GeneSymbol", "label", "value"))
  
  # if(nrow(subset(dat, dat$label=="SRG"))<minInputSleepGeneCuttoff){
  #   print(paste(c("skip as less than 25% sleeep genes found in the data:", ID)))
  #   next
  # }
  
  ef <- calEF_bin(dat,dataType)
  if(nrow(ef) <= 1){
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
    PATH[PATH$ID==ID, "EFscore"] = mean(a[a$evidenceFactors>=quantile(a$evidenceFactors, .9),"evidenceFactors"])
  }
}
dev.off()

write.table(PATH, paste0(curPath, "fList_PATH_",curDate,".txt"),sep="\t", row.names=FALSE)
write.csv(EFtable, paste0(curPath, "fList_EFtable_",curDate,".csv"), row.names=FALSE)

## select feature with high evidence
zz <- apply(EFtable[,3:ncol(EFtable)],2,function(x){sort(na.omit(x),decreasing = TRUE)[1]})
#zz <- apply(a[,3:ncol(a)], 2, function(x){z=x[complete.cases(x)]; mean(z[z>=quantile(z, .9)])})
yy <- zz[zz>=5]
head(yy[order(yy, decreasing = TRUE)], 20)
#write.table(as.data.frame(cbind(names(yy), as.numeric(yy))), file="EFg3.txt", row.names=FALSE)

za <-data.frame(zz) %>%
  magrittr::set_colnames(c("evidence_factors"))

p1 <- ggplot(za, aes(x=evidence_factors)) + 
  geom_histogram(bins=500, fill="darkblue", alpha=0.8) +
  #geom_density() +
  xlab("maxEF") +
  geom_vline(xintercept = 3, linetype="dashed", color="#777777") + 
  scale_x_continuous(trans='log2') +
  theme_bw()
p1
ggsave(p1, filename = "./bin_newFullGeneList/EF_distribution.pdf",width = 2.5, height = 2)

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

## with GWAS in predictions
cleanRawValuetable <- read.csv("./bin_newFullGeneList/fList_cleanRawValue_withJI_20201229.csv", header=T)
GWAS <- read.csv("./GWAS_mapped_gene.csv", header=T)
colnames(GWAS)[2:5] <- paste(colnames(GWAS),"GWAS",sep="_")[2:5]
aa <- left_join(geneListBlank, GWAS, by="GeneSymbol") %>%
  dplyr::mutate_if(is.numeric, funs(replace_na(., 0)))
cleanRawValuetable <- cleanRawValuetable %>% 
  left_join(aa, by="GeneSymbol")
write.csv(cleanRawValuetable, paste0("fList_cleanRawValuewithGWAS_",curDate,".csv"), row.names=FALSE)

## process machine learning prediction
## summarize gene prediction
curPath <- "./allPC_w92.6/"
prefix <- "n100_92.6_20210411_allPC"
curDate <- '20210411'
ratio <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
combineGeneScore <- data.frame(GeneSymbol=character(),
                               label=character(), 
                               RFsum=numeric(),
                               RFsumF=numeric(),
                               ratio=numeric(),
                               stringsAsFactors=FALSE)

for(r in ratio){
  print(r)
  traintesttable <- read.csv(paste0(curPath, "traintesttable_", prefix, "_", r,".csv"),header=T)
  RFtable <- read.csv(paste0(curPath, "RFtable_", prefix, "_", r,".csv"),header=T)
  
  geneScoreTable <- data.frame(RFtable[,c("GeneSymbol","label")]) %>%
    dplyr::mutate(RFsum=0)
  
  for(i in 1:100){
    test_gene <- traintesttable[,c(1:3,i+3)]
    RFcount <- RFtable %>%
      dplyr::select(c(1,2,3,i+3)) %>%
      inner_join(test_gene, by="GeneSymbol")
    #geneScoreTable$RFsum <- geneScoreTable$RFsum + apply(RFcount, 1, function(x){ifelse(x[7]==1,ifelse(x[4]==1,1,0),0)})
    geneScoreTable$RFsum <- geneScoreTable$RFsum + apply(RFcount, 1, function(x){as.numeric(ifelse(x[7]==1,ifelse(x[4]>=0.1, x[4], 0),0))})
  }
  
  geneFreq <- apply(traintesttable[,-c(1:3)],1,sum)
  geneScoreTable <- geneScoreTable %>%
    dplyr::mutate(RFsumF=RFsum/geneFreq) %>%
    dplyr::select(GeneSymbol, label, RFsum, RFsumF) %>%
    dplyr::mutate(ratio=r)
  
  combineGeneScore <- bind_rows(combineGeneScore, geneScoreTable)
}

#class1cutoff <- as.numeric(quantile(subset(combineGeneScore, combineGeneScore$GeneSymbol%in%summarizeGeneScore$GeneSymbol)$RFsumF, c(0.925)))
class1cutoff <- 0.45

summarizeGeneScore <- combineGeneScore %>%
  dplyr::filter(RFsumF>0.1) %>%
  #dplyr::filter(ratio %in% c(0.4,0.5,0.6,0.7,0.8)) %>%
  #dplyr::mutate(RFsumF=ifelse(RFsumF>0.1, RFsumF, 0)) %>%
  group_by(GeneSymbol, label) %>%
  #dplyr::summarise(RFscore=sum(RFsumF*ratio)) %>%
  dplyr::summarise(rawScore=mean(RFsumF), RFscore=max(ratio*10)+mean(RFsumF)) %>%
  ungroup() %>%
  dplyr::arrange(-RFscore) %>%
  dplyr::mutate(RFrank=as.numeric(rownames(.))) %>%
  #dplyr::mutate(tier=round(RFscore,0))
  #dplyr::mutate(tier=ifelse((RFscore-floor(RFscore))>0.423&(RFscore-floor(RFscore))<=0.5, floor(RFscore)+0.5, round(RFscore,0))) %>%
  #dplyr::mutate(class=10-ifelse(rawScore>=0.5, round(RFscore), ifelse(rawScore>=0.419, floor(RFscore)+0.5, round(RFscore,0))))
  #dplyr::mutate(class=10-ifelse(rawScore>=0.5, round(RFscore), ifelse(rawScore>=class1cutoff, floor(RFscore)+0.5, round(RFscore,0))))
  dplyr::mutate(class=10-ifelse(rawScore>=class1cutoff, ceiling(RFscore),floor(RFscore)))
table(summarizeGeneScore$label, summarizeGeneScore$class)
write.table(summarizeGeneScore,paste0(curPath, "RFsummarizeOutput",prefix,".txt"), sep="\t", row.names=FALSE)

# summarizeGeneScore <- combineGeneScore %>%
#   dplyr::filter(ratio %in% c(0.2,0.3,0.4,0.5,0.6,0.7,0.8)) %>%
#   dplyr::mutate(RFsumF=ifelse(RFsumF>=0.2, RFsumF, 0)) %>%
#   group_by(GeneSymbol, label) %>%
#   dplyr::summarise(RFscore=mean(RFsumF), class=9-max(ratio[RFsumF>0])*10) %>%
#   ungroup() %>%
#   dplyr::arrange(-RFscore) %>%
#   dplyr::mutate(RFrank=as.numeric(rownames(.))) %>%
#   dplyr::filter(RFscore>0)
# table(summarizeGeneScore$label, summarizeGeneScore$class)

geneScore <- combineGeneScore %>%
  dplyr::select(GeneSymbol, label, RFsumF, ratio) %>%
  dplyr::mutate(RFsumF=ifelse(RFsumF>=0.1, RFsumF, NA)) %>%
  #dplyr::filter(RFsumF>0.1) %>%
  drop_na() %>%
  dplyr::filter(GeneSymbol %in% summarizeGeneScore$GeneSymbol) %>%
  left_join(inputsleepGene, by="GeneSymbol") %>%
  dplyr::mutate(sublabel=paste(label,Tier, sep=".")) %>%
  dplyr::mutate(GeneSymbol=factor(GeneSymbol, levels=summarizeGeneScore$GeneSymbol),
                ratio=factor(ratio, levels=unique(ratio)),
                sublabel=factor(sublabel, levels=c("nonSRG.NA", "SRG.3", "SRG.2", "SRG.1")))
  
tt <- table(summarizeGeneScore$label, summarizeGeneScore$class)
classsize <- cumsum(tt[1,]+tt[2,])
p2 <- ggplot(geneScore, mapping=aes(x=RFsumF, y=factor(GeneSymbol))) +
  geom_point(color="grey", size=0.25, alpha=0.8) +
  #geom_point(geneScore, mapping=aes(x=GeneSymbol, y=RFsumF, color=ratio), size=0.25, alpha=0.8) +
  #geom_point(subset(geneScore, geneScore$label=="SRG"), mapping=aes(x=GeneSymbol, y=RFsumF), color="black", size=0.1, alpha=0.5) +
  stat_summary(fun="mean", geom="point", aes(color=sublabel), alpha=0.9, size=0.8) +
  stat_summary(subset(geneScore,geneScore$label=="SRG"&GeneSymbol %in% inputsleepGene[inputsleepGene$Tier==1,"GeneSymbol"]), 
               mapping=aes(x=RFsumF, 
                           y=factor(GeneSymbol), 
                           label=GeneSymbol), 
               fun=mean, geom="text", angle=45, size=2, color="#800000", hjust=-0.2, vjust=0) +
  geom_hline(yintercept=classsize, linetype="dashed", alpha=0.6, size=0.3) +
  scale_color_manual(values=c("#223E88", "#5BC0D2", "#0275D8", "#A9534F")) +
  scale_y_discrete("Gene ranking",
                   breaks=summarizeGeneScore$GeneSymbol[classsize], 
                   labels=as.character(classsize)) +
  xlab("Prediction score") +
  theme_classic() +
  theme(axis.text.y = element_text(size=10))

p2
ggsave(paste0(curPath,"prediction_score_",curDate,".pdf"), p2, width=4.5, height=4)

#ggsave(paste0(curPath,"prediction_score.pdf"), p2, width=9, height=2.5)


cutoffpoint <- 8.43
GWAS <- read.csv("./GWAS_mapped_gene.csv", header=T)
predictedSRG <- summarizeGeneScore %>%
  dplyr::filter(class==1) %>%
  left_join(GWAS, by="GeneSymbol") %>%
  dplyr::mutate(GWASann = apply(.[c("chronotype", "sleepDuration", "sleepiness", "insomnia")],1,function(x){sum(!is.na(x)&x>0)})) %>%
  dplyr::select(GeneSymbol, label, RFscore, RFrank, GWASann, everything()) %>%
  dplyr::arrange(-GWASann, RFrank)
write.table(predictedSRG,paste0("RFpredicted_",prefix,"_",nrow(predictedSRG),"_genes.txt"), sep="\t", row.names=FALSE)


## predicted genes and enriched pathway
f <- function(x) factor(x, levels = rev(unique(x)))
curDate <- '20210309'
dd <- read.table(paste0(curPath,"david_",curDate,"_combined_ann.txt"),header=T, sep="\t") %>%
#dd <- read.table("./bin_newFullGeneList/DAVID_combined_reactome5FAnn_202gene_20210106.txt", header=T, sep="\t") %>%
  dplyr::filter(FDR<=0.01, Term!='') %>%
  group_by(Group) %>%
  top_n(1, -FDR) %>%
  ungroup() %>%
  dplyr::arrange(PValue) %>%
  dplyr::mutate(log2pval=-log2(PValue),
                label = paste0(Term, " (", signif(PValue, digits=3), ")"))

pGene <- dd %>%
  dplyr::select(Term, label, Genes) %>%
  separate_rows(Genes) %>%
  dplyr::mutate(Genes=gsub(" ","", Genes)) %>%
  left_join(inputsleepGene, by=c("Genes"="GeneSymbol")) %>%
  left_join(select(summarizeGeneScore,"GeneSymbol", "RFrank"), by=c("Genes"="GeneSymbol")) %>%
  #dplyr::filter(Genes %in% predictedSRG$GeneSymbol) %>%
  #dplyr::mutate(Term=factor(Term, levels = termlevel)) %>%
  #dplyr::arrange(Term, Genes) %>%
  dplyr::mutate(Genes=factor(Genes, levels=unique(.$Genes))) %>%
  dplyr::mutate(color=ifelse(is.na(Tier),"blue","#333333"))
  #dplyr::mutate(color=ifelse(is.na(Tier)|Tier!=1,ifelse(is.na(Tier),"blue","black"),"#165178"))

pGene$group <- as.numeric(factor(pGene$Term))
for(i in 1:10){
  for(g in unique(pGene$Genes)){pGene[pGene$Genes==g,"group"]=min(pGene[pGene$Genes==g,"group"])}
  for(tt in unique(pGene$Term)){pGene[pGene$Term==tt,"group"]=min(pGene[pGene$Term==tt,"group"])}
}

GWAS <- read.csv("./GWAS_mapped_gene.csv", header=T) %>%
  dplyr::filter(GeneSymbol %in% pGene$Genes) %>%
  left_join(pGene, by=c("GeneSymbol"="Genes")) %>%
  dplyr::mutate(GWASann = apply(.[c("chronotype", "sleepDuration", "sleepiness", "insomnia")],1,function(x){sum(!is.na(x)&x>0)}))

selectedGWASterm <- GWAS %>% 
  group_by(Term) %>% 
  summarize(gene=n()) %>% 
  dplyr::filter(gene>1)
# subDavid <- subset(dd, dd$Term %in% GWAS$Term) %>%
#   distinct()
#subGene <- subset(pGene, pGene$Term %in% GWAS$Term) %>%
subDavid1 <- dd %>%
  dplyr::filter(Term %in% selectedGWASterm$Term) %>%
  dplyr::arrange(PValue) %>%
  dplyr::mutate(Term = factor(Term, levels=.$Term))
termlevel <- unique(subDavid1$Term)

subGene1 <- pGene %>%
  dplyr::filter(Term %in% selectedGWASterm$Term) %>%
  dplyr::filter(group==1, Term != "G alpha (q) signalling events", Term != "G alpha (i) signalling events") %>%
  dplyr::mutate(Term = factor(Term, levels=subDavid1$Term)) %>%
  left_join(select(GWAS,"GeneSymbol","GWASann"), by=c("Genes"="GeneSymbol")) %>%
  dplyr::arrange(Term, RFrank) %>%
  dplyr::mutate(label=factor(label, levels=unique(.$label)),
                GWASann=ifelse(is.na(GWASann),"No","Yes"), 
                Genes=factor(Genes, levels=unique(.$Genes))) %>%
  distinct()
  
textCol1 <- subGene1[match(levels(subGene1$Genes)[levels(subGene1$Genes)%in%subGene1$Genes], subGene1$Genes),]$color

p4 <- ggplot(subGene1) + 
  #geom_tile(aes(x=Genes, y=f(Term)), color="white", fill="orange") +
  geom_tile(aes(x=Genes, y=f(label), fill=GWASann), color="white") +
  scale_fill_manual(values=c("#6DA8DF", "#FF575B")) + 
  theme_minimal() +
  theme(axis.title = element_blank(),
        #axis.text.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.4, 
                                   colour=textCol1),
        legend.position = "right")
p4

subGene <- pGene %>%
  dplyr::filter(Term %in% selectedGWASterm$Term) %>%
  dplyr::filter(group>1 | Term == "G alpha (q) signalling events" | Term == "G alpha (i) signalling events") %>%
  dplyr::mutate(Term = factor(Term, levels=subDavid1$Term)) %>%
  left_join(select(GWAS,"GeneSymbol","GWASann"), by=c("Genes"="GeneSymbol")) %>%
  dplyr::arrange(Term, RFrank) %>%
  dplyr::mutate(label=factor(label, levels=unique(.$label)),
                GWASann=ifelse(is.na(GWASann),"No","Yes"), 
                Genes=factor(Genes, levels=unique(.$Genes))) %>%
  distinct()

textCol1 <- subGene[match(levels(subGene$Genes)[levels(subGene$Genes)%in%subGene$Genes], subGene$Genes),]$color

p5 <- ggplot(subGene) + 
  #geom_tile(aes(x=Genes, y=f(Term)), color="white", fill="orange") +
  geom_tile(aes(x=Genes, y=f(label), fill=GWASann), color="white") +
  scale_fill_manual(values=c("#6DA8DF", "#FF575B")) + 
  theme_minimal() +
  theme(axis.title = element_blank(),
        #axis.text.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.4, 
                                   colour=textCol1),
        legend.position = "right")
p5

plots <- align_plots(p4,
                     p5 + theme(legend.position="none"), 
                     align='v', 
                     axis='l')

prow <- plot_grid(
  plots[[1]],
  plots[[2]],
  ncol = 1,
  rel_heights = c(1.1,1)
)

pdf(paste0(curPath,"reactome_pathway_GWAS_",curDate,".pdf"), width=10,height=3.5)
plot_grid(prow)
dev.off()

## supplementary figures for pathway

#dd <- read.table("./bin_newFullGeneList/DAVID_combined_reactome5FAnn_202gene_20210106.txt", header=T, sep="\t") %>%
dd <- read.table(paste0(curPath,"david_",curDate,"_combined_ann.txt"),header=T, sep="\t") %>%
  dplyr::filter(FDR<=0.01, Term!='') %>%
  dplyr::arrange(PValue) %>%
  dplyr::mutate(log2pval=-log2(PValue))
dd$nMembers <- 0
for(i in 1:nrow(dd)){dd[i,"nMembers"]=nrow(dd[dd$Group==dd$Group[i],])}
dd <- dd %>%
  dplyr::filter(nMembers>=2) %>%
  dplyr::arrange(-EASE, -log2pval) %>%
  dplyr::mutate(Group = factor(Group, levels=unique(.$Group)),
                label = paste0(Term, " (", signif(PValue, digits=3), ")"))
termlevel <- unique(dd$Term)

pGene <- dd %>%
  dplyr::select(Term, label,Genes, Group) %>%
  separate_rows(Genes) %>%
  dplyr::mutate(Genes=gsub(" ","", Genes)) %>%
  left_join(inputsleepGene, by=c("Genes"="GeneSymbol")) %>%
  left_join(select(summarizeGeneScore,"GeneSymbol", "RFrank"), by=c("Genes"="GeneSymbol")) %>%
  #dplyr::filter(Genes %in% predictedSRG$GeneSymbol) %>%
  dplyr::mutate(Term=factor(Term, levels = termlevel)) %>%
  dplyr::arrange(Term, RFrank) %>%
  dplyr::mutate(Genes=factor(Genes, levels=unique(.$Genes))) %>%
  dplyr::mutate(color=ifelse(is.na(Tier),"blue","#333333"))

GWAS <- read.csv("./GWAS_mapped_gene.csv", header=T) %>%
  dplyr::filter(GeneSymbol %in% pGene$Genes) %>%
  left_join(pGene, by=c("GeneSymbol"="Genes")) %>%
  dplyr::mutate(GWASann = apply(.[c("chronotype", "sleepDuration", "sleepiness", "insomnia")],1,function(x){sum(!is.na(x)&x>0)}))

subGene1 <- pGene %>%
  #dplyr::filter(group==1) %>%
  left_join(select(GWAS,"GeneSymbol","GWASann"), by=c("Genes"="GeneSymbol")) %>%
  dplyr::mutate(GWASann=ifelse(is.na(GWASann),"No","Yes"), 
                Genes=factor(Genes, levels=unique(pGene$Genes))) %>%
  distinct()
subDavid1 <- dd %>%
  dplyr::filter(Term %in% subGene1$Term)

textCol1 <- subGene1[match(levels(subGene1$Genes)[levels(subGene1$Genes)%in%subGene1$Genes], subGene1$Genes),]$color
p8 <- ggplot(subGene1) + 
  #geom_tile(aes(x=Genes, y=f(Term)), color="white", fill="orange") +
  geom_tile(aes(x=Genes, y=f(label), fill=GWASann), color="white") +
  scale_fill_manual(values=c("#6DA8DF", "#FF575B")) + 
  #facet_wrap(.~Group, scale="free", ncol = 2) +
  facet_grid(Group~., scales = "free", space='free') +
  theme_minimal() +
  theme(axis.title = element_blank(),
        #axis.text.y = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.4, size=9, 
                                   colour=textCol1),
        strip.text = element_blank())
p8

#dd <- read.table("./bin_newFullGeneList/DAVID_combined_reactome5FAnn_202gene_20210106.txt", header=T, sep="\t") %>%
dd <- read.table(paste0(curPath,"david_",curDate,"_combined_ann.txt"),header=T, sep="\t") %>%
  dplyr::filter(FDR<=0.01, Term!='') %>%
  dplyr::arrange(PValue) %>%
  dplyr::mutate(log2pval=-log2(PValue))
dd$nMembers <- 0
for(i in 1:nrow(dd)){dd[i,"nMembers"]=nrow(dd[dd$Group==dd$Group[i],])}
dd <- dd %>%
  dplyr::filter(nMembers<2) %>%
  dplyr::arrange(-log2pval) %>%
  dplyr::mutate(Group = factor(Group, levels=unique(.$Group)),
                label = paste0(Term, " (", signif(PValue, digits=3), ")"))
termlevel <- unique(dd$Term)

pGene <- dd %>%
  dplyr::select(Term, label,Genes, Group) %>%
  separate_rows(Genes) %>%
  dplyr::mutate(Genes=gsub(" ","", Genes)) %>%
  left_join(inputsleepGene, by=c("Genes"="GeneSymbol")) %>%
  left_join(select(summarizeGeneScore,"GeneSymbol", "RFrank"), by=c("Genes"="GeneSymbol")) %>%
  #dplyr::filter(Genes %in% predictedSRG$GeneSymbol) %>%
  dplyr::mutate(Term=factor(Term, levels = termlevel)) %>%
  dplyr::arrange(Term, RFrank) %>%
  dplyr::mutate(Genes=factor(Genes, levels=unique(.$Genes))) %>%
  dplyr::mutate(color=ifelse(is.na(Tier),"blue","#333333"))

GWAS <- read.csv("./GWAS_mapped_gene.csv", header=T) %>%
  dplyr::filter(GeneSymbol %in% pGene$Genes) %>%
  left_join(pGene, by=c("GeneSymbol"="Genes")) %>%
  dplyr::mutate(GWASann = apply(.[c("chronotype", "sleepDuration", "sleepiness", "insomnia")],1,function(x){sum(!is.na(x)&x>0)}))

subGene1 <- pGene %>%
  #dplyr::filter(group==1) %>%
  left_join(select(GWAS,"GeneSymbol","GWASann"), by=c("Genes"="GeneSymbol")) %>%
  dplyr::mutate(GWASann=ifelse(is.na(GWASann),"No","Yes"), 
                Genes=factor(Genes, levels=unique(pGene$Genes))) %>%
  distinct()
subDavid1 <- dd %>%
  dplyr::filter(Term %in% subGene1$Term)

textCol1 <- subGene1[match(levels(subGene1$Genes)[levels(subGene1$Genes)%in%subGene1$Genes], subGene1$Genes),]$color
p9 <- ggplot(subGene1) + 
  #geom_tile(aes(x=Genes, y=f(Term)), color="white", fill="orange") +
  geom_tile(aes(x=Genes, y=f(label), fill=GWASann), color="white") +
  scale_fill_manual(values=c("#6DA8DF", "#FF575B")) + 
  #facet_wrap(.~Group, scale="free", ncol = 2) +
  facet_grid(Group~., scales = "free", space='free') +
  theme_minimal() +
  theme(axis.title = element_blank(),
        #axis.text.y = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.4, size=9, 
                                   colour=textCol1),
        strip.text = element_blank())

plots <- align_plots(p8 + theme(legend.position="none"),
                     p9, 
                     align='v', 
                     axis='l')

prow <- plot_grid(
  plots[[1]],
  plots[[2]],
  ncol = 1,
  rel_heights = c(2,0.8)
)

prow

#pdf("reactome_pathway_202gene_allTerms_20210306.pdf", width=18, height=9)
pdf(paste0(curPath,"reactome_pathway_GWAS_allTerms_",curDate,".pdf"), width=18,height=9)
prow
dev.off()

## plot figure for evidence factors
PATH <- read.table("fList_forPlotting.txt",header=F) %>%
  magrittr::set_colnames(c("Source", "ID","Input","Models","DataType")) %>%
  dplyr::mutate(rID=gsub("-",".",ID)) %>%
  dplyr::select(-ID) %>%
  dplyr::rename("ID"="rID") %>%
  #dplyr::mutate(maxEF=NA, upper99.EF=NA)
  dplyr::mutate(EFscore=NA)

p <- list()
for(i in 1:nrow(PATH)){
  source <- PATH[i,"Source"]
  ID <- PATH[i,"ID"]
  fn <- paste0(pdir,PATH[i,"Input"])
  model <- PATH[i,"Models"]
  dataType <- PATH[i,"DataType"]
  dat <- readFN(fn)
  dat[dat$value==Inf,"value"] <- max(dat[dat$value<Inf,"value"])+0.1
  dat <- convertGeneNameToHGNC(dat, model)
  dat <- labelData(dat, sleepGeneList=inputsleepGeneList)
  
  ef <- calEF_bin(dat,dataType)
  datEF <- matchEF_bin(dat,ef)
  dat <- full_join(dat, datEF, by="GeneSymbol")
  
  p1 <- ggplot(dat) +
    geom_density_ridges2(aes(x=value,y=label,color=label,fill=label)) +
    scale_color_manual(values=c('#0560B0','#A05021')) +
    scale_fill_manual(values=c('dodgerblue','coral')) +
    ylab("density") +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position="none",
          plot.margin = unit(c(-0.5, 3, 3, 3), "point")
    ) +
    coord_cartesian(xlim = c(0,ceiling(max(dat$value)))) 
  p2 <- ggplot(dat) +
    geom_line(aes(x=value, y=evidenceFactors), color="#23395D") +
    #geom_rug(subset(dat,dat$label=="nonSRG"), mapping=aes(x=value, y=evidenceFactors), sides='t', alpha=0.8, color="#A0A0A0") +
    #geom_rug(subset(dat,dat$label=="SRG"), mapping=aes(x=value, y=evidenceFactors), sides='b', alpha=0.8, color="#E36927") +
    ylim(0,4) +
    ylab("evidence\nfactors") +
    xlab("-log2pval(JTK_cycle)") +
    theme_minimal() +
    theme(legend.position="none",
          plot.margin = unit(c(5, 3, 3, 3), "point")
    ) +
    coord_cartesian(xlim = c(0,ceiling(max(dat$value)))) 
  
  p[[i]] <- plot_grid(plotlist=list(p1,p2), ncol=1, align='v', rel_heights = c(1,1))
}

pdf("EF_plot.pdf",6,2.5)
ggarrange(p[[1]], p[[2]], ncol = 2, nrow = 1)
dev.off()

## which feature enriched for sleep genes
## second round of feature selection - after remove to same number of genes
# cleanRawValuetable <- read.csv("./bin_newFullGeneList/fList_cleanRawValue_withJI_20201229.csv", header=T)
# newcleanEFtable <- cleanRawValuetable[,c("GeneSymbol","label")]
# pdf(paste0(pdir,curPath,"EFbin_selectedDat_test3_",curDate,".pdf"),8,3)
# p <- list()
# for(i in 3:ncol(cleanRawValuetable)){
#   ID <- colnames(cleanRawValuetable)[i]
#   if(i%%100==0){print(c("processing column: ", i, ID))}
#   
#   dataType <- PATH[PATH$ID==ID,"DataType"]
#   dat <- cleanRawValuetable[,c(1,2,i)] %>%
#     drop_na() %>%
#     magrittr::set_colnames(c("GeneSymbol", "label", "value"))
#   
#   if(nrow(subset(dat, dat$label=="SRG"))<minInputSleepGeneCuttoff){
#     print(paste(c("skip as less than 25% sleeep genes found in the data:", ID)))
#     next
#   }
#   
#   ef <- calEF_bin(dat,dataType)
#   if(nrow(ef) <= 1){
#     print(paste(c("skip as equal data value:", ID)))
#     next
#   }
#   
#   # p[[1]] <- drawDensityPlot(dat, nCol=3)
#   # p[[2]] <- drawEFPlot(ef, name = ID)
#   # do.call("grid.arrange", c(p, ncol=2))
#   
#   datEF <- matchEF_bin(dat,ef)
#   
#   newcleanEFtable <- newcleanEFtable %>%
#     left_join(datEF, by=c("GeneSymbol")) %>%
#     dplyr::rename(!!ID := evidenceFactors)
#   
# }
# dev.off()
# 
# zz <- apply(newcleanEFtable[,3:ncol(newcleanEFtable)],2,function(x){sort(na.omit(x),decreasing = TRUE)[1]})
# yy <- zz[zz>=3]
# newRawValuetable <- cleanRawValuetable %>%
#   dplyr::select(c("GeneSymbol","label",c(names(yy))))
# 
# newEFtable <- newcleanEFtable %>%
#   dplyr::select(c("GeneSymbol","label",c(names(yy))))
# write.csv(newRawValuetable, paste0(curPath,"fList_clean2RawValue_",curDate,".csv"), row.names=FALSE)
# write.csv(newEFtable, paste0(curPath,"fList_clean2EFtable_",curDate,".csv"), row.names=FALSE)
# 
# newEFtable <- read.csv(paste0(curPath,"fList_clean2EFtable_",curDate,".csv"),header=T)
# colnames(newEFtable) <- gsub(".log2foldchange.txt","",colnames(newEFtable))
# EF_sg <- subset(newEFtable, newEFtable$label=="SRG")
# aa <- apply(EF_sg[,3:ncol(EF_sg)],2,function(x){sort(na.omit(x),decreasing = TRUE)[1]})
# aa <- aa[aa>=3]
# 
# feature_matched_name <- read.csv("./bin_newFullGeneList/EF3_ann.csv",header=TRUE) %>%
#   dplyr::filter(short_name!="gestational.preterm_E-MTAB-4860") %>%
#   dplyr::filter(short_name!="acrolein10mM_E-MEXP-1599") %>%
#   dplyr::filter(short_name!="astrocyte.overexpressing.4.gene_E???MTAB???4771") %>%
#   dplyr::filter(short_name!="GTEx_AdrenalGland")
# 
# selectCol <- feature_matched_name %>%
#   dplyr::filter(datatype != "remove") %>%
#   dplyr::filter(ID %in% gsub(".log2foldchange.txt", "", colnames(newEFtable))) %>%
#   dplyr::mutate(maxEF=aa[.$ID]) %>%
#   dplyr::arrange(-maxEF)
# 
# subFeature <- feature_matched_name %>%
#   #dplyr::filter(ID %in% selectCol$ID) %>%
#   dplyr::mutate(maxEF=aa[.$ID]) %>%
#   dplyr::arrange(-maxEF) %>%
#   group_by(Source, datatype) %>%
#   top_n(3, maxEF) %>%
#   ungroup() %>%
#   as.data.frame()
# subEF <- newEFtable[,c("GeneSymbol", "label", subFeature$ID)]
# subEF_sg <- EF_sg[c("GeneSymbol", "label", subFeature$ID)]
# 
# for(i in 1:nrow(subFeature)){
#   if(subFeature[i,"ID"]!=subFeature[i,"short_name"]){
#     subEF <- subEF %>%
#       dplyr::rename(!!subFeature[i,"short_name"] := subFeature[i,"ID"])
#     subEF_sg <- subEF_sg %>%
#       dplyr::rename(!!subFeature[i,"short_name"] := subFeature[i,"ID"])
#   }
# }
# 
# tmp_EFsg <- subEF_sg[,-c(1,2)]
# rownames(tmp_EFsg) <- subEF_sg[,"GeneSymbol"]
# hc <- hclust(dist(tmp_EFsg), "ward.D2")
# geneOrder <- subEF_sg[hc$order, "GeneSymbol"]
# hc2 <- hclust(dist(t(tmp_EFsg)), "ward.D2")
# featureOrder <- colnames(tmp_EFsg)[hc2$order]
# 
# subEF <- subEF %>%
#   melt(., id.vars=c("GeneSymbol", "label")) %>%
#   left_join(selectCol, by=c("variable"="short_name.2")) %>%
#   drop_na() %>%
#   dplyr::mutate(datatype=factor(datatype, levels=unique(selectCol$datatype))) %>%
#   dplyr::mutate(variable=factor(variable, levels=rev(unique(selectCol$short_name))))
# 
# subEF_sg <- subEF_sg %>%
#   melt(., id.vars=c("GeneSymbol", "label")) %>%
#   left_join(selectCol, by=c("variable"="short_name.2")) %>%
#   dplyr::mutate(datatype=factor(datatype, levels=unique(selectCol$datatype))) %>%
#   #dplyr::arrange(EF, datatype) %>%
#   dplyr::mutate(variable=factor(variable, levels=rev(unique(selectCol$short_name))))
# 
# subFeaturePlot <- subFeature %>%
#   dplyr::mutate(datatype=factor(datatype, levels=unique(selectCol$datatype)),
#                 short_name.2=factor(short_name.2, levels=unique(selectCol$short_name.2)))
# color <- list("human"="orange",
#               "mouse"="dodgerblue")

# p1 <- ggplot() +
#   #geom_density_ridges2(subEF, mapping=aes(x=value, y=variable), color=color[[s]], fill=color[[s]], scale = .85) +
#   geom_boxplot(subEF,
#                mapping=aes(x=value, y=variable, color=datatype),
#                outlier.shape = NA) +
#   geom_jitter(subset(subEF,subEF$label=="SRG"), 
#               mapping=aes(x=value, y=variable, color=datatype),
#               size = 0.2,
#               width = 0.1) +    
#   geom_vline(xintercept = 3, linetype="dashed", color="#777777") +
#   scale_y_discrete(position = "right") +
#   ylab("density") +
#   xlab("evidence factors") +
#   facet_grid(datatype+Source~., scales="free_y", space="free", switch = "y") +
#   #facet_grid(rows=vars(datatype), scales="free_y", space="free", switch = "y") +
#   theme_light() +
#   theme(
#     axis.title.y = element_blank(),
#     strip.background = element_rect(color="white", fill="#EEEEEE"),
#     panel.border = element_rect(color="white"),
#     strip.text.y.left = element_text(color="#333333", angle=0, size=12),
#     strip.placement = "outside",
#     #strip.text.y = element_text(color="#333333", angle=0, hjust = 0, size=12),
#     legend.position = 'none'
#   )
# ggsave(p1, filename = paste0("feature_topGeneEF_20210118.pdf"), width=10, height=10)

#EFtable <- read.csv("./bin_newFullGeneList/fList_EFtable_20201229.csv",header=T)
EFtable <- read.csv(paste0(curPath, "fList_EFtable_", curDate, ".csv"), header=T)
zz <- apply(EFtable[,3:ncol(EFtable)],2,function(x){max(x, na.rm = TRUE)})
yy <- zz[zz>=3]

feature_matched_name <- read.csv("./bin_newFullGeneList/EF3_ann.csv",header=TRUE)

topEF <- as.data.frame(yy) %>%
  dplyr::mutate(ID=gsub(".log2foldchange.txt", "",rownames(.))) %>%
  left_join(feature_matched_name, by="ID") %>%
  dplyr::filter(datatype!="remove") %>%
  dplyr::filter(short_name!="gestational.preterm_E-MTAB-4860") %>%
  dplyr::select(-maxEF) %>%
  dplyr::rename(maxEF=yy)
  
subFeature <- topEF %>%
  dplyr::arrange(-maxEF) %>%
  group_by(Source, datatype) %>%
  top_n(4, maxEF) %>%
  ungroup() %>%
  as.data.frame()

subFeaturePlot <- topEF %>%
  dplyr::filter(ID %in% subFeature$ID) %>%
  dplyr::arrange(-maxEF) %>%
  dplyr::mutate(short_name.2 = factor(short_name.2, levels=unique(subFeature$short_name.2)),
                datatype = factor(datatype, levels=unique(.$datatype)))

p1 <- ggplot(subFeaturePlot) +
  geom_col(mapping=aes(x=short_name.2, y=maxEF, fill=Source)) +
  geom_hline(yintercept = 3, linetype="dashed", color="#333333") +
  scale_fill_manual(values=c("#FC766A","#5B84B1")) +
  ylab("maximum evidence factors") +
  xlab("data metrics") +
  facet_grid(.~datatype, scales = "free_x", space="free") +
  #facet_grid(rows=vars(datatype), scales="free_y", space="free", switch = "y") +
  theme_light() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle=45, hjust = 1, vjust=1, size=10),
    strip.background = element_rect(color="white", fill="white"),
    panel.border = element_rect(color="white"),
    strip.text.x = element_text(color="#333333", angle=0, size=12),
    legend.text = element_text(size=10),
    legend.position = "top"
  )
#pdf("feature_topGeneEF_20210301.pdf", width = 12, height = 5)
pdf(paste0(curPath, "feature_topGeneEF_", curDate, ".pdf"), width = 12, height = 5)
p1 + theme(plot.margin = margin(0.5, 0.5, 0.5, 1, "cm"))
dev.off()


## gene-level with significant features
colnames(EFtable) <- gsub(".log2foldchange.txt", "", colnames(EFtable))

f_3 <- function(x) {ifelse(x>=3, x, NA)}

subEF_sg <- EFtable[c("GeneSymbol", "label", subFeaturePlot$ID)]
tmp_EFsg <- subset(subEF_sg, subEF_sg$label=="SRG") %>%
  dplyr::select(-GeneSymbol, -label)
rownames(tmp_EFsg) <- subset(subEF_sg, subEF_sg$label=="SRG")$GeneSymbol
hc <- hclust(dist(tmp_EFsg), "ward.D2")
geneOrder <- rownames(tmp_EFsg)[hc$order]

bb <- subset(subEF_sg, subEF_sg$label=="SRG") %>%
  dplyr::mutate_if(is.numeric, f_3) %>%
  gather(variable, value, -GeneSymbol, -label ) %>%
  drop_na() %>%
  dplyr::mutate(GeneSymbol=factor(GeneSymbol, levels=rev(geneOrder))) %>%
  left_join(select(subFeaturePlot, datatype, ID, short_name.2, maxEF), by=c("variable"="ID")) %>%
  dplyr::mutate(short_name.2=factor(short_name.2, levels=rev(unique(subFeaturePlot$short_name.2))))
  
p2 <- ggplot(bb) +
  geom_tile(aes(x=GeneSymbol, y=short_name.2, fill=maxEF), color="white") +
  #facet_grid(rows=vars(datatype), cols=vars(Source), scales="free_y", space="free_y", switch = "y") +
  #facet_grid(datatype+Source~., scales="free_y", space="free_y", switch = "y") +
  facet_grid(datatype~., scales = "free", space="free", switch = "y") +
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE, n.dodge=2))+
  #facet_wrap(datatype~., scales = "free", ncol = 2) +
  #scale_y_discrete(position = "right") +
  #scale_x_discrete(position = "top") +
  geom_vline(xintercept = seq(10,length(unique(bb$GeneSymbol)),10)-0.5, linetype="longdash", color="#333333") +
  theme_light() +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle=90, hjust = 1, vjust=0.4, size=9),
    axis.text.y = element_text(size=10),
    strip.background = element_rect(color="white", fill="#EEEEEE"),
    panel.border = element_rect(color="white"),
    strip.text.y.left = element_text(color="#333333", angle=0, size=12),
    #strip.text.x.top = element_text(color="#333333", angle=0, size=12),
    strip.placement = "outside",
    legend.position = "right"
  )
p2 + theme(plot.margin = margin(0.5, 0.5, 0.5, 1, "cm"))

#pdf("selected_feature_heatmap_20210302.pdf", width=16, height=14)
pdf("selected_feature_heatmap_20210303.pdf", width=14, height=8)
p2 + theme(plot.margin = margin(0.5, 0.5, 0.5, 1, "cm"))
dev.off()

## feature importances from random forest
feature_matched_name <- read.csv("./bin_newFullGeneList/EF3_ann.csv",header=TRUE)
PATH <- read.table(paste0(pdir,"combine_fList_updateGeneList.txt"),header=F) %>%
  magrittr::set_colnames(c("Source", "ID","Input","Models","DataType")) %>%
  dplyr::mutate(rID=gsub("-",".",ID)) %>%
  dplyr::mutate(Source = replace(Source, str_detect(Input, "EBI"), "log2FoldChange")) %>%
  dplyr::mutate(Source = replace(Source, str_detect(Input, "GSE"), "log2FoldChange")) %>%
  dplyr::mutate(Source = ifelse(Source == "AllenBrainAtlas", "transcriptAbundance", Source)) %>%
  dplyr::mutate(Source = ifelse(Source == "transcriptAbundance", "transcript\nAbundance", Source)) %>%
  dplyr::select(-ID) %>%
  dplyr::rename("ID"="rID") %>%
  dplyr::mutate(ID=gsub(".log2foldchange.txt","",ID)) %>%
  left_join(select(feature_matched_name, "ID", "full_name", "short_name.2"), by=c("ID")) %>%
  dplyr::mutate(short_name=ifelse(is.na(short_name.2),ID, short_name.2), full_name=ifelse(is.na(full_name),ID,full_name))


nFeatures <- read.table("./bin_newFullGeneList/feature_importance_20201229.txt", header=F) %>%
  separate_rows(., V3, convert = TRUE) %>%
  drop_na() %>%
  magrittr::set_colnames(c("ID", "trainRatio", "score")) %>%
  #dplyr::filter(score>=0.01) %>%
  dplyr::mutate(ID=factor(ID), trainRatio=factor(trainRatio)) %>%
  group_by(ID) %>%
  summarize(count=sum(score>=0.01), mean=median(score)) %>%
  dplyr::filter(count>=8)

features_importance <- read.table("./bin_newFullGeneList/feature_importance_20201229.txt", header=F) %>%
  separate_rows(., V3, convert = TRUE) %>%
  drop_na() %>%
  magrittr::set_colnames(c("ID", "trainRatio", "score")) %>%
  dplyr::mutate(ID=gsub(".log2foldchange.txt","",ID)) %>%
  left_join(PATH, by="ID") %>%
  left_join(nFeatures, by="ID") %>%
  dplyr::filter(ID %in% nFeatures$ID) %>%
  dplyr::select(ID, full_name, short_name, trainRatio, score, Source, Models, count, mean) %>%
  dplyr::mutate(Source=ifelse(is.na(Source), "Gene set collections (JI scoring methods)", Source),
                Models=ifelse(is.na(Models), "human", Models),
                full_name=ifelse(is.na(full_name), ID, full_name),
                short_name=ifelse(is.na(short_name), ID, short_name)) %>%
  dplyr::arrange(Source, Models, -mean) %>%
  dplyr::mutate(Source=factor(Source, levels=unique(as.character(.$Source))), 
                ID=factor(ID, levels=unique(as.character(.$ID))),
                full_name=factor(full_name, levels=rev(unique(as.character(.$full_name)))),
                short_name=factor(short_name, levels=rev(unique(as.character(.$short_name)))),
                trainRatio=factor(trainRatio))

p1 <- ggplot(subset(features_importance, features_importance$Source=="Gene set collections (JI scoring methods)")) + 
  geom_point(mapping=aes(x=score, y=short_name), color="#A0A0A0", position=position_dodge(width=0.5)) +
  geom_point(mapping=aes(x=mean, y=short_name), color="#B72803", size=1.5) +
  #scale_color_brewer(palette = "Blues") +
  facet_grid(.~Source, scales="free", space="free") +
  theme_bw() +
  theme(#axis.text.x = element_text(angle=45, hjust=1),
    axis.text.x = element_text(size=11),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    strip.text = element_text(size=13),
    strip.background = element_rect(color="white", fill="#EEEEEE"))

p2 <- ggplot(subset(features_importance, features_importance$Source!="Gene set collections (JI scoring methods)")) + 
  geom_point(mapping=aes(x=score, y=short_name), color="#A0A0A0", position=position_dodge(width=0.5)) +
  geom_point(mapping=aes(x=mean, y=short_name), color="#B72803", size=1.5) +
  #scale_color_brewer(palette = "Blues") +
  #facet_grid(.~Models+Source, scales="free", space="free") +
  facet_grid(Source~Models, scales="free", space="free_y") +
  theme_bw() +
  theme(#axis.text.x = element_text(angle=45, hjust=1),
    axis.text.x = element_text(size=11),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    strip.text = element_text(size=13),
    strip.text.y = element_blank(),
    strip.background = element_rect(color="white", fill="#EEEEEE"))

plots <- align_plots(p1 + theme(legend.position="none"),
                     p2 + theme(legend.position="none"), 
                     align='v', 
                     axis='l')

prow <- plot_grid(
  #p1 + theme(legend.position="none"),
  plots[[1]],
  p2 + theme(legend.position="none"),
  labels = c("A", "B"),
  ncol = 1,
  rel_heights = c(.8,1)
)

legend <- get_legend(
  # create some space to the left of the legend
  p1 + theme(legend.box.margin = margin(0, 0, 0, 6))
)

pdf("feature_importance_20201229.pdf", 9,8)
plot_grid(prow, legend, rel_widths = c(3, .4))
dev.off()

## plot EF for all selected features
cleanRawValuetable <- read.csv("./bin_newFullGeneList/fList_cleanRawValue_withJI_20201229.csv", header=T)
feature_matched_name <- read.csv("./bin_newFullGeneList/EF3_ann.csv",header=TRUE)

p <- list()
j <- 1
pdf("EF_new.pdf", 16,9)
for(i in 3:ncol(cleanRawValuetable)){
  curFeature <- feature_matched_name[feature_matched_name$ID %in% colnames(cleanRawValuetable)[i], ]
  ID <- curFeature$short_name.2
  model <- curFeature$ID
  dat <- cleanRawValuetable[,c(1,2,i)] %>%
    drop_na()
  colnames(dat)[3] <- "value"
  
  ef <- calEF_bin(dat, ID)
  if(nrow(ef)==0){
    print(paste(ID, "no EF output"))
    next
  }
  datEF <- matchEF_bin(dat,ef)
  dat <- full_join(dat, datEF, by="GeneSymbol")
  
  maxEFrange <- ifelse(max(ef$evidenceFactors)<=5,5, ifelse(max(ef$evidenceFactors)<=8,8,ceiling(max(dat$evidenceFactors))))
  p1 <- ggplot(dat) +
    geom_density_ridges2(aes(x=value,y=label,color=label,fill=label)) +
    scale_color_manual(values=c('#0560B0','#A05021')) +
    scale_fill_manual(values=c('dodgerblue','coral')) +
    ylab("density") +
    ggtitle(ID) +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position="none",
          plot.margin = unit(c(-0.5, 3, 3, 3), "point")
    ) +
    coord_cartesian(xlim = c(0,max(dat$value)))
  p2 <- ggplot(dat) +
    geom_line(aes(x=value, y=evidenceFactors), color="#23395D") +
    geom_hline(yintercept = 3, linetype="dashed", color="#777777") +
    #geom_rug(subset(dat,dat$label=="nonSRG"), mapping=aes(x=value, y=evidenceFactors), sides='t', alpha=0.8, color="#A0A0A0") +
    #geom_rug(subset(dat,dat$label=="SRG"), mapping=aes(x=value, y=evidenceFactors), sides='b', alpha=0.8, color="#E36927") +
    ylim(0,maxEFrange) +
    ylab("evidence\nfactors") +
    theme_minimal() +
    theme(legend.position="none",
          plot.margin = unit(c(5, 3, 3, 3), "point"),
          axis.title.x = element_blank()
    ) +
    coord_cartesian(xlim = c(0,max(dat$value)))
  
  p[[j]] <- plot_grid(plotlist=list(p1,p2), ncol=1, align='v', rel_heights = c(1,1))
  
  if(j%%16==0){
    do.call("grid.arrange", c(p, ncol=4))
    p <- list()
    j <- 0
  }
  j <- j+1
}
if(length(p)>0){
  do.call("grid.arrange", c(p, ncol=4))
}
dev.off()

