## script to process random forest prediction
## summarize gene prediction

curPath <- "./OUT/"
prefix <- "n100_20210401"
class1cutoff <- 0.413

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

summarizeGeneScore <- combineGeneScore %>%
  dplyr::filter(RFsumF>0.1) %>%
  group_by(GeneSymbol, label) %>%
  dplyr::summarise(rawScore=mean(RFsumF), RFscore=max(ratio*10)+mean(RFsumF)) %>%
  #dplyr::summarise(rawScore=mean(RFsumF), RFscore=n(RFsumF[RFsumF>0.1]) + mean(RFsumF)) %>%
  ungroup() %>%
  dplyr::arrange(-RFscore) %>%
  dplyr::mutate(RFrank=as.numeric(rownames(.))) %>%
  dplyr::mutate(class=10-ifelse(rawScore>=class1cutoff, ceiling(RFscore),floor(RFscore)))
table(summarizeGeneScore$label, summarizeGeneScore$class)
write.table(summarizeGeneScore,paste0(curPath, "RFsummarizeOutput",prefix,".txt"), sep="\t", row.names=FALSE)

geneScore <- combineGeneScore %>%
  dplyr::select(GeneSymbol, label, RFsumF, ratio) %>%
  dplyr::mutate(RFsumF=ifelse(RFsumF>=0.1, RFsumF, NA)) %>%
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
  geom_point(color="grey", size=0.25, alpha=0.8, shape=16) +
  stat_summary(fun="mean", geom="point", aes(color=sublabel, shape=sublabel), alpha=0.9, size=1) +
  stat_summary(subset(geneScore,geneScore$label=="SRG"&GeneSymbol %in% inputsleepGene[inputsleepGene$Tier==1,"GeneSymbol"]), 
               mapping=aes(x=RFsumF, 
                           y=factor(GeneSymbol), 
                           label=GeneSymbol), 
               fun=mean, geom="text", angle=45, size=3, color="#800000", hjust=-0.2, vjust=0) +
  geom_hline(yintercept=classsize, linetype="dashed", alpha=0.6, size=0.3) +
  scale_color_manual(values = c("#223E88", "#5BC0D2", "#FFC107", "#A9534F")) +
  scale_fill_manual(values = c("#223E88", "#5BC0D2", "#FFC107", "#A9534F")) +
  scale_shape_manual(values = c(21,15,1,19)) +
  scale_y_discrete("Gene ranking",
                   breaks=summarizeGeneScore$GeneSymbol[classsize], 
                   labels=as.character(classsize)) +
  xlab("Prediction score") +
  theme_classic() +
  theme(axis.text.y = element_text(size=10))
ggsave(paste0(curPath,"prediction_score_",curDate,".pdf"), p2, width=4, height=3.5, useDingbats = FALSE)
