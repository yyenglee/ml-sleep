pdir = "~/workDIR/sleep/"
setwd("~/workDIR/sleep")

library(ggplot2)
library(dplyr)
library(reshape)
library(magrittr)
library(data.table)
library(magrittr)
library(RColorBrewer)
library(cowplot)
library(ROCR)
options(stringsAsFactors = FALSE)

folderName <- "./"
curDate <- "n100_sss_20210506"

algo_name <- c(
  "NaiveBayes" = "Bayes", 
  "LogisticRegression"="LR", 
  "LinearSVM" = "SVM", 
  "DecisionTree" = "DT", 
  "RandomForest" = "RF", 
  "AdaptiveBoosting" = "AB", 
  "NeuralNetwork" = "NN",
  "ensembleNeuralNetwork" = "eNN")

algoLevel <- c("Naive Bayes","Logistic Regression","Linear SVM","Decision Tree",
               "Random Forest","Adaptive Boosting","Neural Network",
               "ensemble Neural Network")
algoLevel <- c("Bayes", "LR", "SVM", "DT", "RF", "AB", "NN", "eNN")
trainTestRatio <- c("0.2|0.8","0.3|0.7","0.4|0.6", "0.5|0.5", "0.6|0.4", "0.7|0.3", "0.8|0.2")
labelLevel <- c("sleep gene", "random label")

ensemblMLSG <- read.csv(paste0("ensembleML_matrixSleepGene_",curDate,".csv"),header=T) %>%
  magrittr::set_colnames(c("n","algorithm","ratio","TN", "FP", "FN","TP")) %>%
  dplyr::filter(TN>0, algorithm %in% names(algo_name)) %>%
  dplyr::mutate(nPredictedPositive=TP+FP, Precision=ifelse(TP>0,TP/(TP+FP),0), Recall=ifelse(TP>0, TP/(TP+FN),0)) %>%
  dplyr::mutate(F1score=ifelse((Precision+Recall)>0, (2*Precision*Recall)/(Precision+Recall), 0)) %>%
  dplyr::mutate(label="sleep gene")
ensemblMLRG <- read.csv(paste0(folderName, "ensembleML_matrixRandomLabel_",curDate,".csv"),header=T) %>% 
  magrittr::set_colnames(c("n","algorithm","ratio","TN", "FP", "FN","TP")) %>%
  dplyr::filter(TN>0, algorithm %in% names(algo_name)) %>%
  dplyr::mutate(nPredictedPositive=TP+FP, Precision=ifelse(TP>0,TP/(TP+FP),0), Recall=ifelse(TP>0, TP/(TP+FN),0)) %>%
  dplyr::mutate(F1score=ifelse((Precision+Recall)>0, 2*(Precision*(Recall))/(Precision+Recall), 0)) %>%
  dplyr::mutate(label="random label")

ensemblML <- bind_rows(ensemblMLSG, ensemblMLRG) %>%
  dplyr::mutate(algo_short=apply(., 1, function(x){algo_name[x[2]]})) %>%
  dplyr::mutate(label=factor(label, levels=labelLevel),
                train=as.numeric(gsub("\\|.*", "", .$ratio))) %>%
  dplyr::mutate(ratio=paste(.$train,1-.$train,sep="|"),
                algo_short=factor(algo_short, levels=algoLevel))

rdEnsemblML <- subset(ensemblML, ensemblML$label %in% labelLevel & ratio %in% trainTestRatio) %>%
  dplyr::select(algorithm, ratio, label, Precision, Recall, F1score, algo_short) %>%
  melt(., id.vars = c("algorithm", "ratio", "label", "algo_short"))

rdSum <- rdEnsemblML %>% 
  group_by(algo_short, label, variable) %>% 
  summarize(value=mean(value))

p1 <- ggplot(rdEnsemblML) + 
  geom_point(aes(x=algo_short, y=value, color=algo_short, shape=ratio), size=0.5, position=position_dodge(width=0.75)) + 
  geom_point(data=rdSum,  mapping=aes(x = algo_short, y = value, fill=algo_short), shape=21, color="white", alpha=0.9, size=3, stroke=0) +
  scale_color_brewer(palette="Set2") +
  scale_fill_brewer(palette="Dark2") +
  scale_shape_manual(values=c(2:8)) +
  #facet_grid(label~variable) +
  facet_grid(variable~label, ) +
  labs(shape="train|test ratio") +
  theme_bw() +
  theme(axis.title = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size=15),
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle = 45, hjust=1, vjust=1, size=12),
        #axis.text.x = element_blank(),
        strip.background = element_rect(color="white", fill="#EEEEEE"),
        strip.text.y.right = element_text(angle = 0), 
        panel.border = element_blank(),
        legend.box = "horizontal",
        legend.position = "none"
  ) +
  guides(shape = guide_legend(override.aes = list(size = 3)),
         color = guide_legend(override.aes = list(size = 3)))
p1
ggsave(paste0(folderName,"SRG_randomlabel_bayesPC.pdf"), p1, width=7, height=7, useDingbats=FALSE)

## make precision-recall plot
curDate <- "sss_20210506"
inputsleepGeneList <- read.table("inputSleepGeneList.txt", header=FALSE)$V1
predScore <- read.csv(paste0(folderName,"predScore_",curDate,".csv"),header=T)
predScore <- labelData(predScore, sleepGeneList = inputsleepGeneList)

seqList <- seq(0.00001,1,0.001)

prcTable <- data.frame(
        algo=character(),
        TruePositiveRate=numeric(), 
        FalsePositiveRate=numeric(),
        cutoff=numeric(),
        stringsAsFactors=FALSE
)

for(i in seqList){
  for(algo in colnames(predScore)[4:ncol(predScore)]){
    if(algo!="LinearSVM"){
      curTable <- predScore[,1:2]
      curTable$score <- ifelse(predScore[, algo]>=i,1,0)
      cm <- table(curTable$score, curTable$label)
      if(!("0" %in% rownames(cm))){cm <- rbind(cm, "0"=cbind(0,0)); rownames(cm)[2]='0'}
      if(!("1" %in% rownames(cm))){cm <- rbind(cm, "1"=cbind(0,0)); rownames(cm)[2]='1'}
      FPR <- cm["1","nonSRG"]/(cm["1","nonSRG"]+cm["0","nonSRG"])
      TPR <- cm["1","SRG"]/(cm["0","SRG"]+cm["1","SRG"])
      prcTable <- prcTable %>%
        add_row(algo=algo, FalsePositiveRate=FPR, TruePositiveRate=TPR, cutoff=i)
    }
  }
}

prcTable <- prcTable %>%
  dplyr::mutate(algo_short=apply(., 1, function(x){algo_name[x[2]]})) %>%
  dplyr::mutate(algo_short=factor(algo_short, levels=algoLevel))
p2 <- ggplot(prcTable) + 
  geom_line(aes(x=FalsePositiveRate, y=TruePositiveRate, color=algo), alpha=0.8) +
  scale_color_brewer(palette="Dark2") +
  geom_abline(intercept = 0, color="#A0A0A0", linetype = "dashed") +
  theme_bw()
p2
ggsave(filename = "AUC.pdf", p2, width = 4.5, height=2)

aucList <- c()
for(i in colnames(predScore)[4:ncol(predScore)]){
  if(i=="LinearSVM"){
    aucList <- c(aucList, NA)
  }else{
    aucList <- c(aucList, round(auc(predScore$label, predScore[,i])[1], 4))
  }
}
names(aucList) <- colnames(predScore)[4:ncol(predScore)]
aucList <- aucList[order(aucList, decreasing = TRUE)]
aucList
