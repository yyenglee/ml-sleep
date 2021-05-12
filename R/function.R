library(homologene)
library(limma)
library(reshape)
library(ggplot2)
library(gridExtra)

## convert to human homolog
convertModelHs <- function(gene, inTax){
  homoD <- homologene(gene, inTax = inTax, outTax = 9606, db = homologeneData2) %>%
    dplyr::select(c(1,2)) %>%
    magrittr::set_colnames(c("mGene", "hGene"))
  return(homoD)
} 

## convert alias to official gene symbol
GeneAliasConversion <- function(gene, REF){
  return(alias2SymbolUsingNCBI(gene,REF)$Symbol)
}

## Read gene table, convert gene alias to gene name and look for human homolog gene
createSleepGeneList <- function(sleepGeneFile,excludemodel=c()){
  sleepGene = readxl::read_excel(sleepGeneFile, sheet = 1)  %>%
    dplyr::filter(Tag == "1", !(Tier %in% excludemodel))  %>%
    dplyr::select(Gene, Model_for_genename, Phenotype_curated, HGNC, Tier)
  sleepGene$officialSymbol <- apply(sleepGene,1,function(x){GeneAliasConversion(gene=as.character(x[1]),REF=geneAliasDict[[x[2]]])})
  geneDict <- data.frame()
  for(sp in unique(sleepGene$Model_for_genename)){
    modelTable <- subset(sleepGene,sleepGene$Model_for_genename==sp)
    geneList <- as.character(modelTable$officialSymbol)
    #geneList <- GeneAliasConversion(gene=geneList,REF=geneAliasDict[[sp]])
    if(sp == "human"){
      humanHom <- data.frame(cbind(geneList,geneList))
      colnames(humanHom) <- c("mGene", "hGene")
      humanHom$model <- sp
      humanHom$phenotype <- apply(humanHom,1,function(x){as.character(modelTable[modelTable$officialSymbol==x[1],"Phenotype_curated"])})
    }else{
      humanHom <- convertModelHs(geneList, inTax=taxid[[sp]])
      humanHom$model <- sp
      humanHom$phenotype <- apply(humanHom,1,function(x){as.character(modelTable[modelTable$officialSymbol==x[1],"Phenotype_curated"])})
    }
    geneDict <- rbind(geneDict, humanHom)
  }
  geneDict <- geneDict %>% 
    left_join(sleepGene, by=c("hGene"="HGNC")) %>%
    dplyr::select(mGene, hGene, model,phenotype, Tier)
  return(geneDict)
}

## return log3
log3 <- function(x) log(x)/log(3)

## return 10^x
x_10 <- function(x) 10^x

## return 3^x
x_3 <- function(x) 3^x

## return 2^x
x_2 <- function(x) 2^x

## return count of NA in each row
count_na <- function(x) sum(!is.na(x))

## return standard error of mean
se <- function(x) sqrt(var(x)/length(x))

## replace evi<3 to 0, to reduce noise
filterEF <- function(x) ifelse(x>=3.0,x,0.0)

## keep if at the top 95 percentile
quantile95 <- function(x) ifelse(x>=quantile(x,na.rm=TRUE,c(0.9)),1,0)

## return average EF for top 99 percentile
top99percentileEF <- function(x){
  x <- x[!is.na(x)]
  q99 <- quantile(x, .99)
  return(mean(x[x>=q99]))
}

## rescale data from 0 to 1  
scale2 <- function(x, na.rm = FALSE) (x - min(x)) / (max(x) - min(x))

## return zscore
zscore <- function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
  
## read file and return table with unique gene id and value (maximum value for each gene)
readFN <- function(fn){
  data <- read.table(fn, header=T, sep="\t", quote="") %>%
    dplyr::select(c(1,2)) %>%
    magrittr::set_colnames(c("SYMBOL", "value")) %>%
    ## if multiple gene name for the same row, keep the first one
    dplyr::mutate(GeneSymbol = stringr::str_extract(SYMBOL, "[A-Za-z0-9-]+[^|]"), value=as.numeric(value)) %>%
    tidyr::drop_na() %>%
    dplyr::group_by(GeneSymbol) %>%
    ## if a gene with multiple value, keep the highest
    dplyr::slice(which.max(value)) %>%
    dplyr::select(GeneSymbol, value) %>%
    ungroup()
  return(data)
}

## Filter gene symbol to gene with human homolog and convert gene symbol to HGNC
##  currently only work with mouse, drosophila and human data
##  additional reference needed for other species
convertGeneNameToHGNC <- function(data, model){
  data <- data %>%
    dplyr::mutate(OfficialGeneName=GeneAliasConversion(gene=data$GeneSymbol,REF=geneAliasDict[[model]])) %>%
    dplyr::filter(!is.na(OfficialGeneName)) %>%
    group_by(OfficialGeneName) %>%
    dplyr::slice(which.max(value)) %>%
    ungroup() %>%
    dplyr::select(OfficialGeneName,value)
  if(model=="human"){
    data <- data %>%
      dplyr::rename(GeneSymbol=OfficialGeneName)
  }else{
    HsGeneName <- convertModelHs(data$OfficialGeneName, inTax=taxid[[model]])
    data <- data %>%
      left_join(HsGeneName, by=c("OfficialGeneName"="mGene")) %>%
      dplyr::filter(!is.na(hGene)) %>%
      group_by(hGene) %>%
      dplyr::slice(which.max(value)) %>%
      ungroup() %>%
      dplyr::select(hGene, value) %>%
      dplyr::rename(GeneSymbol=hGene)
  }
  return(data)
}

## label gene as sleep regulating gene (SRG) or otherwise (nonSRG)
## sleepGeneList must be character
labelData <- function(data, sleepGeneList){
  if(typeof(sleepGeneList) != "character"){
    stop(print("'sleepGeneList' must be a list of character"))
  }
  data <- data %>%
    dplyr::mutate(label = ifelse(GeneSymbol %in% sleepGeneList, "SRG","nonSRG")) %>%
    dplyr::select(GeneSymbol,label,everything())
  return(data)
}

## Capped value to .999
setCap <- function(data, dataType, capLimit=.999){
  ## capped the top 0.001% value to reduce effect of outliers
  if(dataType=="numeric"){
    capped <- as.numeric(quantile(data$value,capLimit))
  }else{
    capped <- sort(data$value)[round(length(data$value) *capLimit)]
  }
  data <- dplyr::mutate(data,value=ifelse(value>capped,capped,value))
  return(data)
}
  
## Calculate evidence factor - 2 states
calEF_2state <- function(data, fn){
  SRG_ecdf <- ecdf(data %>% filter(label=="SRG") %>% .$value)
  nonSRG_ecdf <- ecdf(data %>% filter(label=="nonSRG") %>% .$value)
  
  ## reduce dataset into roughly 100 blocks
  steps <- round(-log10((max(data$value)-min(data$value))/100))

  SRG_X <- knots(SRG_ecdf)
  nonSRG_X <- knots(nonSRG_ecdf)
  rangeX <- seq(round(min(SRG_X,nonSRG_X),steps),round(max(SRG_X,nonSRG_X),steps),1/(10^steps))
  #rangeX <- seq(round(min(SRG_X,nonSRG_X),2),round(max(SRG_X,nonSRG_X),2),length.out=100)
  
  SRG_Y <- round(SRG_ecdf(rangeX),10)
  nonSRG_Y <- round(nonSRG_ecdf(rangeX),10)
  
  # Set value=0 for jaccard index data
  if(grepl("jaccardIndex", fn)|grepl("JIscore", fn)){
    SRG_Y[1] <- 0
    nonSRG_Y[1] <- 0
  }
  
  evFc <- as.data.frame(cbind(rangeX,SRG_Y,nonSRG_Y)) %>%
    dplyr::mutate(evidenceFactors=(1-SRG_Y)/(1-nonSRG_Y), scale=as.character(round(rangeX,steps)))
  ## if nonSRG group reach 100%, use last EF value for remaining SRG group
  lastEF <- tail(subset(evFc,evFc$evidenceFactors<Inf&evFc$evidenceFactors>0&is.numeric(evFc$evidenceFactors))$evidenceFactors,1)
  if(length(table(evFc$evidenceFactors))>1){
    evFc$evidenceFactors <- apply(evFc,1,function(x){ifelse((1-as.numeric(x[3])==0)|(1-as.numeric(x[2])==0),as.numeric(lastEF),as.numeric(x[4]))})
  }
  #dplyr::mutate(evidenceFactors=ifelse(nonSRG_Y==1,lag(evidenceFactors),evidenceFactors))

  return(evFc)
}

## Calculate evidence factor - bin
calEF_bin <- function(data, fn, nbin=NULL){
  new_evFc <- data.frame(start=double(),
                         end=double(),
                         SRG_Y=integer(),
                         nonSRG_Y=integer(),
                         stringsAsFactors=FALSE)
  
  steps <- length(unique(data$value))
  if(steps <= 3){return(new_evFc)}
  if(length(nbin)==0){
    nbin <- floor(ifelse(steps<10, floor(steps/3), 100))
  }
  
  SRGList <- data %>% filter(label=="SRG") %>% .$value
  nonSRGList <- data %>% filter(label=="nonSRG") %>% .$value

  ## spread dataset into roughly 100 blocks
  initStep <- floor(min(data$value))
  endStep <- ceiling(max(data$value)) # increase maximum by 1%, so that it will always include the maximum value
  #endStep <- ceiling(max(data$value)*1.01)
  rangeX <- round(seq(initStep, endStep, length.out = nbin + 1), 2)
  rangeX <- sort(unique(c(quantile(data$value, seq(0,1,length.out = 10)), rangeX)))
  #rangeX <-  c(quantile(data$value, seq(0,1,length.out = 101)))
  
  cur_evFc <- as.data.frame(cbind(c(rangeX[1:length(rangeX)-1]),c(rangeX[2:length(rangeX)]))) %>%
    magrittr::set_colnames(c("start","end")) %>%
    rowwise() %>%
    dplyr::mutate(SRG_Y=sum(SRGList>start&SRGList<=end), 
                   nonSRG_Y=sum(nonSRGList>start&nonSRGList<=end))
    # dplyr::mutate(SRG_Y=sum(SRGList>=start&SRGList<end), 
    #               nonSRG_Y=sum(nonSRGList>=start&nonSRGList<end))
  
  for(j in 1:nrow(cur_evFc)){
    if(j==1){
      curStart <- as.numeric(cur_evFc[j, "start"])
      curEnd <- as.numeric(cur_evFc[j, "end"])
      curSRG <- as.numeric(cur_evFc[j, "SRG_Y"])
      curnonSRG <- as.numeric(cur_evFc[j, "nonSRG_Y"])
    }else{
      ## restricted to at least 1% genes in a bin
      if((curnonSRG>=(length(nonSRGList)*.1) & curSRG>=(length(SRGList)*.01)) |
         (curSRG>=(length(SRGList)*.1) & curnonSRG>=(length(nonSRGList)*.01))){
        start <- curStart
        end <- curEnd
        SRG_Y <- curSRG
        nonSRG_Y <- curnonSRG
        curRow <- cbind(start, end, SRG_Y, nonSRG_Y)
        new_evFc <- rbind(new_evFc, curRow)
        
        ## reset tick
        curStart <- as.numeric(cur_evFc[j, "start"])
        curEnd <- as.numeric(cur_evFc[j, "end"])
        curSRG <- as.numeric(cur_evFc[j, "SRG_Y"])
        curnonSRG <- as.numeric(cur_evFc[j, "nonSRG_Y"])
        
      }else{
        curSRG <- curSRG + as.numeric(cur_evFc[j, "SRG_Y"])
        curnonSRG <- curnonSRG + as.numeric(cur_evFc[j, "nonSRG_Y"])
        curEnd <- as.numeric(cur_evFc[j, "end"])
      }
    }

  }
  if((curnonSRG>=(length(nonSRGList)*.1) & curSRG>=(length(SRGList)*.01)) |
  (curSRG>=(length(SRGList)*.1) & curnonSRG>=(length(nonSRGList)*.01)) |
  nrow(new_evFc)==0){
  #if((curnonSRG>=(length(nonSRGList)*.05) & curSRG>=(length(SRGList)*.01) & curSRG > 2) |
  # (curSRG>=(length(SRGList)*.05) & curSRG > 2 & (curnonSRG>=(length(nonSRGList)*.01)|curnonSRG>=100)) |
  #  nrow(new_evFc)==0){
    start <- curStart
    end <- curEnd
    ## if the last bin has only sleep genes,
    ##  replace non-sleep gene count to 1, to prevent output of inf
    if(curnonSRG == 0){
      nonSRG_Y = 1
    }else{
      nonSRG_Y <- curnonSRG
    }
    SRG_Y <- curSRG
    
    curRow <- cbind(start, end, SRG_Y, nonSRG_Y)
    new_evFc <- rbind(new_evFc, curRow)
    
  }else{
    new_evFc[nrow(new_evFc),"end"] <- curEnd
    new_evFc[nrow(new_evFc),"SRG_Y"] <- new_evFc[nrow(new_evFc),"SRG_Y"] + curSRG
    new_evFc[nrow(new_evFc),"nonSRG_Y"] <- new_evFc[nrow(new_evFc),"nonSRG_Y"] + curnonSRG
  }

  evFc <- new_evFc %>%
    dplyr::mutate(rangeX = (start+end)/2,
                  SRG_Y=SRG_Y/sum(SRG_Y), 
                  nonSRG_Y=nonSRG_Y/sum(nonSRG_Y)) %>%
    dplyr::mutate(evidenceFactors=SRG_Y/nonSRG_Y) %>%
    dplyr::select(rangeX, start, end, evidenceFactors, SRG_Y, nonSRG_Y) %>%
    drop_na()
  return(evFc)
}

## Match evidence factor for each gene
matchEF_bin <- function(data, evFc){
  data$evidenceFactors <- NA
  for(k in 1:nrow(evFc)){
    curStart <- evFc[k,"start"]
    curEnd <- evFc[k, "end"]
    data[data$value>=curStart&data$value<=curEnd, "evidenceFactors"] <- evFc[k,"evidenceFactors"]
  }
  data <- data %>%
    dplyr::select(GeneSymbol, evidenceFactors) %>%
    arrange(desc(evidenceFactors))
  return(data)
}

## Match evidence factor for each gene
matchEF <- function(data, evFc){
  steps <- round(-log10((max(data$value)-min(data$value))/100))
  data <- data %>%
    dplyr::mutate(scale = as.character(round(value,steps))) %>%
    left_join(evFc, by = "scale") %>%
    dplyr::select(GeneSymbol, evidenceFactors) %>%
    arrange(desc(evidenceFactors))
  return(data)
}

## Process a dataset and return table with evidence factor for that file
evalPerDataEF <- function(lineInfo, labelGene=sleepGeneList$hGene){
  source <- lineInfo["Source"]
  ID <- as.character(lineInfo["ID"])
  fn <- paste0(pdir,as.character(lineInfo["Input"]))
  model <- as.character(lineInfo["Models"])
  dataType <- as.character(lineInfo["DataType"])
  dat <- readFN(fn)
  dat <- convertGeneNameToHGNC(dat, model)
  dat <- labelData(dat, sleepGeneList=labelGene)
  dat <- setCap(dat, dataType)
  ef <- calEF(dat,dataType)
  datEF <- matchEF(dat,ef)
  datEF <- labelData(datEF, sleepGeneList=labelGene)
  return(datEF)
}

## Prepare plot
drawDensityPlot <- function(data, nCol, name=NULL){
  subDat <- data %>%
    dplyr::select(c(1,2,nCol)) %>%
    magrittr::set_colnames(c("GeneSymbol","label","value"))
    
  p <- ggplot(subDat,aes(x=value)) + 
    geom_histogram(aes(y=..density.., fill=label), position='identity', alpha=0.4, bins=50) +
    geom_density(aes(color=label),size=0.8,alpha=0.5) +
    labs(title=paste0(name,"\n#gene=",dim(data)[1]," | #srg=",as.numeric(table(data$label)[2]))) +
    scale_color_manual(values=c('dodgerblue','coral')) +
    scale_fill_manual(values=c('dodgerblue','coral')) + 
    theme_bw() + 
    theme(legend.position="none")
  return(p)
}

drawCumDict <- function(ef, name=NULL){
  tmpDat <- melt(ef[,c("rangeX", "SRG_Y", "nonSRG_Y", "evidenceFactors")], id=c("rangeX","evidenceFactors"))
  p <- ggplot(tmpDat, aes(x=rangeX, y=value, color=variable)) +
    geom_point(size=0.5) +
    scale_color_manual(values=c("coral","dodgerblue"),labels=c("SRG","nonSRG")) +
    labs(title=name) +
    theme_bw() + 
    guides(colour = guide_legend(override.aes = list(size=3))) +
    theme(legend.position="top",
          legend.title=element_blank())
  return(p)
}

drawEFPlot <- function(ef, name=NULL){
  # p <- ggplot(ef, aes(x=rangeX, y=evidenceFactors)) +
  #   geom_line(size=0.5) + 
  #   labs(title=name) +
  #   theme_bw()
  p <- ggplot(ef) +
    geom_segment(aes(x=start, xend=end, y=evidenceFactors, yend=evidenceFactors), size=1.2, color="#23395D") +
    labs(title=name) +
    theme_bw()
  return(p)
}

## Add raw value to table
addRawValue <- function(PATH, rawValuetable, inputsleepGeneList){
  subPATH <- PATH[!PATH$ID%in%colnames(rawValuetable),]
  print(paste(nrow(subPATH),"data metrics to be added"))
  print(head(subPATH))
  for(i in 1:nrow(subPATH)){
    source <- subPATH[i,"Source"]
    ID <- subPATH[i,"ID"]
    fn <- paste0(pdir,subPATH[i,"Input"])
    model <- subPATH[i,"Models"]
    dataType <- subPATH[i,"DataType"]
    dat <- readFN(fn)
    dat[dat$value==Inf,"value"] <- max(dat[dat$value<Inf,"value"])+0.1
    dat <- convertGeneNameToHGNC(dat, model)
    dat <- labelData(dat, sleepGeneList=inputsleepGeneList)

    rawValuetable <- rawValuetable %>%
      left_join(dat[,c(1,3)], by=c("GeneSymbol")) %>%
      dplyr::rename(!!ID := value)
  }
  return(rawValuetable)
}

