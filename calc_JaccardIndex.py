## This script contain functions to calculate jaccard index given a gene set annotation file (.gmt) given a set of seed genes.

import os
import math
import re
import yaml
import scipy.stats as sc
import numpy as np
import pandas as pd
import function as SGfn

dirname = os.path.dirname(__file__)
#dirname = "./"

class Gene_JI(SGfn.Gene):
    def __init__(self, name, label):
        super().__init__(name, label)
        self.JIscore = 0
        self.count = 0

    def addJIscore(self, score):
        self.JIscore += score

    def addCount(self):
        self.count += 1

    def resetScore(self):
        self.JIscore = 0
        self.count = 0

def createGeneList(geneList):
    GeneInfo = {}
    for gene in geneList:
        GeneInfo[gene] = Gene_JI(gene, 0)
    return GeneInfo

def labelSeedGenes(GeneInfo, seedGeneList):
    for g in seedGeneList:
        try:
            GeneInfo[g].markGene()
        except KeyError:
            pass

    return GeneInfo

def calcJIforTerms(GenesInfo, seedGeneList, inFileName, geneAliasDict):
    jaccard_index = {}
    nGeneSet = len(GenesInfo)

    with open(inFileName, encoding="utf8", errors='ignore') as fp:
        for line in fp:
            info = line.strip().split("\t")
            term = info[0]
            tGene = []
            for i in range(2, len(info)):
                try:
                    if geneAliasDict[info[i]] in GenesInfo:
                        tGene.append(geneAliasDict[info[i]])
                except KeyError:
                    if info[i] in GenesInfo:
                        tGene.append(info[i])

            n_overlap = len(list(set(seedGeneList) & set(tGene)))
            n_all = len(list(set(seedGeneList + tGene)))
            n_sleep_only = len(seedGeneList) - n_overlap
            n_term_only = len(tGene) - n_overlap
            odd_ratio, pvalue = sc.fisher_exact([[n_overlap, n_sleep_only], [n_term_only, (nGeneSet-n_all)]])
            jaccard_index[term] = {"JI": (n_overlap / n_all), "nGene": n_all, "nOverlap": n_overlap, "FISHER.OR": odd_ratio, "FISHER.Pval": pvalue}
    return jaccard_index

def converttoVector(geneList, geneOrder):
    vList = []
    for gene in geneOrder:
        try:
            vList.append(math.log2(geneList[gene].JIscore + 1))
        except KeyError:
            vList.append(np.nan)
    return vList

def resetGenescore(GeneInfo):
    for gene in GeneInfo.keys():
        GeneInfo[gene].resetScore()

def writeJIterms(jaccard_index, outFilePrefix):
    wf = open(outFilePrefix + "_jaccardIndex.txt", "w")
    wf.write("term\tjaccard_index\tnGene\tnOverlap\tFisher.Pval\tFisher.OR\n")
    for term in jaccard_index.keys():
        wf.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (term,
                                               jaccard_index[term]["JI"],
                                               jaccard_index[term]["nGene"],
                                               jaccard_index[term]["nOverlap"],
                                               jaccard_index[term]["FISHER.Pval"],
                                               jaccard_index[term]["FISHER.OR"]))
    wf.close()

def calcJIforGene(geneAliasDict, GenesInfo, jaccard_index, inFileName):
    with open(inFileName, encoding="utf8", errors='ignore') as fp:
        for line in fp:
            info = line.strip().split("\t")
            term = info[0]
            for i in range(2, len(info)):
                if jaccard_index[term]["FISHER.Pval"] <= 0.05:
                    try:
                        GenesInfo[geneAliasDict[info[i]]].addJIscore(jaccard_index[term]["JI"])
                        GenesInfo[geneAliasDict[info[i]]].addCount()
                    except KeyError:
                        pass

def calcJIforSimTerm(geneAliasDict, inFileName, simCutoff, ovlpcorfCutoff):
    outfn = re.sub('.[a-z]\w$', 'yaml', inFileName)

    with open(inFileName, encoding="utf8", errors='ignore') as fp:
        gmt_dict = {}
        for line in fp:
            info = line.strip().split("\t")
            term = info[0]
            gmt_dict[term] = []
            for i in range(2, len(info)):
                try:
                    gmt_dict[term].append(geneAliasDict[info[i]])
                except KeyError:
                    gmt_dict[term].append(info[i])

    termList = list(gmt_dict.keys())

    ji_coef = {}
    for i in range(len(termList)):
        for j in range(len(termList)):
            if i == j:
                continue

            g1 = gmt_dict[termList[i]]
            g2 = gmt_dict[termList[j]]
            n_overlap = len(list(set(g1) & set(g2)))
            n_all = len(list(set(g1 + g2)))
            nt1 = len(g1)
            nt2 = len(g2)
            n_min = min(nt1, nt2)

            if n_overlap / n_all > simCutoff or n_overlap / n_min > ovlpcorfCutoff:
                try:
                    ji_coef[termList[i]].append(termList[j])
                except KeyError:
                    ji_coef[termList[i]] = [termList[j]]

    with open(outfn, 'w') as file:
        documents = yaml.dump(ji_coef, file)

def calcJIforGeneAverageSimilarTerms(geneAliasDict, GenesInfo, jaccard_index, simCutoff, ovlpcorfCutoff, inFileName):
    yamlFN = re.sub('.[a-z]\w$', 'yaml', inFileName)
    if os.path.exists(yamlFN):
        with open(yamlFN) as file:
            simJIterm = yaml.load(file, Loader=yaml.FullLoader)
    else:
        calcJIforSimTerm(geneAliasDict, inFileName, simCutoff, ovlpcorfCutoff)
        with open(yamlFN) as file:
            simJIterm = yaml.load(file, Loader=yaml.FullLoader)

    with open(inFileName, encoding="utf8", errors='ignore') as fp:
        for line in fp:
            info = line.strip().split("\t")
            term = info[0]
            for i in range(2, len(info)):
                if jaccard_index[term]["FISHER.Pval"] <= 0.05:
                    try:
                        simTerm = simJIterm[term] + [term]
                        maxScore = 0
                        nTerm = 0
                        for tt in simTerm:
                            if jaccard_index[tt]["FISHER.Pval"] <= 0.05:
                                nTerm += 1
                                maxScore = max(jaccard_index[tt]["JI"], maxScore)
                        try:
                            GenesInfo[geneAliasDict[info[i]]].addJIscore(maxScore/nTerm)
                            GenesInfo[geneAliasDict[info[i]]].addCount()
                        except KeyError:
                            pass
                    except KeyError:
                        try:
                            GenesInfo[geneAliasDict[info[i]]].addJIscore(jaccard_index[term]["JI"])
                            GenesInfo[geneAliasDict[info[i]]].addCount()
                        except KeyError:
                            pass

def writeJIforGene(GenesInfo, outFilePrefix):
    wf = open(outFilePrefix + "_JIscore.txt", "w")
    wf.write("GeneSymbol\tJIscore\tnTerms\n")
    for gene in GenesInfo.keys():
        wf.write("%s\t%s\t%s\n" % (gene, math.log2(GenesInfo[gene].JIscore+1), GenesInfo[gene].count))
    wf.close()

def process_database_gene(rName, inFileName, geneAliasDict, GenesInfo, averagesimilarterm, simCutoff, ovlpcorfCutoff, seedGeneList, wkdir):
    outFilePrefix = wkdir + "Database_JI/" + rName
    jaccard_index = calcJIforTerms(GenesInfo, seedGeneList, inFileName, geneAliasDict)
    writeJIterms(jaccard_index, outFilePrefix)
    if averagesimilarterm:
        calcJIforGeneAverageSimilarTerms(geneAliasDict, GenesInfo, jaccard_index, simCutoff, ovlpcorfCutoff, inFileName)
    else:
        calcJIforGene(geneAliasDict, GenesInfo, jaccard_index, inFileName)
    writeJIforGene(GenesInfo, outFilePrefix)
    return GenesInfo

def convertJITermstoGene(geneAliasDict, GenesInfo, jaccard_index):
    for term in jaccard_index.keys():
        try:
            GenesInfo[geneAliasDict[term]].addJIscore(jaccard_index[term]["JI"])
            GenesInfo[geneAliasDict[term]].addCount()
        except KeyError:
            pass

def process_database_term(rName, inFileName, geneAliasDict, GenesInfo, seedGeneList, wkdir):
    outFilePrefix = wkdir + "Database_JI/" + rName
    jaccard_index = calcJIforTerms(GenesInfo, seedGeneList, inFileName, geneAliasDict)
    writeJIterms(jaccard_index, outFilePrefix)
    convertJITermstoGene(geneAliasDict, GenesInfo, jaccard_index)
    return GenesInfo

def process_input(geneSetDataList, GenesInfo, geneAliasDict, seedGeneList, wkdir, suffix, inserttodataframe=True, averagesimilarterm=True, simCutoff=0.5, ovlpcorfCutoff=0.9, geneOrder=[], dataset=False):
    with open(geneSetDataList, "r") as scFn:
        for line in scFn:
            info = line.strip().split("\t")
            rName = info[0]
            fn = info[1]
            dataType = info[2]
            print("Calculating JI for %s" % rName)
            if dataType == 'geneJI':    # use gene level JI as features
                GenesInfo = process_database_gene(rName, fn, geneAliasDict, GenesInfo, averagesimilarterm, simCutoff, ovlpcorfCutoff, seedGeneList, wkdir+suffix)
            elif dataType == 'termJI': # use term level JI as feature, only work if each term is a gene (e.g. protein-protein interactions)
                GenesInfo = process_database_term(rName, fn, geneAliasDict, GenesInfo, seedGeneList, wkdir+suffix)
            if inserttodataframe and len(geneOrder) > 0:
                infoList = converttoVector(GenesInfo, geneOrder)
                dataset.insert(1, rName, infoList, True)
            resetGenescore(GenesInfo)
    return dataset

### Calculate JI score from gene set collections based on a seed gene list and insert this as additional feature to the dataset
def main(dataset, geneSetDataList, seedGeneList, SUFFIX, outPATH=dirname, inserttodataframe=True, averagesimilarterm=True, simCutoff=0.5, ovlpcorfCutoff=0.9):
    ## Create folder to store JI result
    if not os.path.exists(outPATH + SUFFIX + "Database_JI"):
        os.makedirs(outPATH + SUFFIX + "Database_JI")

    ## Create a dictionary to map gene names alias
    geneAliasDict = SGfn.createHumanGeneAliasDict()

    ## Create a class to store JI info for each gene in the dataset
    geneList = list(dataset["GeneSymbol"])
    allGeneInfo = createGeneList(geneList)
    allGeneInfo = labelSeedGenes(allGeneInfo, seedGeneList)

    dataset = process_input(geneSetDataList=geneSetDataList, GenesInfo=allGeneInfo, geneAliasDict=geneAliasDict, seedGeneList=seedGeneList,
                            wkdir=outPATH, suffix=SUFFIX, inserttodataframe=inserttodataframe, averagesimilarterm=averagesimilarterm, simCutoff=simCutoff, ovlpcorfCutoff=ovlpcorfCutoff, geneOrder=geneList, dataset=dataset)

    dataset.to_csv(outPATH + SUFFIX + "_withJI.csv", index=False)

    return dataset

## test run
#dataset = pd.read_csv("./DATA/GeneList.txt")
#seedGeneList = SGfn.readSeedGeneList("./DATA/testseedGeneList.txt")
#dataset = main(dataset=dataset, geneSetDataList="./DATA/gmt_file/geneset2.list", seedGeneList=seedGeneList, SUFFIX='testRun_20210831_averagesimTerm.5JI', outPATH="./OUT/")

