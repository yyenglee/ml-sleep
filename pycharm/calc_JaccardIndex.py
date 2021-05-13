import argparse
import pandas as pd
import numpy as np
import function as SRGfn
import math
import scipy.stats as sc
import os


class Gene_JI(SRGfn.Gene):
    def __init__(self, name, label):
        super().__init__(name, label)
        self.JIscore = 0
        self.count = 0

    def addJIscore(self, score):
        self.JIscore += score

    def maxJIscore(self, score):
        self.JIscore = max(self.JIscore, score)

    def addCount(self):
        self.count += 1

    def cleanScore(self):
        self.JIscore = 0
        self.count = 0


def resetGenescore(geneList):
    for gene in geneList.keys():
        geneList[gene].cleanScore()


def calcJIforTerms(geneDict, geneList, SRG_list, inFileName):
    jaccard_index = {}
    data = open(inFileName, "r")
    nGeneSet = len(geneList)
    for line in data:
        info = line.strip().split("\t")
        term = info[0]
        tGene = []
        for i in range(2, len(info)):
            try:
                if geneDict[info[i]] in geneList:
                    tGene.append(geneDict[info[i]])
            except KeyError:
                if info[i] in geneList:
                    tGene.append(info[i])

        n_overlap = len(list(set(SRG_list) & set(tGene)))
        n_all = len(list(set(SRG_list + tGene)))
        n_sleep_only = len(SRG_list) - n_overlap
        n_term_only = len(tGene) - n_overlap
        odd_ratio, pvalue = sc.fisher_exact([[n_overlap, n_sleep_only], [n_term_only, (nGeneSet-n_all)]])
        jaccard_index[term] = {"JI": (n_overlap / n_all), "nGene": n_all, "nOverlap": n_overlap, "FISHER.OR": odd_ratio, "FISHER.Pval": pvalue}
    data.close()
    return jaccard_index


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


def calcJIforGene(geneDict, geneList, jaccard_index, inFileName):
    data = open(inFileName, "r")
    for line in data:
        info = line.strip().split("\t")
        term = info[0]
        for i in range(2, len(info)):
            if jaccard_index[term]["FISHER.Pval"] <= 0.05:
                try:
                    geneList[geneDict[info[i]]].addJIscore(jaccard_index[term]["JI"])
                    geneList[geneDict[info[i]]].addCount()
                except KeyError:
                    pass


def convertJITermstoGene(geneDict, geneList, jaccard_index):
    for term in jaccard_index.keys():
        try:
            geneList[geneDict[term]].addJIscore(jaccard_index[term]["JI"])
            geneList[geneDict[term]].addCount()
        except KeyError:
            pass

def writeJIforGene(geneList, outFilePrefix):
    wf = open(outFilePrefix + "_JIscore.txt", "w")
    wf.write("GeneSymbol\tJIscore\tnTerms\n")
    for gene in geneList.keys():
        wf.write("%s\t%s\t%s\n" % (gene, math.log2(geneList[gene].JIscore+1), geneList[gene].count))
    wf.close()


def converttoVector(geneList, geneOrder):
    vList = []
    for gene in geneOrder:
        try:
            vList.append(math.log2(geneList[gene].JIscore + 1))
        except KeyError:
            vList.append(np.nan)
    return vList


def createGeneList(geneOrder):
    geneList = {}
    for gene in geneOrder:
        geneList[gene] = Gene_JI(gene, 0)
    return geneList


def process_database_gene(rName, fn, geneDict, geneList, SRG_list, DATAPATH, curwkPATH):
    inFileName = DATAPATH + fn
    outFilePrefix = curwkPATH + "term_JI/" + rName
    jaccard_index = calcJIforTerms(geneDict, geneList, SRG_list, inFileName)
    writeJIterms(jaccard_index, outFilePrefix)
    calcJIforGene(geneDict, geneList, jaccard_index, inFileName)
    writeJIforGene(geneList, outFilePrefix)
    return geneList


def process_database_term(rName, fn, geneDict, geneList, SRG_list, DATAPATH, curwkPATH):
    inFileName = DATAPATH + fn
    outFilePrefix = curwkPATH + "term_JI/" + rName
    jaccard_index = calcJIforTerms(geneDict, geneList, SRG_list, inFileName)
    writeJIterms(jaccard_index, outFilePrefix)
    convertJITermstoGene(geneDict, geneList, jaccard_index)
    return geneList


def process_input(geneSetList, geneList, geneDict, SRG_list, DATAPATH, curwkPATH, inserttodataframe=True, geneOrder=[], dataset=False):
    source = open(geneSetList, "r")
    for line in source:
        info = line.strip().split("\t")
        rName = info[0]
        fn = info[1]
        dataType = info[2]
        print("Calculating JI for %s" % rName)
        if dataType == 'geneJI':    # use gene level JI as features
            geneList = process_database_gene(rName, fn, geneDict, geneList, SRG_list, DATAPATH, curwkPATH)
        elif dataType == 'termJI': # use term level JI as feature, only work if each term is a gene (e.g. protein-protein interactions)
            geneList = process_database_term(rName, fn, geneDict, geneList, SRG_list, DATAPATH, curwkPATH)
        if inserttodataframe and len(geneOrder) > 0:
            infoList = converttoVector(geneList, geneOrder)
            dataset.insert(2, rName, infoList, True)
        resetGenescore(geneList)
    if inserttodataframe:
        return dataset

### actual code to insert JI score from database to the dataset
def runJI(dataset, curwkPATH, geneSetList, SUFFIX, DATAPATH):
    geneDict = SRGfn.createHumanGeneAliasDict()

    if not os.path.exists(curwkPATH + "term_JI"):
        os.makedirs(curwkPATH + "term_JI")

    y = dataset.iloc[:, 1].values
    geneOrder = list(dataset["GeneSymbol"])

    SRG_list = list(dataset["GeneSymbol"][y=='SRG'])
    print(len(SRG_list), SRG_list)

    ## Create a dictionary to store gene list
    geneList = createGeneList(geneOrder)
    dataset = process_input(geneSetList=geneSetList, geneList=geneList, geneDict=geneDict, SRG_list=SRG_list, DATAPATH=DATAPATH, curwkPATH=curwkPATH, inserttodataframe=True, geneOrder=geneOrder, dataset=dataset)

    dataset.to_csv(curwkPATH + "fList_cleanRawValue_withJI_" + SUFFIX + ".csv", index=False)

    return dataset

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('-inTable', type=str, required=True, help="input table with at least 2 columns: genes and labels")
    parser.add_argument('-wkdir', type=str, default="../OUT/")
    parser.add_argument('-gmtList', default="../DATA/gmt_file/geneset.list")
    parser.add_argument('-SUFFIX', type=str, default="")
    parser.add_argument('-DATA_FOLDER', type=str, default="../DATA/")

    args = parser.parse_args()

    dataset = pd.read_csv(args.inTable)
    runJI(dataset=dataset, curwkPATH=args.wkdir, geneSetList=args.gmtList, SUFFIX=args.SUFFIX, DATAPATH=args.DATA_FOLDER)

if __name__ == "__main__":
    main()
