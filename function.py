# Helper libraries

import gzip
import pandas as pd
import numpy as np
import os
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler

dirname = os.path.dirname(__file__)

class Gene:
    def __init__(self, name, label):
        self.name = name
        self.label = label

    def markGene(self):
        self.label = 1

    def unmarkGene(self):
        self.label = 0

def createHumanGeneAliasDict():
    # Create dictionary to map all gene alias to official gene name
    humanREF = dirname + '/REFERENCE/Homo_sapiens.gene_info.gz'
    geneDict = {}
    with gzip.open(humanREF, "rb") as f:
        next(f)
        for line in f:
            line = str(line, 'utf-8')
            info = line.strip().split("\t")
            geneSymbol = info[2]
            geneDict[geneSymbol] = geneSymbol
            geneAlias = info[4].split("|")
            ## if an alias mapped to multiple official gene symbol, the one with smallest entrez id is chosen
            if geneAlias != ['-']:
                for gene in geneAlias:
                    try:
                        geneDict[gene]  ## test if gene already in dictionary
                    except KeyError:
                        geneDict[gene] = geneSymbol
    return geneDict

def readSeedGeneList(filename):
    ## read seed genes from file, line started with # will be skipped, one gene per line
    seedGeneList = []
    with open(filename) as f:
        for line in f:
            if line[0] == "#":
                continue
            else:
                seedGeneList.append(line.strip())
    return seedGeneList

def createSeedGeneBoolean(dataset, seedGeneList):
    ## return numpy array with 1 * number of genes, with values 1 indicate seed genes, 0 if not.
    geneAliasDict = createHumanGeneAliasDict()
    geneList = list(dataset["GeneSymbol"])
    stdseedGeneList = [geneAliasDict[gene] for gene in seedGeneList]
    for gene in stdseedGeneList:
        if gene not in geneList:
            print('%s removed from seedGeneList'%gene)
            stdseedGeneList.remove(gene)
    sgIndex = [geneList.index(sg) for sg in stdseedGeneList]
    sgIndex.sort()

    seedGeneBool = dataset[["GeneSymbol"]].copy()
    seedGeneBool["Y"] = 0
    seedGeneBool.iloc[sgIndex,1] = 1
    y = seedGeneBool.iloc[:,1].values

    return y

def labelSeedGene(indata, seedGeneList):
    ## return dataframe with additional column, 'label' with 1 as seed gene, 0 as non seed gene.
    stdseedGeneList = [1 if gene in seedGeneList else 0 for gene in indata.loc[:,"GeneSymbol"]]
    indata.insert(1, "labels", stdseedGeneList)
    return indata

def extract_info(dataset, seedGeneList):
    ## extract features and labels from dataset
    labeldata = dataset.copy()
    print(labeldata.shape)

    X = labeldata.iloc[:, 1:].values
    y = createSeedGeneBoolean(labeldata, seedGeneList)

    return X, y

def preprocess_feature(X, imputeData=True, scaleFeature=True):
    ## run imputation and scale feature
    if imputeData:
    ## impute missing value
        imputer = SimpleImputer(missing_values=np.nan, strategy="mean")
        imputer = imputer.fit(X)
        X = imputer.transform(X)

    if scaleFeature:
    # scale features
        X = StandardScaler().fit_transform(X)
    return X

def filterTable(dataset, featureList, rowCutOff, exemptGeneList=np.array([])):
    ## select column based on feature list and filter row if proportions of missing value >rowCutOff
    tmpDataset = dataset.loc[:,featureList]
    newDataset = dataset.loc[
        (tmpDataset.isnull().sum(axis=1) < (tmpDataset.shape[1] * rowCutOff)) | (dataset.GeneSymbol.isin(exemptGeneList))]
    fcList = ["GeneSymbol"] + featureList
    newDataset = newDataset.loc[:, fcList]
    return newDataset

def filterbyCor(dataset, maxEFtable, corCutoff=0.8):
    ## run Pearson correlation and return a dataset removing the feature with lower maxEF for feature pairs with correlation > corCutoff (default=0.8)
    if isinstance(maxEFtable.index[0], int):
        maxEFtable.index = maxEFtable.iloc[:, 0]
    corTable = dataset.corr().abs()
    scorTable = corTable.unstack()
    so = scorTable.sort_values(kind="quicksort")
    so = so[so>corCutoff]
    newso = [[s[0],s[1],maxEFtable.loc[s[0],'maxEF'],maxEFtable.loc[s[1],'maxEF']] for s in so.index if s[0]!=s[1]]
    filterlist = list(set([s[0] if s[3]>s[2] else s[1] for s in newso]))
    return dataset.drop(columns=filterlist)

def calcFinalPredScore(finalPredTable, seedGeneList, bgnoisethreshold=0.1, highCthreshold=0.5, predRatioStep=[0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]):
    # calculate and return final prediction score
    ncol = finalPredTable.shape[1]
    ## replace all value less than background noise threshold (default=0.1) with 0
    finalPredTable.iloc[:, 1:] = finalPredTable.iloc[:, 1:].apply(lambda x: [y if (y >= bgnoisethreshold) else None for y in x])
    ## calculate average prediction score for each gene
    RFscore = finalPredTable.iloc[:, 1:].mean(axis=1, numeric_only=True)
    finalPredTable['RFscore'] = RFscore.apply(lambda x: x if x>bgnoisethreshold else 0)

    ## Identify the lowest prediction ratio that lead to a score higher than the bgnoisethreshold
    tmpPredTable = finalPredTable.copy()
    for i in range(1,ncol):
        curW = (predRatioStep[i-1])*10
        tmpPredTable.iloc[:,i] = tmpPredTable.iloc[:,i].apply(lambda x: curW if x>=bgnoisethreshold else None)
    maxRatio = tmpPredTable.iloc[:, 1:ncol].idxmax(axis=1)
    finalPredTable['maxRatio'] = [None if pd.isna(c) else finalPredTable.columns.get_loc(c) for c in maxRatio]
    ## remove genes with all negative prediction results
    finalPredTable = finalPredTable.dropna(thresh=3)

    ## sum up prediction score, weighted by minimum training ratio, and add in additional score if the average prediction score is higher than the high confidence threshold (default=0.5)
    #finalRFsum = finalPredTable.apply(lambda row: row['RFscore']+row['maxRatio'] if row['RFscore']<highCthreshold else row['RFscore']+row['maxRatio']+1, axis=1)
    finalRFsum = finalPredTable.apply(lambda row: row['RFscore'] + row['maxRatio'], axis=1)
    finalPredTable = pd.concat([finalPredTable, finalRFsum], axis=1)
    finalPredTable = finalPredTable.sort_values(0, ascending=False)
    finalPredTable = labelSeedGene(finalPredTable, seedGeneList)
    finalPredTable["sensitivity"] = finalPredTable.loc[:,"labels"].cumsum()/len(seedGeneList)

    return finalPredTable
