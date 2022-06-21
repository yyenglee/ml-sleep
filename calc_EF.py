## This script contain functions to calculate evidence factors (enrichment score) given a set of seed genes.

import os
import math
import time
import numpy as np
import pandas as pd
import function as SGfn

dirname = os.path.dirname(__file__)

def getIndicesofTwoSet(dataset, seedGeneList):
    ## return row indices for seed genes and remaining genes in the dataset
    geneAliasDict = SGfn.createHumanGeneAliasDict()
    nrow, ncol = dataset.shape
    geneList = list(dataset["GeneSymbol"])
    stdseedGeneList = [geneAliasDict[gene] for gene in seedGeneList]

    sgIndex = [geneList.index(sg) for sg in stdseedGeneList]
    sgIndex.sort()
    nonsgIndex = []
    for i in range(0,nrow):
        if i not in sgIndex:
            nonsgIndex.append(i)

    return sgIndex, nonsgIndex

def defineEqualBinMarker(FCvalue, fixEqualbinSize = 100):
    # create markers based on fixed data range
    if len(set(FCvalue)) < 10:
        binSize = int(round(len(set(FCvalue)) / 3, 0))
    else:
        binSize = fixEqualbinSize
    return np.linspace(math.floor(FCvalue.min()), math.ceil(FCvalue.max()), binSize+1)

def defineQuantileBinMarker(FCvalue, fixQuantilebinSize = 10):
    # create markers based on value range
    return np.array([np.quantile(FCvalue, i) for i in np.linspace(0, 1, fixQuantilebinSize)])

def mergeBinMarker(range1, range2):
    # merge and return unique markers from two lists
    rangeX = np.concatenate((range1, range2))
    rangeX = np.unique(rangeX)
    rangeX.sort()
    return rangeX

def defineBinMarkers(FCvalue, EqualBins=True, QuantileBins=True, fixEqualbinSize = 100, fixQuantilebinSize = 10, ebinround=True, digit=2):
    if EqualBins and QuantileBins:
        range1 = defineEqualBinMarker(FCvalue, fixEqualbinSize)
        if ebinround:
            range1 = np.round(range1, digit)
        range2 = defineQuantileBinMarker(FCvalue, fixQuantilebinSize)
        rangeX = mergeBinMarker(range1, range2)
        return rangeX
    elif QuantileBins:
        rangeX = defineQuantileBinMarker(FCvalue, fixQuantilebinSize)
        return rangeX
    elif EqualBins:
        rangeX = defineEqualBinMarker(FCvalue, fixEqualbinSize)
        if ebinround:
            rangeX = np.round(rangeX, digit)
        return rangeX

def createRawEFtable(rangeX, SGdist, nonSGdist):
    EFtable = {'start': rangeX[:-1], 'end': rangeX[1:]}
    EFtable = pd.DataFrame(EFtable)

    SG_count = [sum((SGdist > EFtable['start'][i]) & (SGdist <= EFtable['end'][i])) for i in range(len(EFtable))]
    nonSG_count = [sum((nonSGdist > EFtable['start'][i]) & (nonSGdist <= EFtable['end'][i])) for i in
                   range(len(EFtable))]
    EFtable["SG_count"] = SG_count
    EFtable["nonSG_count"] = nonSG_count
    return EFtable

def mergeSmallBininEFtable(EFtable, nSG, nnSG, minSGcount=2, lowpctcutoff=0.01, highpctcutoff=0.1):
    ## updates on Sep 30, 2021, add in minimum seed genes count to eliminate over-represented by invidual genes when the seed gene size <100>
    newEFtable = {'start': [], 'end': [], 'SG_count': [], 'nonSG_count': []}
    for j in range(len(EFtable)):
        if j == 0:
            curStart = EFtable['start'][j]
            curEnd = EFtable['end'][j]
            curSG = EFtable['SG_count'][j]
            curnonSG = EFtable['nonSG_count'][j]
        elif ((curnonSG >= nnSG * highpctcutoff) and (curSG >= nSG * lowpctcutoff) and (curSG >= minSGcount)) or ((curSG >= nSG * highpctcutoff) and (curSG >= minSGcount) and (curnonSG >= nnSG * lowpctcutoff)):
            newEFtable['start'].append(curStart)
            newEFtable['end'].append(curEnd)
            newEFtable['SG_count'].append(curSG)
            newEFtable['nonSG_count'].append(curnonSG)

            ## reset tick
            curStart = EFtable['start'][j]
            curEnd = EFtable['end'][j]
            curSG = EFtable['SG_count'][j]
            curnonSG = EFtable['nonSG_count'][j]

        else:
            curSG = curSG + EFtable['SG_count'][j]
            curnonSG = curnonSG + EFtable['nonSG_count'][j]
            curEnd = EFtable['end'][j]

    if ((curnonSG >= nnSG * highpctcutoff) and (curSG >= nSG * lowpctcutoff) and (curSG >= minSGcount)) or (
            (curSG >= nSG * highpctcutoff) and (curnonSG >= nnSG * lowpctcutoff) and (curSG >= minSGcount)) or (
            len(newEFtable['start']) == 0):
        newEFtable['start'].append(curStart)
        newEFtable['end'].append(curEnd)
        newEFtable['SG_count'].append(curSG)
        if curnonSG == 0:
            newEFtable['nonSG_count'].append(1)
        else:
            newEFtable['nonSG_count'].append(curnonSG)
    else:
        newEFtable['end'][-1] = curEnd
        newEFtable['SG_count'][-1] = newEFtable['SG_count'][-1] + curSG
        newEFtable['nonSG_count'][-1] = newEFtable['nonSG_count'][-1] + curnonSG
    newEFtable = pd.DataFrame(newEFtable)
    return newEFtable

def processEFcalculation(EFtable):
    EFtable['rangeX'] = (EFtable['start'] + EFtable['end']) / 2
    EFtable['SG_count'] = EFtable['SG_count'] / sum(EFtable['SG_count'])
    EFtable['nonSG_count'] = EFtable['nonSG_count'] / sum(EFtable['nonSG_count'])
    EFtable['evidenceFactors'] = EFtable['SG_count'] / EFtable['nonSG_count']
    return EFtable

def calculateEF(SGdist, nonSGdist, rangeX, minSGcount=2, lowpctcutoff=0.01, highpctcutoff=0.1):
    nSG = len(SGdist)
    nnSG = len(nonSGdist)
    EFtable = createRawEFtable(rangeX, SGdist, nonSGdist)
    newEFtable = mergeSmallBininEFtable(EFtable, nSG, nnSG, minSGcount, lowpctcutoff, highpctcutoff)
    newEFtable = processEFcalculation(newEFtable)
    newEFtable = newEFtable.dropna()
    return newEFtable

def matchValuetoEF(FCvalue, EFtable):
    FCvalue.columns = ['value']
    FCvalue['evidenceFactors'] = float("nan")
    for i in range(len(EFtable)):
        tmpDict = pd.DataFrame({
            '1': FCvalue['value'] >= EFtable.iloc[i,]['start'],
            '2': FCvalue['value'] < EFtable.iloc[i,]['end']
        })
        FCvalue.loc[tmpDict.all(axis='columns'), 'evidenceFactors'] = EFtable.iloc[i,]['evidenceFactors']
    return FCvalue

def main(dataset, seedGeneList, SUFFIX, outPATH=dirname, writeEFtable=True, EqualBins=True, QuantileBins=True,
         fixEqualbinSize = 100, fixQuantilebinSize = 10, minSGcount=2, lowpctcutoff=0.01, highpctcutoff=0.1):
    sgIndex, nonsgIndex = getIndicesofTwoSet(dataset, seedGeneList)
    datasetEF = pd.DataFrame(dataset['GeneSymbol'].copy())
    featureList = list(dataset.columns)[1:]
    maxEFdict = pd.DataFrame({'maxEF':[0 for i in range(len(featureList))]})
    maxEFdict.index = featureList
    nFeature = 0

    for fc in featureList:
        if nFeature % 1000 == 0:
            current_time = time.strftime("%H:%M:%S", time.localtime())
            print("%s - Processing %s features."%(current_time, nFeature))
        #print("Calculate evidence factors for %s"%(fc))
        if len(dataset.loc[sgIndex, fc].dropna()) == 0 or len(list(set(dataset[fc].dropna()))) == 1:
            print("Skipped %s as no useful information available"%(fc))
            nFeature += 1
            continue

        rangeX = defineBinMarkers(dataset[fc].dropna(), EqualBins=EqualBins, QuantileBins=QuantileBins, fixEqualbinSize = fixEqualbinSize, fixQuantilebinSize = fixQuantilebinSize)
        EFtable = calculateEF(SGdist=dataset.loc[sgIndex, fc].dropna(), nonSGdist=dataset.loc[nonsgIndex, fc].dropna(), rangeX=rangeX, minSGcount=minSGcount, lowpctcutoff=lowpctcutoff, highpctcutoff=highpctcutoff)

        if len(EFtable) == 0:
            print("Skipped %s as insufficient information to calculate EF"%(fc))
            nFeature += 1
            continue
        FCvalue = matchValuetoEF(pd.DataFrame(dataset[fc].copy()), EFtable)
        maxEFdict.loc[fc,"maxEF"] = max(EFtable.loc[:,'evidenceFactors'])
        datasetEF.insert(len(datasetEF.columns), fc, FCvalue.loc[:,'evidenceFactors'], True)
        nFeature += 1

    print("%s - Finish calculate EF for %s features."%(current_time, nFeature))

    if writeEFtable:
        datasetEF.to_csv(outPATH + SUFFIX + "_EFtable.csv", index=False)
    maxEFdict.to_csv(outPATH + SUFFIX + "_maxEF.csv", index=True)

    maxEFtable = pd.DataFrame(maxEFdict)
    return maxEFtable

## test run - this script can replicate calc EF result from original R script
#dataset = pd.read_csv("./DATA/fList_rawValuetable_20210901.csv")
#seedGeneList = pd.read_csv("./DATA/testseedGeneList.txt")['SeedGene'].tolist()
#maxEFtable = main(dataset=dataset, seedGeneList=seedGeneList, SUFFIX='testRun_20210831', outPATH="./OUT/")

