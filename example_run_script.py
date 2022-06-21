## Updated protocol to run the model - Sept 23, 2021
## Major modifications - calculate JI first and select only JI features with positive evidence to train the prediction model.

import os
import yaml
import pandas as pd
import numpy as np
import calc_JaccardIndex as calcJI
import calc_EF as calcEF
import function as SGfn
import ML_model as mlm

## parameters - file names
yamlFN = "./DATA/parameter.yaml"
outPATH = "./OUT/"               ## string: directory to store output file
SUFFIX = 'exampleRun'                          ## string: suffix added to out file

seedGeneFN = "./DATA/seedGene.txt"          ## string: filename with seed gene list, one gene(human gene name) per line, line started with # is skipped.
inputTable = "./DATA/testInput.csv"         ## string: filename with genes as row and features as column

model_eval = False           ## boolean:True/False, whether to run model evaluation using seed genes
random_model_eval = False    ## boolean:True/False, whether to run model evaluation using random labels
RF_prediction = True        ## boolean:True/False, whether to run random forest prediction

with open(yamlFN) as file:
    # The FullLoader parameter handles the conversion from YAML
    # scalar values to Python the dictionary format
    parameter = yaml.load(file, Loader=yaml.FullLoader)

    print(parameter)

##### starting of the script ######
## read features and seed genes
dataset = pd.read_csv(inputTable)
seedGeneList = SGfn.readSeedGeneList(seedGeneFN)

## calculate evidence factors
if os.path.exists(outPATH + SUFFIX + "_maxEF.csv"):
    maxEFtable = pd.read_csv(outPATH + SUFFIX + "_maxEF.csv", index_col=0)
else:
    maxEFtable = calcEF.main(dataset=dataset, seedGeneList=seedGeneList, SUFFIX=SUFFIX, outPATH=outPATH,
                             writeEFtable=parameter['calcEF_writeEFtable'], EqualBins=parameter['calcEF_EqualBins'],
                             QuantileBins=parameter['calcEF_QuantileBins'],
                             fixEqualbinSize=parameter['calcEF_fixEqualbinSize'],
                             fixQuantilebinSize=parameter['calcEF_fixQuantilebinSize'],
                             minSGcount=parameter['calcEF_minSGcount'], lowpctcutoff=parameter['calcEF_lowpctcutoff'],
                             highpctcutoff=parameter['calcEF_highpctcutoff'])
selectedFC = list(maxEFtable[maxEFtable["maxEF"]>=parameter['maxEFcutoff']].index)

if os.path.exists(outPATH + SUFFIX + "_cleanDataSet.csv"):
    cleanDataset = pd.read_csv(outPATH + SUFFIX + "_cleanDataSet.csv")
else:
    ## select features and genes for prediction, seed genes might be filtered out in this step if proportion of column with missing value > rowCutOff.
    ## To keep certain genes in the table regardless the cutoff value, list the HGNC symbol in 'exemptGeneList'
    ## e.g. exemptGeneList = np.array(['NPSR1'])
    exemptGeneList = np.array(['NPSR1'])
    cleanDataset = SGfn.filterTable(dataset, selectedFC, rowCutOff=parameter['propCutOffForGeneswithNA'], exemptGeneList=exemptGeneList)

    ## run Pearson correlation, if two features found to have correlation >corCutoff (default=0.8), remove the feature with lower maxEF
    cleanDataset = SGfn.filterbyCor(cleanDataset, maxEFtable, corCutoff=parameter['corCutOffForFeatures'])

    ## add jaccard index feature based on seed gene list
    cleanDataset = calcJI.main(dataset=cleanDataset, geneSetDataList="./DATA/gmt_file/geneset.list", seedGeneList=seedGeneList, SUFFIX=SUFFIX, outPATH=outPATH, inserttodataframe=True, averagesimilarterm=True, simCutoff=0.5)
    cleanDataset.to_csv(outPATH + SUFFIX + "_cleanDataSet.csv", index=False)

## prepare X-features and y-label for model evaluation
X, y = SGfn.extract_info(cleanDataset, seedGeneList)
X = SGfn.preprocess_feature(X, imputeData=True, scaleFeature=True)

if model_eval:
    ## model evaluation
    seedGeneWeight = (100/len(y[y==1])) * 100
    modelEval_rawMatrix, modelEval_valMatrix, modelEval_predScore = mlm.model_eval(X, y, testRatioStep=parameter['RFpred_predRatioStep'], iter=5, seedGeneWeight=seedGeneWeight, seed_value=parameter['RFpred_seed_value'],
                   LR_max_iter=1000, DT_max_leaf_nodes=12, SVC_C=0.025, RF_max_leaf_nodes=12, ADB_max_leaf_nodes=4, eNN_iter=20, eNN_cutoff=0.5)
    modelEval_predScore.insert(0, "GeneSymbol", list(cleanDataset["GeneSymbol"]), True)
    modelEval_rawMatrix.to_csv(outPATH + SUFFIX + "_ensembleML.csv")
    modelEval_valMatrix.to_csv(outPATH + SUFFIX + "_val.ensembleML.csv")
    modelEval_predScore.to_csv(outPATH + SUFFIX + "_predScore.csv")

if random_model_eval:
    # evaluate model with random labels
    seedGeneWeight = (100/len(y[y==1])) * 100
    modelEval_randomLabelMatrix, modelEval_randomLabelvalMatrix = mlm.random_label_eval(X, y, testRatioStep=parameter['RFpred_predRatioStep'], nRandomLabelSet=10, iter=10, seedGeneWeight=seedGeneWeight, seed_value=parameter['RFpred_seed_value'],
                   LR_max_iter=1000, DT_max_leaf_nodes=12, SVC_C=0.025, RF_max_leaf_nodes=12, ADB_max_leaf_nodes=4, eNN_iter=20, eNN_cutoff=0.5)
    modelEval_randomLabelMatrix.to_csv(outPATH + SUFFIX + "_ensembleML_RandomLabel.csv")
    modelEval_randomLabelvalMatrix.to_csv(outPATH + SUFFIX + "_val.ensembleML_RandomLabel.csv")

if RF_prediction:
    ## predict with RF model
    seedGeneWeight = (100/len(y[y==1])) * 100   # define weight of the seed genes
    RFPredTable = mlm.RF_pred(cleanDataset, seedGeneList, predRatioStep=parameter['RFpred_predRatioStep'], iter_n=parameter['RFpred_iter_n'], seedGeneWeight=seedGeneWeight, seed_value=parameter['RFpred_seed_value'], dynamicRFmaxLeafNodes=parameter['RFpred_dynamicRFmaxLeafNodes'], RFmaxLeaftoLabelsRatio=parameter['RFpred_RFmaxLeaftoLabelsRatio'], SUFFIX=SUFFIX, outPATH=outPATH)
    process_score = SGfn.calcFinalPredScore(RFPredTable, seedGeneList, bgnoisethreshold=parameter['RFpred_bgnoisethreshold'], highCthreshold=parameter['RFpred_highCthreshold'], predRatioStep=parameter['RFpred_predRatioStep'])
    process_score.to_csv(outPATH + SUFFIX + "_finalScore.csv")
