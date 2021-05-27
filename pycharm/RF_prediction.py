# Helper libraries
import numpy as np
import pandas as pd
import random
import calc_JaccardIndex as calcJI
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.ensemble import RandomForestClassifier

PATH = "../OUT/"
DATE = "20210503"
INPUT_TABLE = "../OUT/fList_cleanRawValue_" + DATE + ".csv"
runJI = True    ## Do you want to add in features from jaccard index in this script?

if runJI == True:
    inData = pd.read_csv(INPUT_TABLE)
    dataset = calcJI.runJI(dataset = inData,
                          curwkPATH = PATH,
                          geneSetList = "../DATA/gmt_file/geneset.list",
                          SUFFIX = DATE,
                          DATAPATH="../DATA/")
else:
    dataset = pd.read_csv(INPUT_TABLE)

labeldata = dataset.copy()
print("Input dimensions: (rows, columns)", labeldata.shape)

X = labeldata.iloc[:, 2:].values
y = labeldata.iloc[:, 1].values

## impute missing value
imputer = SimpleImputer(missing_values=np.nan, strategy="mean")
imputer = imputer.fit(X)
X = imputer.transform(X)

## transfrom label to digit
labelencoder_y = LabelEncoder()
y = labelencoder_y.fit_transform(y)
y_rename = abs(1-y)  # label sleep gene as 1 and non sleep gene as 0

## Create correlation table between sleep genes
testRatio = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
sleepGeneList = list(dataset["GeneSymbol"][y==0])
weight = (10000/labeldata[labeldata["label"] == "SRG"].shape[0])
all_indices = list(range(X.shape[0]))
feature_list = dataset.columns[2:]
seed_value = 18
nStep = 100
SUFFIX = "n" + str(nStep) + "_" + DATE + "_"

featureDict = {}
for f in dataset.columns.values[2:]:
    featureDict[f] = {}

for i in range(len(testRatio)):
    curRatio = str(round(1-testRatio[i], 1)) + '|' + str(testRatio[i])
    sss = StratifiedShuffleSplit(n_splits=nStep, test_size=testRatio[i], random_state=seed_value)
    print(curRatio)

    new_weight = weight * (1.2-testRatio[i])
    genepredictTableRFpos = pd.DataFrame(0, index=sleepGeneList, columns=sleepGeneList)
    genepredictTableRFneg = pd.DataFrame(0, index=sleepGeneList, columns=sleepGeneList)

    traintesttable = dataset[['GeneSymbol', 'label']].copy()
    RFtable = dataset[['GeneSymbol', 'label']].copy()

    j = 0
    for train_index, test_index in sss.split(X, y_rename):
        j += 1
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y_rename[train_index], y_rename[test_index]

        testcase = np.array([1 if i in test_index else 0 for i in range(len(y))])
        traintesttable.insert(2, "testcase_"+str(j), testcase, True)

        trainSleepGene = list(dataset.loc[train_index, "GeneSymbol"][y_train == 1])
        testSleepGene = list(dataset.loc[test_index, "GeneSymbol"][y_test == 1])

        name = "RandomForest"
        rf = RandomForestClassifier(max_leaf_nodes=12, class_weight={1: weight, 0: 1})
        rf.fit(X_train, y_train)
        predictions = rf.predict_proba(X_test)
        y_RF = rf.predict_proba(X)
        y_predictSRG = [round(i[1], 4) for i in y_RF]
        RFtable.insert(2, "pred__"+str(j), y_predictSRG, True)

        # feature importances
        importances = list(rf.feature_importances_)
        # List of tuples with variable and importance
        feature_importances = [(feature, round(importance, 4)) for feature, importance in zip(feature_list, importances)]
        # Sort the feature importances by most important first
        feature_importances = sorted(feature_importances, key=lambda x: x[1], reverse=True)

        # Record features have effect >0
        for k in range(len(feature_importances)):
            if feature_importances[k][1] > 0:
                try:
                    featureDict[feature_importances[k][0]][testRatio[i]].append(feature_importances[k][1])
                except KeyError:
                    featureDict[feature_importances[k][0]][testRatio[i]] = [feature_importances[k][1]]

    traintesttable.to_csv(PATH + "traintesttable_" + SUFFIX + str(testRatio[i])+".csv")
    RFtable.to_csv(PATH + "RFtable_" + SUFFIX + str(testRatio[i])+".csv")

wf = open(PATH + "feature_importance_" + DATE + ".txt", "w")
for f in featureDict.keys():
    #wf.write('%s\n' % f)
    for j in featureDict[f].keys():
        wf.write('%s\t%s\t'%(f,j))
        for k in featureDict[f][j]:
            wf.write(',%s'%k)
        wf.write('\n')
wf.close()