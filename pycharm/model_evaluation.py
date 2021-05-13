# TensorFlow and tf.keras
import tensorflow as tf

# Helper libraries
import random
import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.naive_bayes import GaussianNB
from keras.models import Sequential
from keras.layers import Dense
from sklearn.metrics import confusion_matrix

print(tf.__version__)

def fit_NNmodel(X_train, y_train, n_cols, weight):
    NNmodel = Sequential()
    NNmodel.add(Dense(units=12, activation='relu', input_dim=n_cols))
    NNmodel.add(Dense(units=6, activation='relu'))
    NNmodel.add(Dense(units=1, activation='sigmoid'))
    NNmodel.compile(loss='binary_crossentropy', optimizer='sgd')
    NNmodel.fit(X_train, y_train, validation_split=0.2, class_weight={0: 1, 1: weight}, epochs=100,	batch_size=int(len(y_train) / 5), verbose=0)
    return NNmodel

## Parameter
PATH = "~/OUT/"
INPUT_TABLE = "~/fList_cleanRawValue_withJI_20201229.csv"
dataset = pd.read_csv(INPUT_TABLE)
DATE = "20210512"

labeldata = dataset.copy()
print(labeldata.shape)

X = labeldata.iloc[:, 2:].values
y = labeldata.iloc[:, 1].values

## impute missing value
imputer = SimpleImputer(missing_values=np.nan, strategy="mean")
imputer = imputer.fit(X)
X = imputer.transform(X)

## transfrom label to digit
labelencoder_y = LabelEncoder()
y = labelencoder_y.fit_transform(y)
y_rename = abs(1-y) # label sleep gene as 1 and non sleep gene as 0

# scale features
X = StandardScaler().fit_transform(X)

testRatio = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
nStep = 20

weight = 10000/labeldata[labeldata["label"] == "SRG"].shape[0]
namesCLF = ["LogisticRegression", "NaiveBayes", "DecisionTree", "LinearSVM", "RandomForest", "AdaptiveBoosting"]
classifiers = [
    LogisticRegression(class_weight={0: 1, 1: weight}, max_iter=1000),
    GaussianNB(),
    DecisionTreeClassifier(max_leaf_nodes=12, class_weight={0: 1, 1: weight}),
    SVC(kernel="linear", C=0.025, class_weight={0: 1, 1: weight}),
    RandomForestClassifier(max_leaf_nodes=12, class_weight={0: 1, 1: weight}),
    AdaBoostClassifier(base_estimator=DecisionTreeClassifier(max_leaf_nodes=4, class_weight={0: 1, 1: weight}), algorithm="SAMME")]

totalStep = len(testRatio)*(len(classifiers)+2)*nStep
rawMatrix = [[0 for x in range(6)] for y in range(totalStep)]

all_indices = list(range(X.shape[0]))
seed_value = 42
tf.random.set_seed(seed_value)

predScore = pd.DataFrame(0, index=all_indices, columns=namesCLF + ['NeuralNetwork','ensembleNeuralNetwork'])
predCount = pd.DataFrame(0, index=all_indices, columns=namesCLF + ['NeuralNetwork','ensembleNeuralNetwork'])

k = 0
for i in range(len(testRatio)):
    curRatio = str(round(1-testRatio[i], 1)) + '|' + str(testRatio[i])
    sss = StratifiedShuffleSplit(n_splits=nStep, test_size=testRatio[i], random_state=seed_value)

    j = 0
    for train_index, test_index in sss.split(X, y_rename):
        j += 1
        print("Cycle", j)
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y_rename[train_index], y_rename[test_index]

        # iterate over classifiers
        for name, clf in zip(namesCLF, classifiers):

            try:
                model = clf.fit(X_train, y_train)
                predictions = clf.predict(X_test)
                cm = confusion_matrix(y_test, predictions)
                print(name, curRatio, 'P:', cm[0][1] + cm[1][1], 'FN:', cm[1][0], 'TP:', cm[1][1])
                rawMatrix[k] = [name, curRatio, cm[0][0], cm[0][1], cm[1][0], cm[1][1]]
                k += 1

                if name != 'LinearSVM':
                    predictions = clf.predict_proba(X_test)
                    for zz in range(len(X_test)):
                        predCount.loc[test_index[zz], name] += 1
                        predScore.loc[test_index[zz], name] += predictions[zz][1]

            except ValueError:
                print(name, 'ValueError')
                rawMatrix[k] = [name, curRatio, -1,  -1,  -1,  -1]
                k += 1

        # neural network
        name = 'NeuralNetwork'
        n_cols = X_train.shape[1]
        model = fit_NNmodel(X_train, y_train, n_cols, weight)
        predictions = model.predict_classes(X_test)
        cm = confusion_matrix(y_test, predictions)
        print(name, curRatio, 'P:', cm[0][1] + cm[1][1], 'FN:', cm[1][0], 'TP:', cm[1][1])
        rawMatrix[k] = [name, curRatio, cm[0][0], cm[0][1], cm[1][0], cm[1][1]]
        k += 1

        predictions = model.predict_proba(X_test)
        predictions = predictions.T[0]
        for zz in range(len(X_test)):
            predCount.loc[test_index[zz], name] += 1
            predScore.loc[test_index[zz], name] += predictions[zz]

        ## ensemble neural network
        name = 'ensembleNeuralNetwork'
        n_members = 20
        cutoff = 0.5
        members = [fit_NNmodel(X_train, y_train, n_cols, weight) for _ in range(n_members)]
        yhats = [model.predict_classes(X_test) for model in members]
        yThats = [yh.T[0] for yh in yhats]
        y_pred = [np.zeros(len(yhats[0]))]
        ## combine predictions from all network
        for ii in range(len(yhats[0])):
            for jj in range(n_members):
                y_pred[0][ii] += yThats[jj][ii]
        y_final = np.array([1 if ii > (cutoff * n_members) else 0 for ii in y_pred[0]])
        cm = confusion_matrix(y_test, y_final)
        print(name, curRatio, 'P:', cm[0][1] + cm[1][1], 'FN:', cm[1][0], 'TP:', cm[1][1])
        rawMatrix[k] = [name, curRatio, cm[0][0], cm[0][1], cm[1][0], cm[1][1]]
        k += 1

        predictions = np.array(y_pred[0]/n_members)
        for zz in range(len(X_test)):
            predCount.loc[test_index[zz], name] += 1
            predScore.loc[test_index[zz], name] += predictions[zz]

        tf.keras.backend.clear_session()

for i in range(predScore.shape[1]):
    predScore.iloc[:,i] = predScore.iloc[:,i]/predCount.iloc[:,i]

rawMatrix = pd.DataFrame(rawMatrix)
rawMatrix.to_csv(PATH + "ensembleML_matrixSleepGene_rawValue_sss_" + DATE + ".csv")

predScore.insert(0, "GeneSymbol", list(dataset["GeneSymbol"]), True)
predScore.to_csv(PATH + "predScore_sss_" + DATE + ".csv")

## run algorithm using random label gene
# generate data with random label
testRatio = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
nStep = 10
weight = (100/labeldata[labeldata["label"] == "SRG"].shape[0]) * 100
all_indices = list(range(X.shape[0]))

namesCLF = ["LogisticRegression", "NaiveBayes", "DecisionTree", "LinearSVM",
         "RandomForest", "AdaptiveBoosting"]

classifiers = [
    LogisticRegression(class_weight={0: 1, 1: weight}, max_iter=1000),
    GaussianNB(),
    DecisionTreeClassifier(max_leaf_nodes=12, class_weight={0: 1, 1: weight}),
    SVC(kernel="linear", C=0.025, class_weight={0: 1, 1: weight}),
    RandomForestClassifier(max_leaf_nodes=12, class_weight={0: 1, 1: weight}),
    AdaBoostClassifier(base_estimator=DecisionTreeClassifier(max_depth=1, class_weight={0: 1, 1: weight}), algorithm="SAMME", n_estimators=200)]

nRandomLabelSet = 10
totalStep = len(testRatio)*(len(classifiers)+2)*nStep*nRandomLabelSet
rawMatrix_random = [[0 for x in range(6)] for y in range(totalStep)]
k = 0

nonSRG_indices = np.where(y == 1)[0].tolist()
n_label = len(y[y == 0])
for nn in range(nRandomLabelSet):
    ## Generate random labels
    i_label = sorted(random.sample(nonSRG_indices, n_label))
    y_random = np.asarray([0 for i in range(len(y))])
    for h in i_label:
        y_random[h] = 1
    print(nn, list(dataset.loc[y_random==1,"GeneSymbol"]))

    for i in range(len(testRatio)):
        curRatio = str(round(1 - testRatio[i], 1)) + '|' + str(testRatio[i])
        sss = StratifiedShuffleSplit(n_splits=nStep, test_size=testRatio[i], random_state=seed_value)

        for train_index, test_index in sss.split(X, y_random):
            X_train_random, X_test_random = X[train_index], X[test_index]
            y_train_random, y_test_random = y_random[train_index], y_random[test_index]

            # iterate over classifiers
            for name, clf in zip(namesCLF, classifiers):
                try:
                    model = clf.fit(X_train_random, y_train_random)
                    predictions = clf.predict(X_test_random)
                    cm = confusion_matrix(y_test_random, predictions)
                    print(name, curRatio, 'P:', cm[0][1] + cm[1][1], 'FN:', cm[1][0], 'TP:', cm[1][1])
                    rawMatrix_random[k] = [name, curRatio, cm[0][0], cm[0][1], cm[1][0], cm[1][1]]
                    k += 1
                except ValueError:
                    print(name, 'ValueError')
                    rawMatrix_random[k] = [name, curRatio, -1, -1, -1, -1]
                    k += 1

            name = 'NeuralNetwork'
            n_cols = X_train_random.shape[1]

            model = fit_NNmodel(X_train_random, y_train_random, n_cols, weight)
            predictions = model.predict_classes(X_test_random)
            cm = confusion_matrix(y_test_random, predictions)
            print(name, curRatio, 'P:', cm[0][1] + cm[1][1], 'FN:', cm[1][0], 'TP:', cm[1][1])
            rawMatrix_random[k] = [name, curRatio, cm[0][0], cm[0][1], cm[1][0], cm[1][1]]
            k += 1

            ## ensemble neural network
            name = 'ensembleNeuralNetwork'
            n_members = 20
            cutoff = 0.5
            members = [fit_NNmodel(X_train_random, y_train_random, n_cols, weight) for _ in range(n_members)]
            yhats = [model.predict_classes(X_test_random) for model in members]
            yThats = [yh.T[0] for yh in yhats]
            y_pred = [np.zeros(len(yhats[0]))]
            ## combine predictions from all network
            for ii in range(len(yhats[0])):
                for jj in range(n_members):
                    y_pred[0][ii] += yThats[jj][ii]
            y_final = np.array([1 if ii > (cutoff * n_members) else 0 for ii in y_pred[0]])
            cm = confusion_matrix(y_test_random, y_final)
            print(name, curRatio, 'P:', cm[0][1] + cm[1][1], 'FN:', cm[1][0], 'TP:', cm[1][1])
            rawMatrix_random[k] = [name, curRatio, cm[0][0], cm[0][1], cm[1][0], cm[1][1]]
            k += 1

            tf.keras.backend.clear_session()

rawMatrix_random = pd.DataFrame(rawMatrix_random)
rawMatrix_random.to_csv(PATH + "ensembleML_matrixRandomLabel_rawValue_sss_" + DATE + ".csv")
