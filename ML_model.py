# TensorFlow and tf.keras
import tensorflow as tf

# Helper libraries
import os
import random
import math
import time
import numpy as np
import pandas as pd
import function as SGfn
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.naive_bayes import GaussianNB
from keras.models import Sequential
from keras.layers import Dense
from sklearn.metrics import confusion_matrix
from sklearn.decomposition import PCA

print('Tensorflow version:', tf.__version__)
dirname = os.path.dirname(__file__)

def fit_NNmodel(X_train, y_train, n_cols, weight):
    NNmodel = Sequential()
    NNmodel.add(Dense(units=12, activation='relu', input_dim=n_cols))
    NNmodel.add(Dense(units=6, activation='relu'))
    NNmodel.add(Dense(units=1, activation='sigmoid'))
    NNmodel.compile(loss='binary_crossentropy', optimizer='sgd')
    NNmodel.fit(X_train, y_train, validation_split=0.2, class_weight={0: 1, 1: weight}, epochs=100,	batch_size=int(len(y_train) / 5), verbose=0)
    return NNmodel

def model_eval(X, y, testRatioStep=[0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8], iter=20, seedGeneWeight=100, seed_value=42, LR_max_iter=1000,
               DT_max_leaf_nodes=12, SVC_C=0.025, RF_max_leaf_nodes=12, runRFoob=True, ADB_max_leaf_nodes=4, cutOffForPosLabel=0.5, eNN_iter=20, eNN_cutoff=0.5):

    ## principal component analysis - for Naive Bayes input
    pca = PCA()
    X_PC = pca.fit_transform(X)

    ## run model evaluation using 6 machine learning algorithm
    namesCLF = ["LogisticRegression", "NaiveBayes", "DecisionTree", "LinearSVM", "RandomForest", "AdaptiveBoosting"]
    classifiers = [
        LogisticRegression(class_weight={0: 1, 1: seedGeneWeight}, max_iter=LR_max_iter, random_state=seed_value),
        GaussianNB(),
        DecisionTreeClassifier(max_leaf_nodes=DT_max_leaf_nodes, class_weight={0: 1, 1: seedGeneWeight}, random_state=seed_value),
        SVC(kernel="linear", C=SVC_C, class_weight={0: 1, 1: seedGeneWeight}, random_state=seed_value),
        RandomForestClassifier(max_leaf_nodes=RF_max_leaf_nodes, oob_score=runRFoob, class_weight={0: 1, 1: seedGeneWeight}, random_state=seed_value),
        AdaBoostClassifier(base_estimator=DecisionTreeClassifier(max_leaf_nodes=ADB_max_leaf_nodes, class_weight={0: 1, 1: seedGeneWeight}, random_state=seed_value), algorithm="SAMME", random_state=seed_value)]

    tf.random.set_seed(seed_value)

    totalStep = len(testRatioStep)*(len(classifiers)+2)*iter
    rawMatrix = [[0 for x in range(6)] for y in range(totalStep)]
    valMatrix = [[0 for x in range(6)] for y in range(totalStep)]

    all_indices = list(range(X.shape[0]))

    predScore = pd.DataFrame(0, index=all_indices, columns=namesCLF + ['NeuralNetwork','ensembleNeuralNetwork'])
    predCount = pd.DataFrame(0, index=all_indices, columns=namesCLF + ['NeuralNetwork','ensembleNeuralNetwork'])

    k = 0
    for i in range(len(testRatioStep)):
        curRatio = str(round(1-testRatioStep[i], 1)) + '|' + str(testRatioStep[i])
        sss = StratifiedShuffleSplit(n_splits=iter, test_size=testRatioStep[i], random_state=seed_value)

        j = 0
        for train_index, test_index in sss.split(X, y):
            j += 1
            print("Cycle", j)
            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = y[train_index], y[test_index]

            # iterate over classifiers
            for name, clf in zip(namesCLF, classifiers):
                try:
                    if name == "NaiveBayes":
                        ## for naive bayes, used principal componenet transformed input to train the model
                        X_PC_train, X_PC_test = X_PC[train_index], X_PC[test_index]

                        model = clf.fit(X_PC_train, y_train)
                        predictions = clf.predict(X_PC_test)
                        cm = confusion_matrix(y_test, predictions)
                        print(name, curRatio, 'P:', cm[0][1] + cm[1][1], 'FN:', cm[1][0], 'TP:', cm[1][1])
                        rawMatrix[k] = [name, curRatio, cm[0][0], cm[0][1], cm[1][0], cm[1][1]]

                        predictions = clf.predict_proba(X_PC_test)
                        for zz in range(len(X_PC_test)):
                            predCount.loc[test_index[zz], name] += 1
                            predScore.loc[test_index[zz], name] += predictions[zz][1]

                        predictions = clf.predict(X_PC_train)
                        cm = confusion_matrix(y_train, predictions)
                        print(name, curRatio, 'P:', cm[0][1] + cm[1][1], 'FN:', cm[1][0], 'TP:', cm[1][1])
                        valMatrix[k] = [name, curRatio, cm[0][0], cm[0][1], cm[1][0], cm[1][1]]
                        k += 1

                    else:
                        model = clf.fit(X_train, y_train)
                        predictions = clf.predict(X_test)
                        cm = confusion_matrix(y_test, predictions)
                        print(name, curRatio, 'P:', cm[0][1] + cm[1][1], 'FN:', cm[1][0], 'TP:', cm[1][1])
                        rawMatrix[k] = [name, curRatio, cm[0][0], cm[0][1], cm[1][0], cm[1][1]]

                        if name != 'LinearSVM':
                            predictions = clf.predict_proba(X_test)
                            for zz in range(len(X_test)):
                                predCount.loc[test_index[zz], name] += 1
                                predScore.loc[test_index[zz], name] += predictions[zz][1]

                        predictions = clf.predict(X_train)
                        cm = confusion_matrix(y_train, predictions)
                        valMatrix[k] = [name, curRatio, cm[0][0], cm[0][1], cm[1][0], cm[1][1]]
                        k += 1

                except ValueError:
                    print(name, 'ValueError')
                    rawMatrix[k] = [name, curRatio, -1,  -1,  -1,  -1]
                    valMatrix[k] = [name, curRatio, -1,  -1,  -1,  -1]
                    k += 1

            # neural network
            name = 'NeuralNetwork'
            n_cols = X_train.shape[1]
            model = fit_NNmodel(X_train, y_train, n_cols, seedGeneWeight)
            predict_x = model.predict(X_test)
            predictions = np.array([1 if px>=cutOffForPosLabel else 0 for px in predict_x])

            cm = confusion_matrix(y_test, predictions)
            print(name, curRatio, 'P:', cm[0][1] + cm[1][1], 'FN:', cm[1][0], 'TP:', cm[1][1])
            rawMatrix[k] = [name, curRatio, cm[0][0], cm[0][1], cm[1][0], cm[1][1]]

            predictions = predict_x.T[0]
            for zz in range(len(X_test)):
                predCount.loc[test_index[zz], name] += 1
                predScore.loc[test_index[zz], name] += predictions[zz]

            predict_x = model.predict(X_train)
            predictions = np.array([1 if px>=cutOffForPosLabel else 0 for px in predict_x])
            cm = confusion_matrix(y_train, predictions)
            valMatrix[k] = [name, curRatio, cm[0][0], cm[0][1], cm[1][0], cm[1][1]]

            k += 1

            ## ensemble neural network
            name = 'ensembleNeuralNetwork'
            members = [fit_NNmodel(X_train, y_train, n_cols, seedGeneWeight) for _ in range(eNN_iter)]
            yhats = [model.predict(X_test) for model in members]
            yThats = [np.where(yh>=cutOffForPosLabel,1,0) for yh in yhats]     ## predict class by larger or less than 0.5 probability
            yThats = [yh.T[0] for yh in yThats]         ## transform the array
            y_pred = np.sum(yThats, axis=0)
            y_final = np.array([1 if ii > (eNN_cutoff * eNN_iter) else 0 for ii in y_pred])

            cm = confusion_matrix(y_test, y_final)
            print(name, curRatio, 'P:', cm[0][1] + cm[1][1], 'FN:', cm[1][0], 'TP:', cm[1][1])
            rawMatrix[k] = [name, curRatio, cm[0][0], cm[0][1], cm[1][0], cm[1][1]]

            predictions = np.array(y_pred/eNN_iter)
            for zz in range(len(X_test)):
                predCount.loc[test_index[zz], name] += 1
                predScore.loc[test_index[zz], name] += predictions[zz]

            yhats = [model.predict(X_train) for model in members]
            yThats = [np.where(yh>=cutOffForPosLabel,1,0) for yh in yhats]     ## predict class by larger or less than 0.5 probability
            yThats = [yh.T[0] for yh in yThats]         ## transform the array
            y_pred = np.sum(yThats, axis=0)
            y_final = np.array([1 if ii > (eNN_cutoff * eNN_iter) else 0 for ii in y_pred])

            cm = confusion_matrix(y_train, y_final)
            valMatrix[k] = [name, curRatio, cm[0][0], cm[0][1], cm[1][0], cm[1][1]]

            k += 1

            tf.keras.backend.clear_session()

    for i in range(predScore.shape[1]):
        predScore.iloc[:,i] = predScore.iloc[:,i]/predCount.iloc[:,i]

    rawMatrix = pd.DataFrame(rawMatrix)
    valMatrix = pd.DataFrame(valMatrix)
    predScore = pd.DataFrame(predScore)

    return rawMatrix, valMatrix, predScore


def random_label_eval(X, y, testRatioStep=[0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8], nRandomLabelSet=10, iter=10, seedGeneWeight=100, seed_value=42, LR_max_iter=1000, DT_max_leaf_nodes=12, SVC_C=0.025, RF_max_leaf_nodes=12, runRFoob=True, ADB_max_leaf_nodes=4, cutOffForPosLabel=0.5, eNN_iter=20, eNN_cutoff=0.5):
    ## principal component analysis - for Naive Bayes input
    pca = PCA()
    X_PC = pca.fit_transform(X)

    ## run algorithm using random label gene
    namesCLF = ["LogisticRegression", "NaiveBayes", "DecisionTree", "LinearSVM", "RandomForest", "AdaptiveBoosting"]
    classifiers = [
        LogisticRegression(class_weight={0: 1, 1: seedGeneWeight}, max_iter=LR_max_iter, random_state=seed_value),
        GaussianNB(),
        DecisionTreeClassifier(max_leaf_nodes=DT_max_leaf_nodes, class_weight={0: 1, 1: seedGeneWeight}, random_state=seed_value),
        SVC(kernel="linear", C=SVC_C, class_weight={0: 1, 1: seedGeneWeight}, random_state=seed_value),
        RandomForestClassifier(max_leaf_nodes=RF_max_leaf_nodes, oob_score=runRFoob, class_weight={0: 1, 1: seedGeneWeight}, random_state=seed_value),
        AdaBoostClassifier(base_estimator=DecisionTreeClassifier(max_leaf_nodes=ADB_max_leaf_nodes, class_weight={0: 1, 1: seedGeneWeight}, random_state=seed_value), algorithm="SAMME", random_state=seed_value)]

    tf.random.set_seed(seed_value)

    # generate data with random label
    totalStep = len(testRatioStep)*(len(classifiers)+2)*iter*nRandomLabelSet
    rawMatrix_random = [[0 for x in range(6)] for y in range(totalStep)]
    valMatrix_random = [[0 for x in range(6)] for y in range(totalStep)]

    k = 0

    nonSRG_indices = np.where(y == 0)[0].tolist()
    n_label = len(y[y == 1]) # number of seed genes
    for nn in range(nRandomLabelSet):
        ## Generate random labels
        current_time = time.strftime("%H:%M:%S", time.localtime())
        print("%s - Running set %s of random labels"%(current_time, nn))

        i_label = sorted(random.sample(nonSRG_indices, n_label))
        y_random = np.asarray([0 for i in range(len(y))])
        for h in i_label:
            y_random[h] = 1

        for i in range(len(testRatioStep)):
            curRatio = str(round(1 - testRatioStep[i], 1)) + '|' + str(testRatioStep[i])
            sss = StratifiedShuffleSplit(n_splits=iter, test_size=testRatioStep[i], random_state=seed_value)
            j = 0

            for train_index, test_index in sss.split(X, y_random):
                j += 1
                print("Random sets:%s, test-ratio:%s, Cycle:%s"%(nn,testRatioStep[i], j))

                X_train_random, X_test_random = X[train_index], X[test_index]
                y_train_random, y_test_random = y_random[train_index], y_random[test_index]

                # iterate over classifiers
                for name, clf in zip(namesCLF, classifiers):
                    try:
                        if name == "NaiveBayes":
                            ## for naive bayes, used principal componenet transformed input to train the model
                            X_PC_train, X_PC_test = X_PC[train_index], X_PC[test_index]

                            model = clf.fit(X_PC_train, y_train_random)
                            predictions = clf.predict(X_PC_test)
                            cm = confusion_matrix(y_test_random, predictions)
                            print(name, curRatio, 'P:', cm[0][1] + cm[1][1], 'FN:', cm[1][0], 'TP:', cm[1][1])
                            rawMatrix_random[k] = [name, curRatio, cm[0][0], cm[0][1], cm[1][0], cm[1][1]]

                            predictions = clf.predict(X_train_random)
                            cm = confusion_matrix(y_train_random, predictions)
                            valMatrix_random[k] = [name, curRatio, cm[0][0], cm[0][1], cm[1][0], cm[1][1]]

                            k += 1

                        else:
                            model = clf.fit(X_train_random, y_train_random)
                            predictions = clf.predict(X_test_random)
                            cm = confusion_matrix(y_test_random, predictions)
                            print(name, curRatio, 'P:', cm[0][1] + cm[1][1], 'FN:', cm[1][0], 'TP:', cm[1][1])
                            rawMatrix_random[k] = [name, curRatio, cm[0][0], cm[0][1], cm[1][0], cm[1][1]]

                            predictions = clf.predict(X_train_random)
                            cm = confusion_matrix(y_train_random, predictions)
                            valMatrix_random[k] = [name, curRatio, cm[0][0], cm[0][1], cm[1][0], cm[1][1]]

                            k += 1

                    except ValueError:
                        print(name, 'ValueError')
                        rawMatrix_random[k] = [name, curRatio, -1, -1, -1, -1]
                        valMatrix_random[k] = [name, curRatio, -1, -1, -1, -1]
                        k += 1

                name = 'NeuralNetwork'
                n_cols = X_train_random.shape[1]
                model = fit_NNmodel(X_train_random, y_train_random, n_cols, seedGeneWeight)
                predict_x = model.predict(X_test_random)
                predictions = np.array([1 if px >= cutOffForPosLabel else 0 for px in predict_x])

                cm = confusion_matrix(y_test_random, predictions)
                print(name, curRatio, 'P:', cm[0][1] + cm[1][1], 'FN:', cm[1][0], 'TP:', cm[1][1])
                rawMatrix_random[k] = [name, curRatio, cm[0][0], cm[0][1], cm[1][0], cm[1][1]]

                predict_x = model.predict(X_train_random)
                predictions = np.array([1 if px >= cutOffForPosLabel else 0 for px in predict_x])
                cm = confusion_matrix(y_train_random, predictions)
                valMatrix_random[k] = [name, curRatio, cm[0][0], cm[0][1], cm[1][0], cm[1][1]]

                k += 1

                name = 'ensembleNeuralNetwork'
                members = [fit_NNmodel(X_train_random, y_train_random, n_cols, seedGeneWeight) for _ in range(eNN_iter)]
                yhats = [model.predict(X_test_random) for model in members]
                yThats = [np.where(yh >= cutOffForPosLabel, 1, 0) for yh in yhats]  ## predict class by larger or less than 0.5 probability
                yThats = [yh.T[0] for yh in yThats]  ## transform the array
                y_pred = np.sum(yThats, axis=0)
                y_final = np.array([1 if ii > (eNN_cutoff * eNN_iter) else 0 for ii in y_pred])

                cm = confusion_matrix(y_test_random, y_final)
                print(name, curRatio, 'P:', cm[0][1] + cm[1][1], 'FN:', cm[1][0], 'TP:', cm[1][1])
                rawMatrix_random[k] = [name, curRatio, cm[0][0], cm[0][1], cm[1][0], cm[1][1]]

                yhats = [model.predict(X_train_random) for model in members]
                yThats = [np.where(yh >= cutOffForPosLabel, 1, 0) for yh in yhats]  ## predict class by larger or less than 0.5 probability
                yThats = [yh.T[0] for yh in yThats]  ## transform the array
                y_pred = np.sum(yThats, axis=0)
                y_final = np.array([1 if ii > (eNN_cutoff * eNN_iter) else 0 for ii in y_pred])

                cm = confusion_matrix(y_train_random, y_final)
                valMatrix_random[k] = [name, curRatio, cm[0][0], cm[0][1], cm[1][0], cm[1][1]]
                k += 1

                tf.keras.backend.clear_session()

    rawMatrix_random = pd.DataFrame(rawMatrix_random)
    valMatrix_random = pd.DataFrame(valMatrix_random)

    return rawMatrix_random, valMatrix_random

def writeFeatureImportance(featureDict, filename="feature_importance.csv"):
    wf = open(filename, "w")
    for f in featureDict.keys():
        for j in featureDict[f].keys():
            wf.write('%s\t%s\t' % (f, j))
            for k in featureDict[f][j]:
                wf.write(',%s' % k)
            wf.write('\n')
    wf.close()


def model_eval_ttv(X, y, testRatioStep=[0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8], iter=20, seedGeneWeight=100, seed_value=42, LR_max_iter=1000,
               DT_max_leaf_nodes=12, SVC_C=0.025, RF_max_leaf_nodes=12, runRFoob=True, ADB_max_leaf_nodes=4, cutOffForPosLabel=0.5, eNN_iter=20, eNN_cutoff=0.5):

    ## principal component analysis - for Naive Bayes input
    pca = PCA()
    X_PC = pca.fit_transform(X)

    ## run model evaluation using 6 machine learning algorithm
    namesCLF = ["LogisticRegression", "NaiveBayes", "DecisionTree", "LinearSVM", "RandomForest", "AdaptiveBoosting"]
    classifiers = [
        LogisticRegression(class_weight={0: 1, 1: seedGeneWeight}, max_iter=LR_max_iter, random_state=seed_value),
        GaussianNB(),
        DecisionTreeClassifier(max_leaf_nodes=DT_max_leaf_nodes, class_weight={0: 1, 1: seedGeneWeight}, random_state=seed_value),
        SVC(kernel="linear", C=SVC_C, class_weight={0: 1, 1: seedGeneWeight}, random_state=seed_value),
        RandomForestClassifier(max_leaf_nodes=RF_max_leaf_nodes, oob_score=runRFoob, class_weight={0: 1, 1: seedGeneWeight}, random_state=seed_value),
        AdaBoostClassifier(base_estimator=DecisionTreeClassifier(max_leaf_nodes=ADB_max_leaf_nodes, class_weight={0: 1, 1: seedGeneWeight}, random_state=seed_value), algorithm="SAMME", random_state=seed_value)]

    tf.random.set_seed(seed_value)

    totalStep = len(testRatioStep)*(len(classifiers)+2)*iter
    rawMatrix = [[0 for x in range(6)] for y in range(totalStep)]
    valMatrix = [[0 for x in range(6)] for y in range(totalStep)]

    all_indices = list(range(X.shape[0]))

    predScore = pd.DataFrame(0, index=all_indices, columns=namesCLF + ['NeuralNetwork','ensembleNeuralNetwork'])
    predCount = pd.DataFrame(0, index=all_indices, columns=namesCLF + ['NeuralNetwork','ensembleNeuralNetwork'])

    k = 0
    for aa in range(1):
        ssss = StratifiedShuffleSplit(n_splits=1, test_size=0.2, random_state=seed_value)
        for sel_index, val_index in ssss.split(X, y):
            X_selected, X_val = X[sel_index], X[val_index]
            y_selected, y_val = y[sel_index], y[val_index]

        for i in range(len(testRatioStep)):
            curRatio = str(round(1-testRatioStep[i], 1)) + '|' + str(testRatioStep[i])
            sss = StratifiedShuffleSplit(n_splits=iter, test_size=testRatioStep[i], random_state=seed_value)

            j = 0
            for train_index, test_index in sss.split(X_selected, y_selected):
                j += 1
                print("Cycle", j)
                X_train, X_test = X_selected[train_index], X_selected[test_index]
                y_train, y_test = y_selected[train_index], y_selected[test_index]

                # iterate over classifiers
                for name, clf in zip(namesCLF, classifiers):
                    try:
                        if name == "NaiveBayes":
                            ## for naive bayes, used principal componenet transformed input to train the model
                            X_PC_train, X_PC_test = X_PC[sel_index][train_index], X_PC[sel_index][test_index]

                            model = clf.fit(X_PC_train, y_train)
                            predictions = clf.predict(X_PC_test)
                            cm = confusion_matrix(y_test, predictions)
                            print(name, curRatio, 'P:', cm[0][1] + cm[1][1], 'FN:', cm[1][0], 'TP:', cm[1][1])
                            rawMatrix[k] = [name, curRatio, cm[0][0], cm[0][1], cm[1][0], cm[1][1]]

                            predictions = clf.predict_proba(X_PC_test)
                            for zz in range(len(X_PC_test)):
                                predCount.loc[test_index[zz], name] += 1
                                predScore.loc[test_index[zz], name] += predictions[zz][1]

                            predictions = clf.predict(X_PC_train)
                            cm = confusion_matrix(y_train, predictions)
                            print(name, curRatio, 'P:', cm[0][1] + cm[1][1], 'FN:', cm[1][0], 'TP:', cm[1][1])
                            valMatrix[k] = [name, curRatio, cm[0][0], cm[0][1], cm[1][0], cm[1][1]]
                            k += 1

                        else:
                            model = clf.fit(X_train, y_train)
                            predictions = clf.predict(X_test)
                            cm = confusion_matrix(y_test, predictions)
                            print(name, curRatio, 'P:', cm[0][1] + cm[1][1], 'FN:', cm[1][0], 'TP:', cm[1][1])
                            rawMatrix[k] = [name, curRatio, cm[0][0], cm[0][1], cm[1][0], cm[1][1]]

                            if name != 'LinearSVM':
                                predictions = clf.predict_proba(X_test)
                                for zz in range(len(X_test)):
                                    predCount.loc[test_index[zz], name] += 1
                                    predScore.loc[test_index[zz], name] += predictions[zz][1]

                            predictions = clf.predict(X_val)
                            cm = confusion_matrix(y_val, predictions)
                            valMatrix[k] = [name, curRatio, cm[0][0], cm[0][1], cm[1][0], cm[1][1]]
                            k += 1

                    except ValueError:
                        print(name, 'ValueError')
                        rawMatrix[k] = [name, curRatio, -1,  -1,  -1,  -1]
                        valMatrix[k] = [name, curRatio, -1,  -1,  -1,  -1]
                        k += 1

                # neural network
                name = 'NeuralNetwork'
                n_cols = X_train.shape[1]
                model = fit_NNmodel(X_train, y_train, n_cols, seedGeneWeight)
                predict_x = model.predict(X_test)
                predictions = np.array([1 if px>=cutOffForPosLabel else 0 for px in predict_x])

                cm = confusion_matrix(y_test, predictions)
                print(name, curRatio, 'P:', cm[0][1] + cm[1][1], 'FN:', cm[1][0], 'TP:', cm[1][1])
                rawMatrix[k] = [name, curRatio, cm[0][0], cm[0][1], cm[1][0], cm[1][1]]

                predictions = predict_x.T[0]
                for zz in range(len(X_test)):
                    predCount.loc[test_index[zz], name] += 1
                    predScore.loc[test_index[zz], name] += predictions[zz]

                predict_x = model.predict(X_val)
                predictions = np.array([1 if px>=cutOffForPosLabel else 0 for px in predict_x])
                cm = confusion_matrix(y_val, predictions)
                valMatrix[k] = [name, curRatio, cm[0][0], cm[0][1], cm[1][0], cm[1][1]]

                k += 1

                ## ensemble neural network
                name = 'ensembleNeuralNetwork'
                members = [fit_NNmodel(X_train, y_train, n_cols, seedGeneWeight) for _ in range(eNN_iter)]
                yhats = [model.predict(X_test) for model in members]
                yThats = [np.where(yh>=cutOffForPosLabel,1,0) for yh in yhats]     ## predict class by larger or less than 0.5 probability
                yThats = [yh.T[0] for yh in yThats]         ## transform the array
                y_pred = np.sum(yThats, axis=0)
                y_final = np.array([1 if ii > (eNN_cutoff * eNN_iter) else 0 for ii in y_pred])

                cm = confusion_matrix(y_test, y_final)
                print(name, curRatio, 'P:', cm[0][1] + cm[1][1], 'FN:', cm[1][0], 'TP:', cm[1][1])
                rawMatrix[k] = [name, curRatio, cm[0][0], cm[0][1], cm[1][0], cm[1][1]]

                predictions = np.array(y_pred/eNN_iter)
                for zz in range(len(X_test)):
                    predCount.loc[test_index[zz], name] += 1
                    predScore.loc[test_index[zz], name] += predictions[zz]

                yhats = [model.predict(X_val) for model in members]
                yThats = [np.where(yh>=cutOffForPosLabel,1,0) for yh in yhats]     ## predict class by larger or less than 0.5 probability
                yThats = [yh.T[0] for yh in yThats]         ## transform the array
                y_pred = np.sum(yThats, axis=0)
                y_final = np.array([1 if ii > (eNN_cutoff * eNN_iter) else 0 for ii in y_pred])

                cm = confusion_matrix(y_val, y_final)
                valMatrix[k] = [name, curRatio, cm[0][0], cm[0][1], cm[1][0], cm[1][1]]

                k += 1

                tf.keras.backend.clear_session()

    for i in range(predScore.shape[1]):
        predScore.iloc[:,i] = predScore.iloc[:,i]/predCount.iloc[:,i]

    rawMatrix = pd.DataFrame(rawMatrix)
    valMatrix = pd.DataFrame(valMatrix)
    predScore = pd.DataFrame(predScore)

    return rawMatrix, valMatrix, predScore


def RFmodel(X, y, RFtable, featureList, predratio, featureDict, iter_n=100, seedGeneWeight=100, seed_value=18, dynamicRFmaxLeafNodes=False, RFmaxLeaftoLabelsRatio=3, RF_max_leaf_nodes=12, runRFoob=False):
    ## run random forest prediction model for n times, stratified shuffle spilt with certain ratio.
    sss = StratifiedShuffleSplit(n_splits=iter_n, test_size=predratio, random_state=seed_value)
    j = 0
    for train_index, pred_index in sss.split(X, y):
        j += 1
        X_train, X_pred = X[train_index], X[pred_index]
        y_train, y_pred = y[train_index], y[pred_index]

        # additional test line - set maximum nodes as (number of labels in y_train)/RFmaxLeaftoLabelsRatio
        if dynamicRFmaxLeafNodes:
            RF_max_leaf_nodes = math.floor(len(y_train[y_train==1])/RFmaxLeaftoLabelsRatio)

        # run random forest predictions
        rf = RandomForestClassifier(max_leaf_nodes=RF_max_leaf_nodes, class_weight={1: seedGeneWeight, 0: 1}, oob_score=runRFoob, random_state=seed_value)
        rf.fit(X_train, y_train)
        if runRFoob:
            oob_error = 1 - rf.oob_score_
        y_RF = rf.predict_proba(X)
        y_RFpredict = [round(i[1], 4) for i in y_RF]  # convert array to vector of y prediction value
        #y_RFpredict = [None if ii in pred_index else y_RFpredict[ii] for ii in range(len(y_RFpredict))]  # only keep predictions from prediction set
        y_RFpredict = [y_RFpredict[ii] if ii in pred_index else None for ii in range(len(y_RFpredict))]
        curRname = "pred_" + str(j)
        RFtable[curRname] = y_RFpredict

        # feature importances
        importances = list(rf.feature_importances_)
        # List of tuples with variable and importance
        feature_importances = [(feature, round(importance, 4)) for feature, importance in zip(featureList, importances)]
        # Sort the feature importances by most important first
        feature_importances = sorted(feature_importances, key=lambda x: x[1], reverse=True)

        # Record features have effect >0
        for k in range(len(feature_importances)):
            if feature_importances[k][1] > 0:
                try:
                    featureDict[feature_importances[k][0]][predratio].append(feature_importances[k][1])
                except KeyError:
                    featureDict[feature_importances[k][0]][predratio] = [feature_importances[k][1]]

    return RFtable, featureDict


def RF_pred(dataset, seedGeneList, predRatioStep=[0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8], iter_n=100, seedGeneWeight=100, seed_value=18, runRFoob=False,
            dynamicRFmaxLeafNodes=False, RFmaxLeaftoLabelsRatio=3, RF_max_leaf_nodes=12, bgnoisethreshold=0.1, SUFFIX='testRun_20210831', outPATH="./OUT/"):
    ## RF_max_leaf_nodes and dynamic RFmaxLeafNodes are mutually exclusive, if dynamicRFmaxLeafNodes set into True, will overide the fix RF_max_leaf_nodes
    ## if dynamicRFmaxLeafNodes is true, the max leaf nodes for RF is set with RFmaxLeaftoLabelsRatio.

    ## extract X for features and y for labels
    X, y = SGfn.extract_info(dataset, seedGeneList)
    X = SGfn.preprocess_feature(X, imputeData=True, scaleFeature=True)

    ## create feature dict
    featureDict = {}
    for f in dataset.columns.values[1:]:
        featureDict[f] = {}

    geneList = list(dataset["GeneSymbol"])
    finalPredTable = {'GeneSymbol':geneList}

    current_time = time.strftime("%H:%M:%S", time.localtime())
    print("%s - Start RF model training for %s." % (current_time, SUFFIX))

    ## run different test ratio step
    for i in range(len(predRatioStep)):
        RFtable = {'GeneSymbol':geneList}
        curRatio = "predwithTrainPredR_" + str(round(1 - predRatioStep[i], 1)) + '|' + str(predRatioStep[i])
        print(curRatio)

        RFtable, featureDict = RFmodel(X, y, RFtable=RFtable, featureList=dataset.columns[1:], predratio=predRatioStep[i], featureDict=featureDict, iter_n=iter_n, seedGeneWeight=seedGeneWeight, seed_value=seed_value, dynamicRFmaxLeafNodes=dynamicRFmaxLeafNodes, RFmaxLeaftoLabelsRatio=RFmaxLeaftoLabelsRatio, RF_max_leaf_nodes=RF_max_leaf_nodes, runRFoob=runRFoob)
        # replace prediction less than background noise with 0 (default background noise threshold=0.1)
        RFtable = pd.DataFrame(RFtable)
        RFtable.iloc[:, 1:] = RFtable.iloc[:, 1:].apply(lambda x: [y if (y >= bgnoisethreshold or math.isnan(y)) else 0 for y in x])
        curYpred = RFtable.mean(axis=1, skipna=True, numeric_only=True)
        finalPredTable[curRatio] = curYpred

    finalPredTable = pd.DataFrame(finalPredTable)
    finalPredTable.to_csv(outPATH + SUFFIX + "_finalPredictionScore.csv", index=False)
    writeFeatureImportance(featureDict, filename=outPATH+SUFFIX+"_feature_importance.csv")

    current_time = time.strftime("%H:%M:%S", time.localtime())
    print("%s - RF model training - Done!" % (current_time))
    return finalPredTable
