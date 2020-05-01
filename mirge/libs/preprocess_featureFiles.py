import os
import sys
import pandas as pd
import matplotlib
matplotlib.use('agg')
from sklearn.preprocessing import OneHotEncoder
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold, chi2, f_classif
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC
from sklearn.metrics import make_scorer, confusion_matrix, matthews_corrcoef, roc_auc_score, roc_curve, auc
from sklearn.pipeline import Pipeline
import matplotlib.pyplot as plt
from pathlib import Path
import time
from matplotlib.backends.backend_pdf import PdfPages
import joblib

def preprocess_featureFiles(outputdir2, files, infTmp, feature_namelist_file):
    # infTmp is unmapped_mirna_JH-29_Pros_SMC_vs_representative_seq_modified_selected_sorted_features_updated_stableClusterSeq_15.tsv
    readCountLimit = 10
    seqCountLimit = 3
    tmpName = str(Path(outputdir2)/(files + '_dataset_15.csv'))
    #tmpName = '_'.join(os.path.basename(infTmp).split('_')[:-9])+'_dataset_'+os.path.basename(infTmp).split('_')[-1].split('.')[0]
    #with open(os.path.join(os.path.dirname(infTmp), tmpName+'.csv'), 'w') as outf:
    with open(tmpName, 'w') as outf:
        with open(infTmp, 'r') as inf:
            line = inf.readline()
        content = line.strip().split('\t')
        retainedIndexList = [0, 1, 5] + list(range(10, len(content)))
        retainedColNameList = [content[i] for i in range(len(content)) if i in retainedIndexList]
        readCountLimitLabel = content.index('readCountSum')
        seqCountLimitLabel = content.index('seqCount')

        with open(infTmp, 'r') as inf:
            line = inf.readline()
            content = line.strip().split('\t')
            outf.write(','.join([content[i] for i in range(len(content)) if i in retainedIndexList]))
            outf.write('\n')
            line = inf.readline()
            while line != '':
                content = line.strip().split('\t')
                if int(content[readCountLimitLabel]) >= readCountLimit and int(content[seqCountLimitLabel]) >= seqCountLimit and ('N' not in content[seqCountLimitLabel:]):
                    outf.write(','.join([content[i] for i in range(len(content)) if i in retainedIndexList]))
                    outf.write('\n')
                line = inf.readline()
            line = inf.readline()

    tmpName2 = str(Path(outputdir2)/(files + '_dataset_15_refined_tmp.csv'))
    outf = open(tmpName2, 'w')
    with open(tmpName, 'r') as inf:
        line = inf.readline()
        outf.write(line)
        content = line.strip().split(',')
        pruned_precusor_strLabel = content.index('pruned_precusor_str')
        line = inf.readline()
        while line != '':
            content = line.strip().split(',')
            if content[pruned_precusor_strLabel] != 'None':
                outf.write(line)
            else:
                pass
            line = inf.readline()
    outf.close()

    nameList = []
    with open(feature_namelist_file, 'r') as inf:
        for line in inf:
            if line.strip() not in nameList:
                nameList.append(line.strip())

    tmpName3 = str(Path(outputdir2)/(files + '_dataset_15_refined_features.csv'))
    outf = open(tmpName3, 'w')
    with open(tmpName2, 'r') as inf:
        retainedIndexList = []
        line = inf.readline()
        content = line.strip().split(',')
        reatinedContent = []
        for i, name in enumerate(content):
            if name in nameList:
                retainedIndexList.append(i)
                reatinedContent.append(name)
        outf.write(','.join(reatinedContent))
        outf.write('\n')
        line = inf.readline()
        while line != '':
            content = line.strip().split(',')
            reatinedContent = []
            for j in range(len(content)):
                if j in retainedIndexList:
                    reatinedContent.append(content[j])
            outf.write(','.join(reatinedContent))
            outf.write('\n')
            line = inf.readline()
    outf.close()


def model_predict(outputdir2, files, modelFile):
    inputFile = str(Path(outputdir2)/(files+"_dataset_15_refined_features.csv"))
    cutoff_picked = 0.8
    sc, clf, selectFeatureNameList = joblib.load(modelFile)
    data = pd.read_csv(inputFile)
    data_raw = data.iloc[:,:].values
    featureList_raw = data.columns.values.tolist()[:3]
    featureList_raw.insert(0, 'probability')
    featureList_raw.insert(0, 'predictedValue')
    # Encode categorical features using one-hot-coding after deleting the 2nd and 3rd column in data
    data.drop(data.columns[[0,1,2]], axis=1, inplace=True)
    data = pd.get_dummies(data)
    # Check the dummy varibles, such as 'pair_state_Yes', 'pair_state_No', etc.
    # If some dummy varibles are missing, they will be added with value of 0 across all of the samples.
    totalfeatureListTmp = data.columns.values.tolist()
    missedFeatureList = []
    for item in selectFeatureNameList:
        if item[2] not in totalfeatureListTmp:
            missedFeatureList.append(item[2])
    for item in missedFeatureList:
        data[item] = np.array([0]*(data.shape[0]))
    # Subselect data with the selected features and reorder the feature sequence according to the feature sequence in selectFeatureNameList
    totalfeatureList = data.columns.values.tolist()
    subIndexList = [totalfeatureList.index(item[2]) for item in selectFeatureNameList]
    
    data_x = data.iloc[:,subIndexList].values
    x_test_std_new = sc.transform(data_x)
    #featureLabel_new = subIndexList
    x_test_std_new2 = x_test_std_new
    predictedValue = clf.predict(x_test_std_new2)
    probability = [max(value) for value in clf.predict_proba(x_test_std_new2)]
    # Output the predicted result of the test data set into a file.
    array_predicted = np.column_stack((predictedValue, probability, data_raw[:,:3]))
    data_predicted_test = pd.DataFrame(array_predicted, columns=featureList_raw)
    data_predicted_test.to_csv(str(Path(outputdir2)/('predicted_'+files+'_detailed.csv')), index=False)
    outf = open(str(Path(outputdir2)/(files+'_novel_miRNAs_miRge2.0.csv')),'w')
    contentList_tmp = []
    with open(str(Path(outputdir2)/('predicted_'+files+'_detailed.csv')), 'r') as inf:
        line = inf.readline()
        outf.write(line)
        line = inf.readline()
        while line != '':
            content = line.strip().split(',')
            if content[0] == '1' and float(content[1]) >= cutoff_picked:
                contentList_tmp.append([float(content[1]), line])
            line = inf.readline()
    contentList_tmp.sort(reverse=True)
    for item in contentList_tmp:
        outf.write(item[1])
    outf.close()
