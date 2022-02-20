# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 23:51:11 2020

@author: Amaya

"""
import pandas as pd
import numpy as np
import pyHSICLasso as hsic

import os

from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPClassifier

from sklearn.svm import SVC

import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score

# %%

X_files = ['df_cd4.csv', 'df_cd8.csv', 'df_cd14.csv','df_pbmc.csv', 'df_wb.csv' ]
y_files = ['y_cd43class.csv', 'y_cd83class.csv', 'y_cd143class.csv','y_pbmc3class.csv', 'y_wb3class.csv' ]
y_files = ['cd42class.csv', 'cd82class.csv', 'cd142class.csv','pbmc2class.csv', 'wb2class.csv' ]
y_files = ['cd4orig.csv', 'cd8orig.csv', 'cd14orig.csv','pbmcorig.csv', 'wborig.csv' ]
set_names = ['cd4', 'cd8', 'cd14', 'pbmc', 'wb']

n_splits = 5

for name_idx, files in enumerate(zip(X_files, y_files)):

    dfx = pd.read_csv(files[0], )
    dfy = pd.read_csv(files[1], )
    df_merge = pd.merge(dfy, dfx, left_index=True, right_index=True)
    df_merge = df_merge.rename(columns={"outcome": "class"})
    df_merge.drop(df_merge.filter(regex="Unname"),axis=1, inplace=True)
    dfList = [(df_merge[df_merge.index % n_splits != i], df_merge[df_merge.index % n_splits == i]) for i in range(n_splits)]
    
    for idx, (dftrain,dftest) in enumerate(dfList):
        dftrain.to_csv(set_names[name_idx] + '_train_' + str(idx) + '.csv', index = False )
        dftest.to_csv(set_names[name_idx] + '_test_' + str(idx) + '.csv', index = False)
    

#%% feature Selection with MaxRel, mRMR and  HSIC Lasso - runs for a couple of hours at least

n_feat = 120

celltypes = ['cd4', 'cd8', 'cd14', 'pbmc', 'wb']
feat_sel_types = ['maxRel', 'mRMR', 'HSIC']
n_selected = [1,2,3,4,5,7,10,12,15,17,20,22,25,30,32,35,40,50,60,70,80, 100,120 ]

result_dict = {}
selected_genes = {}
for celltype in celltypes:
    
    input_datasets = [celltype + '_train_' + str(i) + '.csv' for i in range(n_splits)]
    test_datasets = [celltype + '_test_' + str(i) + '.csv' for i in range(n_splits)]
    
    for key in feat_sel_types:
        selected_genes[key] = {} 
      
    results = {}
    
    for ds in input_datasets:
        
        command = 'mrmr_win32 -i ' + ds + ' -t 0.5 -n '+ str(n_feat)+ ' -m MID -v 18600 > mrmr_' + ds
        os.system('cmd /c'  + command ) 
        f = open('mrmr_' + ds , 'r')
        lines = f.readlines()
        selected_genes['maxRel'][ds] = [line.split("\t")[2].strip() for line in lines[6:(n_feat + 6)]]
        selected_genes['mRMR'][ds]  = [line.split("\t")[2].strip() for line in lines[(n_feat + 9):(2 * n_feat + 9)]]
        
        hl = hsic.HSICLasso()
        hl.input(ds)
        hl.classification(n_feat)
        selected_genes['HSIC'][ds] = [hl.featname[i] for i in hl.get_index()]
        
        
        
    result_dict[celltype] = {'mRMR'  : selected_genes['mRMR'],
                             'maxRel': selected_genes['maxRel'],
                             'HSIC'  : selected_genes['HSIC'] }
    
    
#%% Classification with Logistic, SVM and Multi Layer Perceptron - runs for a couple of hours at least

#more could be added here
    
classifiers = {'SVM': SVC(kernel="rbf"),
               'Logistic': LogisticRegression(),
               'MLPC': MLPClassifier(max_iter=2000)}    

classification = {}

outcomes = {}

    
for celltype in celltypes:
    
    print('processing', celltype)
    classification[celltype] = {}
    input_datasets = [celltype + '_train_' + str(i) + '.csv' for i in range(n_splits)]
    test_datasets = [celltype + '_test_' + str(i) + '.csv' for i in range(n_splits)]
    
    
    for method in feat_sel_types:
        
        classification[celltype][method] = {}
        print('processing', celltype,method )    
        
        for inputds, testds in zip(input_datasets, test_datasets):
            
            print('processing', celltype,method, inputds  )  
            classification[celltype][method][inputds] = {}
            df = pd.read_csv(inputds)
            y = df['class'].to_numpy()
            
            if celltype not in outcomes:
                
                outcomes[celltype] = y
            df_test = pd.read_csv(testds)
            y_test = df_test['class'].to_numpy()
            good_genes = result_dict[celltype][method][inputds]
            temp_result = {clf:[] for clf in classifiers.keys()}
            
            for n_feat in n_selected:
                
                filtered_df = df[good_genes[0:n_feat]]
                filtered_df_test = df_test[good_genes[0:n_feat]]
                X = filtered_df.to_numpy()
                X = StandardScaler().fit_transform(X)     
                X_test = filtered_df_test.to_numpy()
                X_test = StandardScaler().fit_transform(X_test)     
                
                for clf in classifiers.keys():
                    
                    classifiers[clf].fit(X, y)
                    y_predict = classifiers[clf].predict(X_test)
                    acc = accuracy_score(y_test, y_predict)
                    temp_result[clf].append(acc)
                    print(n_feat, acc, clf)
                    
                
 
            classification[celltype][method][inputds]  = temp_result    
       
#%% Build Plots

avgResultsDict = {}
savefigs =  True

for celltype in celltypes:
    
    max_freq =  max(np.unique(outcomes[celltype], return_counts=True)[1])/len(y)
    avgResultsDict[celltype]={}
    
    for clf in classifiers.keys():
        
        avgResultsDict[celltype][clf]={}
        plt.figure()
        
        for method in feat_sel_types:
            
            list_of_CVresults = [classification[celltype][method][ds][clf] for ds in [celltype + '_train_' + str(i) + '.csv' for i in range(n_splits)]]
            avgResultsDict[celltype][clf][method] = np.average(list_of_CVresults, axis=0)
            plt.plot(n_selected,avgResultsDict[celltype][clf][method], label=method)
        plt.axhline(y=max_freq, color='r', linestyle='dashed', label = 'Mode Frequency')
        plt.legend(fontsize=8)
        plt.xlabel('# of Selected Features')
        plt.ylabel('Average CV Accuracy')
        plt.title('Prediction accuracy - ' + clf + ' classifier - cell type ' + celltype.upper())
       
        plt.ylim(0.,0.6)
        #plt.show()
        if savefigs:
            plt.savefig(celltype + '_' + clf + '.pdf')
    


            