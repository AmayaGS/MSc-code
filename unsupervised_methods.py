# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 19:53:58 2020

@author: Amaya
"""
# %%

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#plt.style.use("ggplot")

%matplotlib inline

import seaborn as sns

#import pyHSICLasso as hsic

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

import umap
from sklearn.cluster import KMeans
from sklearn import metrics
from sklearn.cluster import AgglomerativeClustering
from sklearn import linear_model
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier, GradientBoostingClassifier, ExtraTreesClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn.model_selection import train_test_split, GridSearchCV, StratifiedKFold
from sklearn.metrics import roc_curve, precision_recall_curve, auc, make_scorer, recall_score, accuracy_score, precision_score, confusion_matrix
from sklearn.metrics import classification_report
from sklearn.metrics import mean_squared_error, r2_score
import hdbscan
import sklearn.cluster as cluster
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score
from scipy.stats import gaussian_kde
from pygam import LinearGAM


# %%

# loading the baseline gene expression data
data = pd.read_table('baseline_expression_all.txt', index_col= 0)

# %%

# loading all non gene features and outcome labels - original labels, DAS28 2 class and DAS28 3 class labels
outcomes = pd.read_table('baseline_expr_all_pheno.csv', index_col = 0, sep= ',')

# %%

# concatenating whole data set and phenotype (cell type) information to run tSNE on it looking and check if each cell type clusters together homogeneously. 
data_transp = data.T
cell_type = outcomes['pheno']
df = pd.concat([data_transp, cell_type], axis = 1, join = 'inner', sort=False, verify_integrity= True)
df = df.dropna()
df_labels = df['pheno']
X_df = df.drop('pheno', axis=1)


# %%

y_class = outcomes[['ID', 'pheno', 'DAS28 6M 3 class outcome']]
y_class = y_class.dropna()
y_class = y_class.rename(columns={'DAS28 6M 3 class outcome': 'outcome'})


# %%
# definition for preprocessing the dataset into data and labels for each cell type

def data_preproc(data, outcome, celltype):

    df = data.T

    outcome_df = pd.concat([df, outcome], axis = 1, join = 'inner', sort=False, verify_integrity= True)
    outcome_df = outcome_df.dropna()

    outcome_df['outcome'] = pd.factorize(outcome_df['outcome'])[0]

    outcome_df = outcome_df[outcome_df['pheno'] == celltype]

    X = outcome_df.drop(['ID','pheno', 'outcome'], axis=1)
 
    y_labels = pd.DataFrame(outcome_df['outcome'])

    return X, y_labels

# %%    

# subsetting the dataset into different cell types

df_cd14, y_cd14 = data_preproc(data, y_class, 'CD14')
df_cd4, y_cd4 = data_preproc(data, y_class, 'CD4')
df_cd8, y_cd8 = data_preproc(data, y_class, 'CD8')
df_pbmc, y_pbmc = data_preproc(data, y_class, 'PBMC')
df_wb, y_wb = data_preproc(data, y_class, 'WB')

# %%

labels = y_labels.reset_index(drop=True)
finalDf = pd.concat([principalDf, labels], axis = 1)

# %%

X_pca = df_cd8
labels = y_cd8

# To getter a better understanding of interaction of the dimensions
# plot the first three PCA dimensions
fig = plt.figure(1, figsize=(8, 6))
ax = Axes3D(fig, elev=-150, azim=210)
pca = PCA(n_components= 20)
X = pca.fit_transform(X_pca)

#y_pca_2C_WB = np.choose(y_pca_2C_WB, [0, 1]).astype(np.float)
sc = ax.scatter(np.reshape(X[:, 0], -1), np.reshape(X[:, 1], -1), np.reshape(X[:, 2], -1), c = np.reshape(labels, -1),  cmap="brg", edgecolor='k', s=40)

#ax.set_title()
ax.set_xlabel("PC1")
ax.w_xaxis.set_ticklabels([])
ax.set_ylabel("PC2")
ax.w_yaxis.set_ticklabels([])
ax.set_zlabel("PC3")
ax.w_zaxis.set_ticklabels([])

colors = [sc.cmap(sc.norm(i)) for i in [0, 1, 2]]
custom_lines = [plt.Line2D([],[], ls="", marker='.', 
                mec='k', mfc=c, mew=.1, ms=20) for c in colors]

plt.show()


# %%

# saving as csv files for later use

df_cd14.to_csv('df_cd14.csv', sep=',', index = True)
y_cd14.to_csv('y_cd14.csv', sep=',', index = True)
df_cd4.to_csv('df_cd4.csv', sep=',', index = True)
y_cd4.to_csv('y_cd4.csv', sep=',', index = True)
df_cd8.to_csv('df_cd8.csv', sep=',', index = True)
y_cd8.to_csv('y_cd8.csv', sep=',', index = True)
df_pbmc.to_csv('df_pbmc.csv', sep=',', index = True)
y_pbmc.to_csv('y_pbmc.csv', sep=',', index = True)
df_wb.to_csv('df_wb.csv', sep=',', index = True)
y_wb.to_csv('y_wb.csv', sep=',', index = True)


# %%

#Calculate estimates of variance, coefficient of variation

subset = [X_df, df_cd4, df_cd8, df_cd14, df_pbmc, df_wb]
sub_names = ['Whole dataset', 'CD4', 'CD8', 'CD14', 'PBMC', 'WB']

for celltype, name in zip(subset, sub_names):
    
    means = celltype.mean(axis = 0)
    var = celltype.var(axis = 0)
    means2 = means.pow(2, axis = 0)
    CV = var.div(means2)
    
    log_means = np.log2(means)
    log_cv2 = np.log2(CV )
    
    # Calculate the point density # https://stackoverflow.com/a/20107592
    xy = np.vstack([log_means,log_cv2])
    z = gaussian_kde(xy)(xy)
    
    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    log_means, log_cv2, z = log_means[idx], log_cv2[idx], z[idx]

    # fitting GAM model to the data
    gam = LinearGAM(n_splines=25).gridsearch(log_means.values.reshape(18561,1), log_cv2.values.reshape(18561,1))

    # plotting
    fig = plt.figure(figsize = (12,8))
    ax = fig.add_subplot(1,1,1)
    
    ax.plot(XX, gam.predict(XX), 'r--')
    ax.plot(XX, gam.confidence_intervals(XX, width=.95), color='c', ls='--')
    
    ax.set_xlabel('log(mean) of gene expression', fontsize = 14)
    ax.set_ylabel("log(CV)", fontsize = 14)
    
    plt.scatter(log_means, log_cv2, c=z, s=10)
    plt.title('' + str(name), fontsize = 15)
    
    plt.savefig('mean_cv_' + str(name) +'.pdf')

# %%

# Import High variance genes from R and run corr plot
df_cd4_hv = pd.read_table('cd4_hv_expr.csv', index_col= 0, sep= ',')
df_cd8_hv = pd.read_table('cd8_hv_expr.csv', index_col= 0, sep= ',')
df_cd14_hv = pd.read_table('cd14_hv_expr.csv', index_col= 0, sep= ',')
df_pbmc_hv = pd.read_table('pbmc_hv_expr.csv', index_col= 0, sep= ',')
df_wb_hv = pd.read_table('wb_hv_expr.csv', index_col= 0, sep= ',')


# %%

# correlation plots for each subset
subset_hv = [X_df, df_cd4_hv, df_cd8_hv, df_cd14_hv, df_pbmc_hv, df_wb_hv]
sub_names = ['Dataset', 'CD4', 'CD8', 'CD14', 'PBMC', 'WB']

for celltype, name in zip(subset_hv, sub_names):
    
    sum_corr = celltype.corr().abs().sum().sort_values(ascending=True).index.values
    df = celltype[sum_corr]
    
    corr = df.corr()
    
    fig = plt.figure()
    ax = fig.add_subplot()
    ax = sns.heatmap(
        corr, 
        vmin=-1, vmax=1, center=0,
        cmap= "vlag",
        square=True,
        xticklabels=False,
        yticklabels=False, rasterized=True
    )
    
    plt.title('' + str(name), fontsize = 10)
    
    #fig.savefig('corr_' + str(name) + '.pdf')

# %%

# correlation plots for a random dataset

import random

random.seed(32)
random_genes = np.random.random([200, 800])

df_rand_genes = pd.DataFrame(random_genes)

sum_corr = df_rand_genes.corr().abs().sum().sort_values(ascending=True).index.values
df = df_rand_genes[sum_corr]

corr = df.corr()

fig = plt.figure()
ax = fig.add_subplot()
ax = sns.heatmap(
    corr, 
    vmin=-1, vmax=1, center=0,
    cmap= "vlag",
    square=True,
    xticklabels=False,
    yticklabels=False, rasterized=True
)

plt.title('Random data', fontsize = 10)

fig.savefig('corr_random.pdf')

# %%
# First I want to see if the whole data set clusters in accordance with the different cell types, hence indicating that each cell type gives off its on distinctive signal. I'll use tSNE for this as its a non linear transform and hence more sensitive to complex high dimensional structure. 

# t_SNE for different perplexities levels - going from local structure to global structure. Takes some time to run. 8min 20s on my computer. No visible batch effect. 

%%time

perplexities = [4, 15, 30, 50]

for p in perplexities:
    
    tsne = TSNE(n_components=2, verbose=1, perplexity= p, n_iter=5000, learning_rate=100)
    tsne_results = tsne.fit_transform(ss_df)
    
    df_tSNE = pd.DataFrame(tsne_results, columns = ['tSNE1', 'tSNE2'])
    labels = df_labels.reset_index(drop=True)
    final_df_tSNE = pd.concat([df_tSNE, labels], axis = 1)
    
    fig = plt.figure(figsize = (8,8))
    ax = fig.add_subplot(1,1,1) 
    ax.set_xlabel('tSNE1', fontsize = 15)
    ax.set_ylabel('tSNE2', fontsize = 15)
    ax.set_title('tSNE on whole dataset, all celltypes - perplexity = %s' % p, fontsize = 20)
    targets = ['CD4', 'CD14', 'CD8', 'PBMC', 'WB']
    colors = ['r', 'g', 'c', 'm', 'b']
    
    for target, color in zip(targets, colors):
        indicesToKeep = final_df_tSNE['pheno'] == target
        ax.scatter(final_df_tSNE.loc[indicesToKeep, 'tSNE1']
                   , final_df_tSNE.loc[indicesToKeep, 'tSNE2']
                   , c = color
                   , s = 10)
        
    ax.legend(targets)
    ax.grid()
    # ADD savefig here



# %%
# UMAP

# Next I want to see if I replicate these results using UMAP, which also appliea a non linear transform to the data reducing it to a low dimensional manifold, seeking to preserve both local and global structure.

# %%

# transforming the whole dataset
# the parameter n_neighbors controls the number of neighbours taken into account to perform the reduction. 15 is the default, which tries to balance local/global structure, whilst 30 tends more to the global structure. For use of UMAP for clustering with HDBSCAN 30 is the recommended value. 

%%time
embedding = umap.UMAP(n_neighbors= 30, metric='cosine', transform_seed = 34, random_state = 44, n_components = 2).fit_transform(X_df)

# %%

# passsing cell type to numerical labels
cell_code = pd.factorize(cell_type)[0]
classes = ['CD14', 'CD8', 'PBMC', 'CD4', 'WB']

# %%

# plotting the result of the UMAP reduction
fig, ax = plt.subplots(1, figsize=(14, 10))
plt.scatter(*embedding.T, s=10, c= cell_code, cmap="brg", alpha=1.0)
plt.setp(ax, xticks=[], yticks=[])
cbar = plt.colorbar(boundaries=np.arange(6)-0.5)
cbar.set_ticks(np.arange(5))
cbar.set_ticklabels(classes)
plt.title('UMAP reduction on the whole dataset')

plt.savefig("umap_wholedata_30.pdf") 


# %%

# I now want to reduce the feature set for each cell type, as the tSNE and UMAP analysis showed that each cell type has a clear signal. The result of the feature reduction can then be clustered using the hbscan algorithm to obtain labels for each sample. I do this in the hope of identifying patient sub-phenotypes and usable labels. 

# NOTE: all the predefined distance metrics where tried, but I only plot Manhattan, Cosine, Minkowski, Euclidean and Chebyshev distances which showed signs of potential clusters. Chebyshev distance showed  very clear clusters in all cell type subsets compared to other metrics, which I find strange and I'm not sure what it means


subset = [df_cd4, df_cd8, df_cd14, df_pbmc, df_wb]
sub_names = ['CD4', 'CD8', 'CD14', 'PBMC', 'WB']
distance = ['euclidean', 'cosine', 'manhattan', 'minkowski', 'chebyshev']

for celltype, name in zip(subset, sub_names):
    
    for metric in distance:
        
            embedding = umap.UMAP(n_neighbors= 18, metric=metric, transform_seed = 34, random_state = 44, n_components = 50, min_dist = 0.0).fit_transform(celltype)
            hdbscan_labels = hdbscan.HDBSCAN(min_cluster_size=4).fit_predict(embedding)
            
            np.savetxt('hdbscan_labels' + str(name) + str(metric) + '.txt', hdbscan_labels, delimiter = ',')
    
            clustered = (hdbscan_labels >= 0)
            plt.figure()
            plt.scatter(embedding[clustered, 0],
                    embedding[clustered, 1],
                    c= hdbscan_labels[clustered],
                    s= 15,
                    cmap='brg')
                
            plt.title('%s' % name + ' ' + metric, fontsize = 15)
            
            plt.savefig('umap_hdbscan_' + str(name) + str(metric) +'.pdf')


# %%

# Aplying UMAP and HDBSCAN in same conditions to random uniform data. 

import random

random.seed(32)
random_data = np.random.random([200, 18561])

rand_embedding = umap.UMAP(n_neighbors= 18, metric='chebyshev', transform_seed = 34, random_state = 44, n_components = 50, min_dist = 0.0).fit_transform(random_data)
rand_hdbscan_labels = hdbscan.HDBSCAN(min_cluster_size=4).fit_predict(rand_embedding)

rand_clustered = (rand_hdbscan_labels >= 0)
plt.scatter(rand_embedding[rand_clustered, 0],
            rand_embedding[rand_clustered, 1],
            c= rand_hdbscan_labels[rand_clustered],
            s=15,
            cmap='brg')

plt.title('Random data')

plt.savefig("umap_random_data.pdf")

# %%

# We'll now divide the dataset in training and testing set to see if the testing set can recover the labels based on the embedded feature reduction. 

# %%


from sklearn.model_selection import KFold
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score


subset = [df_cd4, df_cd8, df_cd14, df_pbmc, df_wb]
sub_names = ['CD4', 'CD8', 'CD14', 'PBMC', 'WB']
distances = ['chebyshev', 'cosine', 'manhattan', 'minkowski', 'euclidean']

results = pd.DataFrame(columns=['celltype', 'distance', 'n_val', 'classifier', 'accuracy' ])

for celltype, name in zip(subset, sub_names):

    for metric in distances:
        
        ss = KFold(n_splits= 5, shuffle=False)
        ss.get_n_splits(celltype)
        temp_results = {}
        
        for idx, (train, test) in enumerate(ss.split(celltype)):
            
            df = celltype
            
            trans = umap.UMAP(n_neighbors= 30, metric= metric, transform_seed = 34, random_state = 44, n_components = 50, min_dist = 0.0).fit(df.iloc[train])
            test_embedding = trans.transform(df.iloc[test])
            
            hdbscan_obj =  hdbscan.HDBSCAN(min_cluster_size=4, prediction_data=True)            
            train_hdbscan_labels = hdbscan_obj.fit_predict(trans.embedding_)
            test_hdbscan_labels = hdbscan.prediction.approximate_predict(hdbscan_obj, test_embedding)[0]
                         
            models = ["Decision Tree", "rbf SVM C11", "Random Forest", "AdaBoost"]
            
            classifiers = [DecisionTreeClassifier(max_depth=100),
               #KNeighborsClassifier(3),
               SVC(kernel="rbf", C=10),
               RandomForestClassifier(max_depth=25, n_estimators=15, max_features=10),
               AdaBoostClassifier(),
               #GaussianNB()
              ]
            
            for model, clf in zip(models, classifiers):
            
                   clf.fit(df.iloc[train], train_hdbscan_labels)
                   y_predict = clf.predict(df.iloc[test])
                   
                   # obtain the classification accuracy on the test data
                   acc = accuracy_score(test_hdbscan_labels, y_predict)
                   print(model, name, metric)    
                   print ('Accuracy:', acc)
                   results.loc[len(results)] = [name, metric, idx, model, acc]
       
# ADD PLOTS to visualize results?

# %%

results.to_csv('notransform_cluster_classification_results.csv', sep=',', index = True)

# %%

non_gene_features = outcomes[['ID', 'pheno', 'FATIQUE.0M', 'PAIN.0M', 'TOTAL.SWOLLEN.0M', 'TOTAL.TENDER.0M']]

# %%


embedding = umap.UMAP(n_neighbors= 18, metric='manhattan', transform_seed = 34, random_state = 44, n_components = 50, min_dist = 0.0).fit_transform(df_cd4)
hdbscan_labels = hdbscan.HDBSCAN(min_cluster_size=4).fit_predict(embedding)

np.savetxt('hdbscan_labels' + str(name) + str(metric) + '.txt', hdbscan_labels, delimiter = ',')



hdbscan_labels_cd4 = pd.read_table('hdbscan_labelsCD4chebyshev.txt', delimiter =',', header=None)

# %%

feat_ng = non_gene_features[non_gene_features['pheno'] == 'CD4']

outcome = pd.DataFrame(hdbscan_labels)
#outcome = hdbscan_labels.rename(columns={0: "outcome"})
outcome.reset_index(drop=True, inplace=True)
feat_ng.reset_index(drop=True, inplace=True)

# %%

outcome_df = pd.concat([feat_ng, outcome], axis = 1, join = 'outer')
outcome_df = outcome_df.dropna()

# %%

X = outcome_df.drop(['ID','pheno', 0], axis=1)
   
y = pd.DataFrame(outcome_df[0])

y = y.to_numpy()

y = y.ravel()



# %%

X = StandardScaler().fit_transform(X)

# %%


names = ["Nearest Neighbors", "rbf SVM C11",
        "Decision Tree", "Random Forest",
        "Logistic","Neural Net"]

classifiers = [
   KNeighborsClassifier(3),
   SVC(kernel="rbf", C=10),
   DecisionTreeClassifier(max_depth=50),
   RandomForestClassifier(max_depth=50, n_estimators=8),
   LogisticRegression(max_iter=10000),
   MLPClassifier(alpha=1, max_iter=10000, hidden_layer_sizes=(10,10))
   ]

n_validations = 7

'X es un array con todos los datos. y un array con todos los outcomes'

for name, clf in zip(names, classifiers):
       print(name)
       print('cross validation Score:', cross_val_score(clf, X, y, cv=n_validations).mean())
       print()

# NO DA

# %%

df_cd4_hv.reset_index(drop=True, inplace=True)
#hdbscan_labels.reset_index(drop=True, inplace=True)

clusters = pd.concat([pd.DataFrame(hdbscan_labels), df_cd4_hv], axis = 1, join = 'outer')

#cluster0 = clusters[clusters[0] == -1]
#cluster1 = clusters[clusters[0] == 0]
#cluster2 = clusters[clusters[0] == 1]

clusters = clusters.set_index(0)

sns_plot = sns.clustermap(clusters.T, method='average', cmap='seismic', metric = 'cityblock', row_cluster = True, col_cluster = True, z_score = 0)

sns_plot.savefig('clustermap_manhattan_cd4.pdf')

# %%
df_cd4_hv.reset_index(drop=True, inplace=True)
#hdbscan_labels.reset_index(drop=True, inplace=True)

clusters = pd.concat([pd.DataFrame(hdbscan_labels), pd.DataFrame(embedding)], axis = 1, join = 'outer')


cluster0 = clusters[clusters[0] == -1]
cluster1 = clusters[clusters[0] == 0]
cluster2 = clusters[clusters[0] == 1]

clusters = clusters.set_index(0)

sns_plot = sns.clustermap(clusters.T, method='complete', cmap='seismic', metric = 'cityblock', row_cluster = True, col_cluster = True, z_score = 0)

sns_plot.savefig('clustermap_manhattan_cd4.pdf')

binary = {0: "navy", 1: "yellow"}
cluster_df[cat_columns].apply(lambda x: x.map(binary))

redpal = sns.light_palette('red', len(cluster_df.ESR_mm_hr.unique()))
lut1 = dict(zip(cluster_df.ESR_mm_hr.unique(), redpal))
col_colors = cluster_df.ESR_mm_hr.map(lut1)