# -*- coding: utf-8 -*-
"""
Gene expression analysis

Created on Fri Sep 30 08:32:11 2022

@author: Dimitrios Samaras
"""

#Import libraries needed
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from scipy import stats as st
import statistics
import math
from statsmodels.stats.multitest import multipletests
from bioinfokit import analys, visuz
import seaborn

###############################################################################
###1. Standardize the values since gene expression values may widely vary 

#Upload data
df = pd.read_csv("GeneExpressionData.txt",  sep='\t', index_col=False)

#By creating a file with statistics summary of the first markers you can notice their difference in expression
dfsubset = df.iloc[:, 0:9]
statist = dfsubset.describe()

#Save the group variable before excluding it from the standardization process
groups = df["group"]
df1 = df.loc[:, df.columns!='group']

#Standardize data with StandardScaler()
scaler = StandardScaler()
standardized_data = scaler.fit_transform(df1)

#Save the standardize data and retrieve the names of the genes from the original dataset
dfcolnames = df1.columns
df3 = pd.DataFrame(standardized_data)
df3.columns = dfcolnames
df3["group"] = groups

###############################################################################
###2. PCA analysis

#calculate the principal components with sklearn
pca = PCA(n_components=2)
principalComponents = pca.fit_transform(df3)
principalDf = pd.DataFrame(data = principalComponents, columns = ['PC1', 'PC2'])
finalDf = pd.concat([principalDf, df["group"]], axis = 1)

#Visualize PCA
plt.scatter(x = finalDf["PC1"], y = finalDf["PC2"], c =finalDf["group"])
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.title('PCA analysis')
plt.colorbar()
plt.savefig("PCA.png")

###############################################################################
###3. T-test to find the markers that constitute a statistically important difference between the 2 groups of people 

#Subset the original and the standardized datasets
group0 = df3[df3["group"] == 0]
group1 = df3[df3["group"] == 1]

group0_unstand = df[df["group"] == 0]
group1_unstand = df[df["group"] == 1]

#A loop to retrieve the p-values from the t-tests and the log2fc needed for the volcano plot, which is the next step
pvalue = []
log2fc = []
for i in dfcolnames:
    
    #Part of calculating p-values
    groupzero = group0[i]
    groupone = group1[i]
    
    results = st.ttest_ind(a=groupzero, b=groupone, equal_var=False)
    pvalue.append(results.pvalue)
    

    #Part of calculating log2fc values
    groupzero_unstand = group0_unstand[i]
    groupone_unstand = group1_unstand[i]
    
    mean0 = statistics.mean(groupzero_unstand)
    mean1 = statistics.mean(groupone_unstand)
    logvalue = math.log(mean1/mean0, 2)
    log2fc.append(logvalue)

#Create a dataset with markers' names, p-values and log2fc values
Pvalue = pd.DataFrame(data = pvalue, index = dfcolnames)
Log2fc = pd.DataFrame(data = log2fc, index = dfcolnames)
Pvalue.columns = ["pvalue"]
Log2fc.columns = ["log2fc"]
stats_results = Pvalue.join(Log2fc)

###############################################################################
###4. P-value adjustment 

#First method of adjustment
#Bonferroni correction by adjusting the 5% error rate
a = 0.05/len(stats_results["pvalue"])

#Second method of adjustment
#Adjusting p-value by using the bonferroni method from multipletests
p_adj = multipletests(pvalue, alpha= 0.05, method = "bonferroni")
padj = p_adj[1]
padj = pd.DataFrame(padj)
padj.index = Pvalue.index
padj.columns = ["p_adj"]
stats_results2 = pd.merge(stats_results, padj, left_index=True, right_index=True)

#Third way of p-value adjustment 
#FDR correction method
pvalue.sort()
rank = 1
len_pvalue = len(pvalue)
adjpvalue = []

for i in pvalue:
    fdr_pvalue = i * len_pvalue/rank
    rank += 1
    adjpvalue.append(fdr_pvalue)
adjpvalue = pd.DataFrame(adjpvalue)
adjpvalue.columns = ["adj_pvalue"]


pd.options.display.float_format = '{:.10f}'.format
stats_results2.head()

stats_results = stats_results.sort_values(by=["pvalue"])
adjpvalue.index = stats_results.index
stats_results = stats_results.join(adjpvalue) 

#Important markers based on adjusted p-values
only_important_markers_initial = stats_results[stats_results["pvalue"] < 0.05]

#Importnat markers based on adjusted values by FDR
only_important_markers_fdr = stats_results[stats_results["adj_pvalue"] < 0.05]

#Important values based on different error level calculated by Bonferroni correction
only_important_markers_a = stats_results[stats_results["pvalue"] < a]

#Important markers based on Bonferroni adjusted values by statsmodels
only_important_markers_bonf = stats_results2[stats_results2["p_adj"] < 0.05]

###############################################################################
###5. Volcano plot

visuz.GeneExpression.volcano(df=stats_results2, lfc='log2fc', pv='p_adj', lfc_thr = (0,0))

###############################################################################
###6. Heatmap

#Subset the dataset wit the standardized values of gene expression 
#so that it only has the genes with low p-values adjusted by Bonferroni method 
subset = df3.loc[:,list(only_important_markers_bonf.index)]

#Create and save the heatmap
seaborn.heatmap(subset)
plt.savefig("heatmap.png")
