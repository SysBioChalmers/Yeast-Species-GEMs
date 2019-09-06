# -*- coding: utf-8 -*-
# -*- python 3 -*-
# -*- hongzhong Lu -*-

'''
kegg web service annotation of pan-genome
estimate the homolog gene number from this method
'''

import os    ##for directory
import pandas as pd
import sys

# set the directory
sys.path.append(r"/Users/luho/PycharmProjects/3D_model/find_homolog_for_panID_py/code")
os.chdir('/Users/luho/PycharmProjects/3D_model/find_homolog_for_panID_py/code')
from mainFunction import *

f1_kegg = pd.read_excel('../data/fasta1_kegg.xlsx', header=None)
f2_kegg = pd.read_excel('../data/fasta2_kegg.xlsx', header=None)
f3_kegg = pd.read_excel('../data/fasta3_kegg.xlsx', header=None)

# combine the dataframe
# combine the interface and remove the duplicated one
annotation_kegg = pd.concat([f1_kegg, f2_kegg, f3_kegg], axis=0, join='outer')
annotation_kegg.columns = ['panID', 'KO']
# remove the panID with no KO
annotation_kegg0 = annotation_kegg.dropna()
annotation_kegg0 = annotation_kegg0.reset_index()
annotation_kegg0 = annotation_kegg0[['panID', 'KO']]

'''
# analyse the KO
KO_summary = annotation_kegg0['KO'].value_counts()
KO_summary0 = pd.DataFrame({'KO':KO_summary.index, 'num':KO_summary.values})
# visualization
import seaborn as sns
import matplotlib.pyplot as plt
fig = plt.figure(figsize=(8,6))
ax =sns.distplot(KO_summary0['num'])
ax.set(xlabel='Number', ylabel='Frequency')
ax.set(xlim=(0, 50))
plt.show()
'''

# Group the gene into reference from s288c and other species
ref_gene = annotation_kegg0[annotation_kegg0['panID'].str.contains('Saccharomyces_cerevisiae')]
ref_gene['geneID'] = ref_gene['panID'].str.replace('Saccharomyces_cerevisiae@','')
# next, input all the gene annotation of s288c from kegg database
sce_kegg = pd.read_csv('../data/sce_ko_2019_8.txt', sep="\t", header=None)
sce_kegg.columns = ['geneID','KO']
sce_kegg = sce_kegg.dropna()
# combine the above annotation
gene_288c = list(set(sce_kegg['geneID'])-set(ref_gene['geneID']))
ref_gene2 = sce_kegg[sce_kegg['geneID'].isin(gene_288c)]
ref_gene2['panID'] = [None]*len(ref_gene2['geneID'])
sce_all = pd.concat([ref_gene, ref_gene2], axis=0, join='outer')

# get the non-reference gene not from s288c
nonRefORF = annotation_kegg0[~annotation_kegg0['panID'].isin(sce_all['panID'])]


# then we have tasks
# the task is to find the genes belong to the same KO for genes from references strain s288c and other species
'''
Note: here, it seems that multiple panID could have the same KO of one gene from s288c.
'''
sameKO =  list(set(nonRefORF['KO']) & set(sce_all['KO']))
sameKO_sce_kegg = sce_all[sce_all['KO'].isin(sameKO)]
sameKO_sce_kegg['panID'] = multiMapping(nonRefORF['panID'], nonRefORF['KO'], sameKO_sce_kegg['KO'])
# calculate the number of non reference geneID for each gene from s288c
# it should be noted that different genes from s288c could have the same KO id
list0 = sameKO_sce_kegg['panID'].tolist()
num_count = [x.count(';') +1  for x in list0]
sameKO_sce_kegg['panID_count'] = num_count
saveExcel(sameKO_sce_kegg, '../result/sameKO_sce_kegg.xlsx')
