# -*- coding: utf-8 -*-
# -*- python 3 -*-
# -*- hongzhong Lu -*-

'''
Build pan-model using RAVEN2
estimate the homolog gene number from this method
'''
import os    ##for directory
import pandas as pd
import sys

# set the directory
sys.path.append(r"/Users/luho/PycharmProjects/3D_model/find_homolog_for_panID_py/code")
os.chdir('/Users/luho/PycharmProjects/3D_model/find_homolog_for_panID_py/code')
from mainFunction import *

panYeast2_biocyc = pd.read_excel('../data/panYeast2_biocyc/excelRxns.xlsx')
panYeast2_biocyc0 = panYeast2_biocyc[['ID','GENE ASSOCIATION']]
panYeast2_biocyc1 = splitAndCombine(panYeast2_biocyc0['GENE ASSOCIATION'], panYeast2_biocyc0['ID'], sep0=' or ', moveDuplicate=True)


panYeast2_kegg = pd.read_excel('../data/panYeast2_kegg/excelRxns.xlsx')
panYeast2_kegg0 = panYeast2_kegg[['ID','GENE ASSOCIATION']]
panYeast2_kegg1 = splitAndCombine(panYeast2_kegg0['GENE ASSOCIATION'], panYeast2_kegg0['ID'], sep0=' or ', moveDuplicate=True)


sce_biocyc = pd.read_excel('../data/s288c_biocyc/excelRxns.xlsx')
sce_biocyc0 = sce_biocyc[['ID','GENE ASSOCIATION']]
sce_biocyc1 = splitAndCombine(sce_biocyc0['GENE ASSOCIATION'], sce_biocyc0['ID'], sep0=' or ', moveDuplicate=True)

sce_kegg = pd.read_excel('../data/s288c_kegg/excelRxns.xlsx')
sce_kegg0 = sce_kegg[['ID','GENE ASSOCIATION']]
sce_kegg1 = splitAndCombine(sce_kegg0['GENE ASSOCIATION'], sce_kegg0['ID'], sep0=' or ', moveDuplicate=True)


# In the pan model we remove the gene from s288c and then we get the gene for each rxn only from
# non reference genome
panYeast2_biocyc1 = panYeast2_biocyc1[~panYeast2_biocyc1['V2'].str.contains('Saccharomyces_cerevisiae')]
panYeast2_biocyc1 = panYeast2_biocyc1[~panYeast2_biocyc1['V2'].str.contains('NA')]
panYeast2_biocyc1 = panYeast2_biocyc1.dropna()

panYeast2_kegg1 = panYeast2_kegg1[~panYeast2_kegg1['V2'].str.contains('Saccharomyces_cerevisiae')]
panYeast2_kegg1 = panYeast2_kegg1[~panYeast2_kegg1['V2'].str.contains('NA')]
panYeast2_kegg1 = panYeast2_kegg1.dropna()


# establish the mapping for the common reactions catalyzed by proteins from non refrence genome and s288c
common_rxn_biocyc = list(set(panYeast2_biocyc1['V1']) & set(sce_biocyc1['V1']))
common_annot_biocyc = panYeast2_biocyc1[panYeast2_biocyc1['V1'].isin(common_rxn_biocyc)]

common_rxn_kegg = list(set(panYeast2_kegg1['V1']) & set(sce_kegg1['V1']))
common_annot_kegg = panYeast2_kegg1[panYeast2_kegg1['V1'].isin(common_rxn_kegg)]

# calculate the panID from biocyc and kegg together
panID_kegg = list(common_annot_kegg['V2'].unique())
panID_biocyc =list(common_annot_biocyc['V2'].unique())
panID_all = list(set(common_annot_biocyc['V2']) | set(common_annot_kegg['V2']))
panID_intersect = list(set(common_annot_biocyc['V2']) & set(common_annot_kegg['V2']))

# plot simple bar chart
label = ['kegg', 'biocyc', 'kegg_or_biocyc', 'kegg_and_biocyc']
number = [len(panID_kegg),len(panID_biocyc),len(panID_all),len(panID_intersect)]
simpleBarPlot(x=label, y=number, y_title='Homolog gene number',title='annotation_RAVEN')
