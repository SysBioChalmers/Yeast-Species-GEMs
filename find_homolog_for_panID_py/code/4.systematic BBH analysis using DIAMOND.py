# -*- coding: utf-8 -*-
# -*- python 3 -*-
# -*- hongzhong Lu -*-

'''
How to find the homology genes from pan-genome for genes from s288c
Maybe the bi-directional analysis is a simple and best method!!!
'''

import os    ##for directory
import pandas as pd
import sys
from cobra.io import read_sbml_model


# set the directory
sys.path.append(r"/Users/luho/PycharmProjects/3D_model/find_homolog_for_panID_py/code")
os.chdir('/Users/luho/PycharmProjects/3D_model/find_homolog_for_panID_py/code')
from mainFunction import *

# input the blast analysis result using diamond tool
blast_pan_to_sce = pd.read_csv('../data/blastNonRefORFwithSGD_332_yeasts.txt', sep="\t", header= None)
blast_sce_to_pan = pd.read_csv('../data/blastSGDbasedNonRefORF_332_yeasts.txt', sep="\t", header= None)
# if the pidentity is not used:
result_bbh = bbhFilterWithPident(a_to_b=blast_pan_to_sce, b_to_a= blast_sce_to_pan, pident=60)



# how the pidentity affect the number of homology gene number from BBH analysis
pident0 = []
all_homolog = []
for p in range(25,100,5):
    print(p)
    result_bbh = bbhFilterWithPident(blast_pan_to_sce, blast_sce_to_pan, pident=p)
    homolog_num = len(result_bbh['geneID'].unique())
    pident0.append(p)
    all_homolog.append(homolog_num)

# plot
simpleLineXY(x=pident0, y=all_homolog, x_title='pidentity(%)', y_title='homolog gene number', title='All gene')




# considering the metabolic gene from yeast GEM repo
# load the GEM model from the yeast GEM repo
yeastGEM = read_sbml_model("/Users/luho/Documents/GitHub/YeastMetabolicNetwork-GEM/ModelFiles/xml/yeastGEM.xml")
genes = []
for gene in yeastGEM.genes:
    x = gene.id
    if '__45__' in x:
        x = x.replace('__45__','-')
    print(x)
    genes.append(x)

# filter t he metabolic gene
blast_pan_to_sce0 = blast_pan_to_sce[blast_pan_to_sce[1].isin(genes)]
blast_sce_to_pan0 = blast_sce_to_pan[blast_sce_to_pan[0].isin(genes)]

# how the pidentity affect the number of homology gene number from BBH analysis
pident0 = []
all_homolog = []
for p in range(25,100,5):
    print(p)
    result_bbh = bbhFilterWithPident(blast_pan_to_sce0, blast_sce_to_pan0, pident=p)
    homolog_num = len(result_bbh['geneID'].unique())
    pident0.append(p)
    all_homolog.append(homolog_num)

# plot
simpleLineXY(x=pident0, y=all_homolog, x_title='pidentity(%)', y_title='homolog gene number', title='Metabolic')


# if we only use the single blast analysis from pan to sce 288c
pident0 = []
all_homolog = []
for p in range(25,100,5):
    print(p)
    result_bbh = singleBlastFilterWithPident(blast_pan_to_sce0, pident=p)
    homolog_num = len(result_bbh['geneID'].unique())
    pident0.append(p)
    all_homolog.append(homolog_num)

# plot
simpleLineXY(x=pident0, y=all_homolog, x_title='pidentity(%)', y_title='homolog gene number', title='Metabolic_simple_blast')
