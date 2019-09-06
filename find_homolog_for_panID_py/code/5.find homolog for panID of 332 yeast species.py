# -*- coding: utf-8 -*-
# -*- python 3 -*-
# -*- hongzhong Lu -*-

'''
Using the genome annotation from eggnog,kegg web service and bi-directional blast analysis,
Firstly we need find the best hit for new pan proteins from s288c which connect the reaction in yeastGEM
Then these new pan proteins could be merged with the yeastGEM to formulate the universal panYeast as the first steps.
Next we can add the new reactions catalyzed by all the new pan proteins, with which we get the final universal panYeast model
'''
import os    ##for directory
import pandas as pd
import sys
from cobra.io import read_sbml_model

# set the directory
sys.path.append(r"/Users/luho/PycharmProjects/3D_model/find_homolog_for_panID_py/code")
os.chdir('/Users/luho/PycharmProjects/3D_model/find_homolog_for_panID_py/code')
from mainFunction import *


# load the GEM model from the yeast GEM repo
yeastGEM = read_sbml_model("/Users/luho/Documents/GitHub/YeastMetabolicNetwork-GEM/ModelFiles/xml/yeastGEM.xml")

yeastGEM = correctSomeWrongFormat(yeastGEM)
genes = []
for gene in yeastGEM.genes:
    x = gene.id
    genes.append(x)

# input the blast analysis result using diamond tool
blast_pan_to_sce = pd.read_csv('../data/blastNonRefORFwithSGD_332_yeasts.txt', sep="\t", header= None)
blast_sce_to_pan = pd.read_csv('../data/blastSGDbasedNonRefORF_332_yeasts.txt', sep="\t", header= None)
# if the pidentity is not used:
pident0 = 70
result_bbh = bbhFilterWithPident(a_to_b=blast_pan_to_sce, b_to_a= blast_sce_to_pan, pident=pident0)
homolog_pan = result_bbh[result_bbh['hitID'].isin(genes)]



'''
Note: as a panID could have multiple hit from s288c in BBH analysis, should we further refine it and 
choose a best one based on pidenty and bitscore?
'''
len(homolog_pan['geneID'].unique())
len(homolog_pan['hitID'].unique())

panID0 = homolog_pan['geneID'].unique().tolist()
hitID = findBestHitFromBlast(panID=panID0, blast_inf=homolog_pan)
pan_hit_mapping = pd.DataFrame({'panID':panID0,'hitID':hitID})
outfile0 = "../result/pan_hit_mapping for panYeast_v2_PI@" + str(pident0) + ".xlsx"
saveExcel(pan_hit_mapping, outfile0)


