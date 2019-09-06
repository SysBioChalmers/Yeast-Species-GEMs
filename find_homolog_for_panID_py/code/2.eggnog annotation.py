# -*- coding: utf-8 -*-
# -*- python 3 -*-
# -*- hongzhong Lu -*-

'''
EggNog annotation of pan-genome
estimate the homolog gene number from this method
'''
import os    ##for directory
import pandas as pd
import sys
from cobra.io import read_sbml_model

# set the directory
sys.path.append(r"/Users/luho/PycharmProjects/3D_model/find_homolog_for_panID_py/code")
os.chdir('/Users/luho/PycharmProjects/3D_model/find_homolog_for_panID_py/code')
from mainFunction import *

# interesting columns from EggNOG annotation
headers = ['query_name', 'KEGG_ko',] # ko annotation
new_name = ['panID', 'KO']

dir_input1 = '../data/query_seqs.fa.emapper.annotations1'
f1_eggnog0 = eggnogInput(dir_input1)
f1_eggnog = f1_eggnog0[headers]
f1_eggnog.columns = new_name


dir_input2 = '../data/query_seqs.fa.emapper.annotations2'
f2_eggnog0 = eggnogInput(dir_input2)
f2_eggnog = f2_eggnog0[headers]
f2_eggnog.columns = new_name

dir_input3 = '../data/query_seqs.fa.emapper.annotations3'
f3_eggnog0 = eggnogInput(dir_input3)
f3_eggnog = f3_eggnog0[headers]
f3_eggnog.columns = new_name

# combine the dataframe
# combine the interface and remove the duplicated one
annotation_eggnog = pd.concat([f1_eggnog, f2_eggnog, f3_eggnog], axis=0, join='outer')
annotation_eggnog.columns = ['panID', 'KO']
# remove the panID with no KO
annotation_eggnog0 = annotation_eggnog.dropna()
annotation_eggnog0 = annotation_eggnog0[annotation_eggnog0['KO'] != '']
annotation_eggnog0 = annotation_eggnog0.reset_index()
annotation_eggnog0 = annotation_eggnog0[['panID', 'KO']]
annotation_eggnog0['KO'] = annotation_eggnog0['KO'] .str.replace('ko:','')

# as on panID has multiple KO, we need split it
annotation_eggnog_correct = splitAndCombine(annotation_eggnog0['KO'], annotation_eggnog0['panID'], sep0=',', moveDuplicate=True)
annotation_eggnog_correct.columns = ['panID','KO']
annotation_eggnog0 = annotation_eggnog_correct

# display the result
plotDensityProfile(ss=annotation_eggnog0['KO'], x_title='Number', y_title='Frequency')


# next, input all the gene annotation of s288c from eggnog database
dir_input4 = '../data/sce288c.fa.emapper.annotations'
f4_eggnog0 = eggnogInput(dir_input4)
f4_eggnog = f4_eggnog0[headers]
f4_eggnog.columns = new_name
f4_eggnog = f4_eggnog.dropna()
sce_eggnog = f4_eggnog[f4_eggnog['KO'] != '']
# as on panID has multiple KO, we need split it
sce_eggnog_correct = splitAndCombine(sce_eggnog['KO'], sce_eggnog['panID'], sep0=',', moveDuplicate=True)
sce_eggnog_correct.columns = ['panID','KO']
sce_eggnog = sce_eggnog_correct
sce_eggnog.columns = ['geneID','KO']
sce_eggnog = sce_eggnog.dropna()
sce_eggnog['KO'] = sce_eggnog['KO'].str.replace('ko:','')

# get the non-reference gene not from s288c
nonRefORF = annotation_eggnog0[~annotation_eggnog0['panID'].str.contains('Saccharomyces_cerevisiae')]


# then we have tasks
# the first task is to find the new KO from the non reference genes
# the second task is to find the genes belong to the same KO for genes from references strain s288c and other species

# task1
sameKO =  list(set(nonRefORF['KO']) & set(sce_eggnog['KO']))
sameKO_sce = sce_eggnog[sce_eggnog['KO'].isin(sameKO)]
sameKO_sce['panID'] = multiMapping(nonRefORF['panID'], nonRefORF['KO'], sameKO_sce['KO'])
# calculate the number of non reference geneID for each gene from s288c
# it should be noted that different genes from s288c could have the same KO id
list0 = sameKO_sce['panID'].tolist()
num_count = [x.count(';') +1  for x in list0]
sameKO_sce['panID_count'] = num_count

# check how many genes in yeast GEM
# load the GEM model from the yeast GEM repo
yeastGEM = read_sbml_model("/Users/luho/Documents/GitHub/YeastMetabolicNetwork-GEM/ModelFiles/xml/yeastGEM.xml")
yeastGEM = correctSomeWrongFormat(yeastGEM)
genes = []
for gene in yeastGEM.genes:
    x = gene.id
    genes.append(x)
sameKO_sce_GEM = sameKO_sce[sameKO_sce['geneID'].isin(genes)]
# An exploration: first establish the mapping between reference gene and panID , then get the blast result
test = splitAndCombine(sameKO_sce_GEM['panID'], sameKO_sce_GEM['geneID'], sep0=';', moveDuplicate=True)
test.columns = ['geneID','panID']
saveExcel(test, '../result/sameKO_sce_eggnog.xlsx')

