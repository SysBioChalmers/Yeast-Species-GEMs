# -*- coding: utf-8 -*-
# -*- python 3 -*-
# -*- hongzhong Lu -*-

'''
The pan genome annotation could be divided into the followed three steps
1. EggNog annotation
2. KEGG annotation
3. RAVEN2 based on bioCyc database
4. TCDB transport reaction annotation
5. Blast with present panYeast??? Here we can only focus the 400 reactions not in sce model
   As we will built new version of panYeast based on Yeast8
Once the annotation is obtained, a find_homolog_for_panID_py will be obtained.
'''
import os    ##for directory
import pandas as pd
import sys

# set the directory
sys.path.append(r"/Users/luho/PycharmProjects/3D_model/find_homolog_for_panID_py/code")
os.chdir('/Users/luho/PycharmProjects/3D_model/find_homolog_for_panID_py/code')
from mainFunction import *

# interesting columns from EggNOG annotation
headers = ['query_name', 'KEGG_TC'] # transport protein
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


# Group the gene into reference from s288c and other species
ref_gene = annotation_eggnog0[annotation_eggnog0['panID'].str.contains('Saccharomyces_cerevisiae')]
ref_gene['geneID'] = ref_gene['panID'].str.replace('Saccharomyces_cerevisiae@','')


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
# combine the above annotation
gene_288c = list(set(sce_eggnog['geneID'])-set(ref_gene['geneID']))
ref_gene2 = sce_eggnog[sce_eggnog['geneID'].isin(gene_288c)]
ref_gene2['panID'] = [None]*len(ref_gene2['geneID'])
sce_all = pd.concat([ref_gene, ref_gene2], axis=0, join='outer')
# get the non-reference gene not from s288c
nonRefORF = annotation_eggnog0[~annotation_eggnog0['panID'].isin(sce_all['panID'])]


# then we have tasks
# the first task is to find the new KO from the non reference genes
# the second task is to find the genes belong to the same KO for genes from references strain s288c and other species

# task1
sameKO =  list(set(nonRefORF['KO']) & set(sce_all['KO']))
sameKO_sce = sce_all[sce_all['KO'].isin(sameKO)]
sameKO_sce['panID'] = multiMapping(nonRefORF['panID'], nonRefORF['KO'], sameKO_sce['KO'])
# calculate the number of non reference geneID for each gene from s288c
# it should be noted that different genes from s288c could have the same KO id
list0 = sameKO_sce['panID'].tolist()
num_count = [x.count(';') +1  for x in list0]
sameKO_sce['panID_count'] = num_count
saveExcel(sameKO_sce, '../result/sameKO_sce_eggnog_transporter_protein.xlsx')