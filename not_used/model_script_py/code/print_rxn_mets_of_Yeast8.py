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
sys.path.append(r"/Users/luho/PycharmProjects/3D_model/model_script_py/code")
os.chdir('/Users/luho/PycharmProjects/3D_model/model_script_py/code')
from mainFunction import *


# load the GEM model from the yeast GEM repo
yeastGEM = read_sbml_model("/Users/luho/Documents/GitHub/YeastMetabolicNetwork-GEM/ModelFiles/xml/yeastGEM.xml")

yeastGEM = correctSomeWrongFormat(yeastGEM)
genes = []
for gene in yeastGEM.genes:
    x = gene.id
    genes.append(x)

# obtain the all the metabolite information
met_yeast8 = produceMetaboliteList(yeastGEM)
met_yeast8['chebi'] = met_yeast8['chebi'].str.replace('CHEBI:','chebi:')
met_yeast8['kegg'] = 'kegg:' + met_yeast8['kegg']
# find the mnxid based on keggid and chebiid
# this will a function
mnx_ref = pd.read_csv('../data/chem_xref.tsv', sep="\t")
mnx_ref_kegg = mnx_ref[mnx_ref['XREF'].str.contains('kegg')]
mnx_ref_chebi = mnx_ref[mnx_ref['XREF'].str.contains('chebi')]

met_yeast8['MNXID_kegg'] = multiMapping(mnx_ref_kegg['MNX_ID'],mnx_ref_kegg['XREF'],met_yeast8['kegg'])
met_yeast8['MNXID_new'] = multiMapping(mnx_ref_chebi['MNX_ID'],mnx_ref_chebi['XREF'],met_yeast8['chebi'])

# obtain all the reaction information
rxn_yeast8 = produceRxnList(yeastGEM)

# save the result
saveExcel(met_yeast8, '../result/met_yeast8.xlsx')
saveExcel(rxn_yeast8, '../result/rxn_yeast8.xlsx')