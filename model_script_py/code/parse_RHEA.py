# -*- coding: utf-8 -*-
# -*- python 3 -*-
# -*- hongzhong Lu -*-

# parse the rhea database

import os    ##for directory
import sys
import pandas as pd

# set the directory
sys.path.append(r"/Users/luho/PycharmProjects/3D_model/model_script_py/code")
os.chdir('/Users/luho/PycharmProjects/3D_model/model_script_py/code')
from mainFunction import *


infile = '../data/rhea-kegg.reaction'
line= []
rxn_rhea = open(infile, 'r').readlines()
for i, x in enumerate(rxn_rhea):
    if '///' in x:
        line.append(i)
rheaID = []
formula = []
equation = []
enzyme = []
for i, x in enumerate(line):
    print(i)
    if x ==8:
        s0 = rxn_rhea[0:x]
        r1 = [x for x in s0 if 'ENTRY' in x][0]
        r10 = r1.replace('ENTRY', '').strip()
        r2 = [x for x in s0 if 'DEFINITION' in x][0]
        r20 = r2.replace('DEFINITION', '').strip()
        r3 = [x for x in s0 if 'EQUATION' in x][0]
        r30 = r3.replace('EQUATION', '').strip()
        r4 = [x for x in s0 if 'ENZYME' in x][0]
        r40 = r4.replace('ENZYME', '').strip()
    else:
        j = line[i-1]
        s0 = rxn_rhea[j:x]
        r1 = [x for x in s0 if 'ENTRY' in x][0]
        r10 = r1.replace('ENTRY', '').strip()
        r2 = [x for x in s0 if 'DEFINITION' in x][0]
        r20 = r2.replace('DEFINITION', '').strip()
        r3 = [x for x in s0 if 'EQUATION' in x][0]
        r30 = r3.replace('EQUATION', '').strip()
        r4 = [x for x in s0 if 'ENZYME' in x]
        if len(r4):
            r40 = r4[0]
            r40 = r40.replace('ENZYME', '').strip()
        else:
            r40 = None
    print(r10)
    rheaID.append(r10)
    formula.append(r20)
    equation.append(r30)
    enzyme.append(r40)

rhea_summary = pd.DataFrame({'rheaID': rheaID, 'formula': formula, 'equation': equation, 'EC': enzyme})
saveExcel(rhea_summary, "../result/rhea_reaction_summary.xlsx")







