# -*- coding: utf-8 -*-
# -*- python 3 -*-
# -*- hongzhong Lu -*-

'''
The code is  used to update the gene existence matrix
'''

import os    ##for directory
import pandas as pd
import sys

# set the directory
sys.path.append(r"/Users/luho/PycharmProjects/3D_model/find_homolog_for_panID_py/code")
os.chdir('/Users/luho/PycharmProjects/3D_model/find_homolog_for_panID_py/code')
from mainFunction import *
OG_pan = pd.read_csv('../data/representatives.tsv', sep="\t")

OG_pan_sce = OG_pan[OG_pan["representative"].str.contains("Saccharomyces_cerevisiae@")]
OG_pan_sce["representative"] = OG_pan_sce["representative"].str.replace("Saccharomyces_cerevisiae@", "")

OG_pan_sce.to_csv("../result/OG_pan_sce.csv")

# input the gene list of sce which need to be updated in the gene existence matrix
sce_list = pd.read_excel("../data/Genelist_genestrainMatrix.xlsx")
sce_list["OG"] = singleMapping (OG_pan_sce["ortho_id"], OG_pan_sce["representative"], sce_list["newID"], dataframe=True)

# read the gene existence matrix
gene_table = open("../data/gene_pa_table.csv").readlines()
row_name = gene_table[0]
row_name0 = row_name.split(",")
# fist build dict for the replace
# replace "Q0050" into "Q0055" or "OG7640"
dict1 = {}
dict2 = {}
for index, rows in sce_list.iterrows():
    print(rows["oldID"])
    dict1[rows["oldID"]] = rows["newID"]
    dict2[rows["oldID"]] = rows["OG"]
# method1
row_name1 = []
for ele in row_name0:
    if ele in dict1.keys():
        new_ele = dict1[ele]
        row_name1.append(new_ele)
    else:
        row_name1.append(ele)
row_name_s1 = ",".join(row_name1)
# check
"Q0055" in row_name_s1
"Q0050" in row_name_s1



# method2
row_name2 = []
for ele in row_name0:
    if ele in dict2.keys():
        new_ele = dict2[ele]
        row_name2.append(new_ele)
    else:
        row_name2.append(ele)
row_name_s2 = ",".join(row_name2)

# check
"Q0055" in row_name_s2
"Q0050" in row_name_s2
"OG7640" in row_name_s2

# replace the original sce id with the OG id
gene_table[0] = row_name_s2
out= "../result/gene_pa_table_update.csv"

def write_list_to_file(guest_list, filename):
    """Write the list to csv file."""

    with open(filename, "w") as outfile:
        for entries in guest_list:
            outfile.write(entries)
            outfile.write("\n")

write_list_to_file(guest_list=gene_table, filename=out)