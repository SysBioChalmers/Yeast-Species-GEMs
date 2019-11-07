# -*- coding: utf-8 -*-
# -*- python 3 -*-
# -*- hongzhong Lu -*-

import os
import csv
os.chdir('/Users/luho/Documents/GitHub/Pan_Model_332_Yeast_Species/draft_GEM_all_yeast_species')



strain = os.listdir('/Users/luho/Documents/GitHub/Pan_Model_332_Yeast_Species/draft_GEM_all_yeast_species/strain specific model from RAVEN_biocyc_55_110/')
strain = [x for x in strain if x != '.DS_Store']
for x in strain:
    print(x)
    data = open('/Users/luho/Documents/GitHub/Pan_Model_332_Yeast_Species/draft_GEM_all_yeast_species/strain specific model from RAVEN_biocyc_55_110/' + x +'/excelRxns.txt', 'r').readlines()
    with open('/Users/luho/Documents/GitHub/Pan_Model_332_Yeast_Species/draft_GEM_all_yeast_species/strain specific model from RAVEN_biocyc_55_110/' + x +'/draft_GEM.tsv', mode='w') as file:
        employee_writer = csv.writer(file, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for i, line in enumerate(data):
            line = line.replace("#", "").lstrip()
            line = line.strip('\n')
            line = line.strip('\t')
            #if i != 0:
            line0 = line.split("\t")
            print(line0)
            employee_writer.writerow(line0)


'''
# remove file using python code
my_dir = "/Users/luho/PycharmProjects/3D_model/refine draft strain specific GEM for 332 yeast species/strain specific model from RAVEN_biocyc_55_110/"
for fname in os.listdir(my_dir):
    print(fname)
    if fname !='.DS_Store':
        newfname = os.path.join(my_dir, fname)
    for name2 in os.listdir(newfname):
        print(name2)
        if name2.startswith("csv"):
            print(name2)
            os.remove(os.path.join(newfname, name2))
'''

strain = os.listdir('/Users/luho/Documents/GitHub/Pan_Model_332_Yeast_Species/draft_GEM_all_yeast_species/strain specific model from RAVEN_kegg/')
strain = [x for x in strain if x != '.DS_Store']
for x in strain:
    print(x)
    data = open('/Users/luho/Documents/GitHub/Pan_Model_332_Yeast_Species/draft_GEM_all_yeast_species/strain specific model from RAVEN_kegg/' + x +'/excelRxns.txt', 'r').readlines()
    with open('/Users/luho/Documents/GitHub/Pan_Model_332_Yeast_Species/draft_GEM_all_yeast_species/strain specific model from RAVEN_kegg/' + x +'/draft_GEM.tsv', mode='w') as file:
        employee_writer = csv.writer(file, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for i, line in enumerate(data):
            line = line.replace("#", "").lstrip()
            line = line.strip('\n')
            line = line.strip('\t')
            #if i != 0:
            line0 = line.split("\t")
            print(line0)
            line0 = line0[0:11]
            employee_writer.writerow(line0)