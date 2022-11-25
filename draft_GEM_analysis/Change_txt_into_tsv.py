# -*- coding: utf-8 -*-
# -*- python 3 -*-
# -*- hongzhong Lu -*-

import os
import csv
os.chdir('/Users/luho/PycharmProjects/3D_model/draft_GEM_all_yeast_species')



strain = os.listdir('../ComplementaryData/draft_yeast_GEMs/strain_specific_model_from_RAVEN_biocyc_55_110/')
strain = [x for x in strain if x != '.DS_Store']
for x in strain:
    print(x)
    data = open('../ComplementaryData/draft_yeast_GEMs/strain_specific_model_from_RAVEN_biocyc_55_110/' + x +'/excelRxns.txt', 'r').readlines()
    with open('../ComplementaryData/draft_yeast_GEMs/strain_specific_model_from_RAVEN_biocyc_55_110/' + x +'/draft_GEM.tsv', mode='w') as file:
        employee_writer = csv.writer(file, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for i, line in enumerate(data):
            line = line.replace("#", "").lstrip()
            line = line.strip('\n')
            line = line.strip('\t')
            #if i != 0:
            line0 = line.split("\t")
            print(line0)
            employee_writer.writerow(line0)



strain = os.listdir('../ComplementaryData/draft_yeast_GEMs/strain_specific_model_from_RAVEN_kegg/')
strain = [x for x in strain if x != '.DS_Store']
for x in strain:
    print(x)
    data = open('../ComplementaryData/draft_yeast_GEMs/strain_specific_model_from_RAVEN_kegg/' + x +'/excelRxns.txt', 'r').readlines()
    with open('../ComplementaryData/draft_yeast_GEMs/strain_specific_model_from_RAVEN_kegg/' + x +'/draft_GEM.tsv', mode='w') as file:
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