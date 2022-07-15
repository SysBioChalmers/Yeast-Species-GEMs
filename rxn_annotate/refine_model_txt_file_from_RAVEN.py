# -*- coding: utf-8 -*-
# -*- python 3 -*-
# -*- hongzhong Lu -*-

import os
import csv
os.chdir('/Users/luho/PycharmProjects/3D_model/Reaction_and_metabolite_standardization')



strain = os.listdir('/Users/luho/PycharmProjects/3D_model/Reaction_and_metabolite_standardization/data/strain specific model from RAVEN_biocyc_55_110/')
strain = [x for x in strain if x != '.DS_Store']
for x in strain:
    print(x)
    data = open('/Users/luho/PycharmProjects/3D_model/Reaction_and_metabolite_standardization/data/strain specific model from RAVEN_biocyc_55_110/' + x +'/excelRxns.txt', 'r').readlines()
    with open('/Users/luho/PycharmProjects/3D_model/Reaction_and_metabolite_standardization/data/strain specific model from RAVEN_biocyc_55_110/' + x +'/csvRxns.csv', mode='w') as file:
        employee_writer = csv.writer(file, delimiter='@', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for i, line in enumerate(data):
            line = line.replace("#", "").lstrip()
            if i != 0:
                line0 = line.split("\t")
                print(line0)
                employee_writer.writerow(line0)




strain = os.listdir('/Users/luho/PycharmProjects/3D_model/Reaction_and_metabolite_standardization/data/strain specific model from RAVEN_kegg/')
strain = [x for x in strain if x != '.DS_Store']
for x in strain:
    print(x)
    data = open('/Users/luho/PycharmProjects/3D_model/Reaction_and_metabolite_standardization/data/strain specific model from RAVEN_kegg/' + x +'/excelRxns.txt', 'r').readlines()
    with open('/Users/luho/PycharmProjects/3D_model/Reaction_and_metabolite_standardization/data/strain specific model from RAVEN_kegg/' + x +'/csvRxns.csv', mode='w') as file:
        employee_writer = csv.writer(file, delimiter='@', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for i, line in enumerate(data):
            line = line.replace("#", "").lstrip()
            if i != 0:
                line0 = line.split("\t")
                print(line0)
                employee_writer.writerow(line0)