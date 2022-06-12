# -*- coding: utf-8 -*-
# -*- python 3 -*-
# -*- hongzhong Lu -*-

import os
dir0 = os.getcwd()
os.chdir(dir0)


my_dir = os.path.join(dir0, "strain_specific_model_from_RAVEN_biocyc_55_110/")
strain = os.listdir(my_dir)

strain = [x for x in strain if x != '.DS_Store']
for fname in strain:
    print(fname)
    if fname !='.DS_Store':
        newfname = os.path.join(my_dir, fname)
    for name2 in os.listdir(newfname):
        print(name2)
        if name2.startswith("excelComps"):
            print(name2)
            os.remove(os.path.join(newfname, name2))

for fname in strain:
    print(fname)
    if fname !='.DS_Store':
        newfname = os.path.join(my_dir, fname)
    for name2 in os.listdir(newfname):
        print(name2)
        if name2.startswith("excelModel"):
            print(name2)
            os.remove(os.path.join(newfname, name2))





my_dir = os.path.join(dir0, "strain_specific_model_from_RAVEN_kegg/")
strain = os.listdir(my_dir)

strain = [x for x in strain if x != '.DS_Store']
for fname in strain:
    print(fname)
    if fname !='.DS_Store':
        newfname = os.path.join(my_dir, fname)
    for name2 in os.listdir(newfname):
        print(name2)
        if name2.startswith("excelComps"):
            print(name2)
            os.remove(os.path.join(newfname, name2))

for fname in strain:
    print(fname)
    if fname !='.DS_Store':
        newfname = os.path.join(my_dir, fname)
    for name2 in os.listdir(newfname):
        print(name2)
        if name2.startswith("excelModel"):
            print(name2)
            os.remove(os.path.join(newfname, name2))
