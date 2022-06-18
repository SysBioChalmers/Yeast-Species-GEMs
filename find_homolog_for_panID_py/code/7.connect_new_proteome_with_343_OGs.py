# -*- coding: utf-8 -*-
# -*- python 3 -*-
# -*- hongzhong Lu -*-

'''
The code is  used to connect the new proteome with the OGs groups of 332 yeast species plus 11 fungal outgroup
'''

import os    ##for directory
import pandas as pd
import sys
from Bio import SeqIO
import itertools

# set the directory
sys.path.append(r"/Users/luho/PycharmProjects/3D_model/find_homolog_for_panID_py/code")
os.chdir('/Users/luho/PycharmProjects/3D_model/find_homolog_for_panID_py/code')
from mainFunction import *

#check which file store the blast result
#this is directory to store the blast data and result
blast_data = "/Users/luho/Documents/blast/data/"
blast_result = "/Users/luho/Documents/blast/result/"

# firstly calculate the lenght of each protein
# input the new proteome and 343 proteome


# species1
infile0 = blast_data + 'Diutina_rugosa.faa'
infile = blast_data + '343taxa_proteins.fasta'
new_blast_all_input =blast_result + "blastnewbasedAll_yeasts.txt"
all_blast_new_input = blast_result + "blastAllORFwithnew.txt"
output_file = "../result/" + 'Diutina_rugosa' + '.csv'

# species2
infile0 = blast_data + 'Scheffersomyces_stipitis.fasta'
infile = blast_data + '343taxa_proteins.fasta'
new_blast_all_input =blast_result + "blastnewbasedAll_yeasts2.txt"
all_blast_new_input = blast_result + "blastAllORFwithnew2.txt"
output_file = "../result/" + 'Scheffersomyces_stipitis' + '.csv'



# the main procedure to do the id mapping
records = list(SeqIO.parse(infile0, "fasta"))
records_343 = list(SeqIO.parse(infile, "fasta"))

gene_id_length_new_species = {}
for record in records:
    print(record.id)
    gene_id_length_new_species[record.id] = len(record.seq)

gene_id_length_343 = {}
for record in records_343:
    print(record.id)
    gene_id_length_343[record.id] = len(record.seq)


## input and preprocess the blast result file
new_blast_all = pd.read_csv(new_blast_all_input, sep="\t", header= None)
ss0 = ["geneID", "hitID", "pident", "length", "mismatch", "gapopen",
       "qstart", "qend", "sstart", "send", "evalues", "bitscore"]
new_blast_all.columns = ss0
new_blast_all["combine"] =  new_blast_all["geneID"] + "&&" + new_blast_all["hitID"]
new_blast_all_filter1 = new_blast_all[new_blast_all["evalues"] < 0.00001]
print(len(new_blast_all["geneID"].unique())) # total 5212 genes could find the best hit
print(len(new_blast_all_filter1["geneID"].unique())) # total 5212 genes could find the best hit

## input the reverse blast result
blast_all = pd.read_csv(all_blast_new_input, sep="\t", header= None)
blast_all.columns = ss0
blast_all["combine"] = blast_all["hitID"] + "&&" + blast_all["geneID"]
blast_all_filter1 = blast_all[blast_all["evalues"] < 0.00001]

## get the common pairs
pair1=new_blast_all_filter1["combine"].tolist()
pair2=blast_all_filter1["combine"].tolist()
#pair3 =[x for x in pair1 if x in pair2] # this is low efficent
pair3 = list(set(pair1) & set(pair2))


new_blast_all_f2 = new_blast_all_filter1[new_blast_all_filter1["combine"].isin(pair3)]
print(len(new_blast_all_f2["geneID"].unique())) # total 5205 genes

## match lenght filteration
match_percent = 0.4 # it seems that 0.4 is an optional value
new_blast_all_f2["matched_length"] = new_blast_all_f2["qend"] - new_blast_all_filter1["qstart"] - new_blast_all_filter1["mismatch"]-new_blast_all_filter1["gapopen"]


geneID0 = new_blast_all_f2["geneID"].tolist()
geneID1 = new_blast_all_f2["hitID"].tolist()
geneID0_length = [gene_id_length_new_species[ID] for ID in geneID0]
geneID1_length = [gene_id_length_343[ID] for ID in geneID1]
# get the minimum lenght of the paired
minimum_length = [min([x,y]) for x,y in zip(geneID0_length, geneID1_length)]
new_blast_all_f2["minimum_length"] = minimum_length

new_blast_all_f2["percent_match"] = new_blast_all_f2["matched_length"]/new_blast_all_f2["minimum_length"]
new_blast_all_f3 = new_blast_all_f2[new_blast_all_f2["percent_match"] >= match_percent]
print(len(new_blast_all_f3["geneID"].unique())) # 0.6--2244; 0.5--3373 genes; 0.475--3618; 0.45--3874; 0.425--4111; 0.4--4320;0.35--4629; 0.3--4878; 0.2--5116; 0.1--5187


## then find the OGs based on the id mapping
# here we establish the connection between the ortholog ID and the number order in each ortholog ID
id_mapping_dir = "/Users/luho/Google Drive/yeast model update/300 yeast species/0_332yeast_genomes/"
infile = id_mapping_dir + 'orthomcl_output/orthomcl_clusters.txt'
all_id = open(infile, "r").readlines()
len(all_id)
# here we just change it into a dict
all_id_dict = {}
for element in all_id:
    s1 = element.split(': ')
    file_name = s1[0]
    a_id = s1[1].strip("\n")
    a_id0 = a_id.split(' ')
    all_id_dict[file_name] = a_id0

# establish one to one mapping between OG and indexes
new_dict = {}
new_value = []
for key, value in all_id_dict.items():
    print(key, value)
    value0 = [x +"&&"+ key for x in value]
    new_value.append(value0)


#list2d = [[1,2,3], [4,5,6], [7], [8,9]]
#merged = list(itertools.chain(*list2d))
new_value1 = list(itertools.chain(*new_value)) # flattern the list
#new_value1 =  sum(new_value, []) # flattern the list
new_dict ={s.split("&&")[0]:s.split("&&")[1] for s in new_value1}



# here we connect the number order (connect with ortholog ID) with the protein ID
infile2 = id_mapping_dir  + 'orthomcl_output/orthomcl_SeqIDs_index.txt'
id_mapping1 = open(infile2, "r").readlines()
len(id_mapping1)
id_map_dict1 = dict()
for id in id_mapping1:
    print(id)
    id = id.strip('\n')
    id0 = id.split(': ')
    id_map_dict1[id0[1]] = id0[0]

id_map_dict2 ={x:new_dict[y] for x, y in id_map_dict1.items()}




## find the OG id based on the hit id
interestID = new_blast_all_f3["hitID"].tolist()
OGID = [id_map_dict2[x] for x in interestID]
new_blast_all_f3["OGID"] = OGID
geneID_with_target = new_blast_all_f3["geneID"]
# check the mapping
geneID_OGID_dict = {}
combine0 = [x + "&&" + y for x, y in zip(geneID_with_target, OGID)]
new_blast_all_f3["combine2"] = combine0

combine0_unique = list(set(combine0)) # 4102 unique pair
new_df = pd.DataFrame({"combine":combine0_unique})
new_df2 = pd.DataFrame(new_df["combine"].str.split('&&',1).tolist(),
                                 columns = ['geneID','OGID'])

new_df3 = pd.DataFrame({"geneID": new_df2["geneID"].unique()})
new_df3["OGID"] = multiMapping(new_df2["OGID"], new_df2["geneID"], new_df3["geneID"])
print(len(new_df3["geneID"]))



# It is found that one gene could be mapped onto several OGs
# Thus need to further filter OGs? or not just keep it?
# strategy to do filteration
# 1) max bitscore
# 2) max match percent
from collections import Counter
frequence_pair =Counter(combine0)

OG_unique = []
for x, y in new_df3.iterrows():
    #print(y["geneID"])
    gene00= y["geneID"]
    print(gene00)
    og00= y["OGID"].split(";")
    #print(og00)
    # test
    # gene00 = "XP_034013157.1"
    sub_df = new_blast_all_f3[new_blast_all_f3["geneID"]==gene00]
    # choose the maximum
    bit_score_max = max(sub_df["bitscore"])
    sub_df2 = sub_df[sub_df["bitscore"]==bit_score_max]
    og_optional = sub_df2["OGID"].unique().tolist()[0]
    if len(sub_df2["OGID"].unique().tolist()) >=2:
        match_percent_max = max(sub_df2["percent_match"])
        sub_df3 = sub_df2[sub_df2["percent_match"] == match_percent_max]
        og_optional = sub_df3["OGID"].unique().tolist()[0]
    else:
        og_optional = sub_df2["OGID"].unique().tolist()[0]
    OG_unique.append(og_optional)

new_df3["optional_OG"]= OG_unique
new_df3.to_csv(output_file)

