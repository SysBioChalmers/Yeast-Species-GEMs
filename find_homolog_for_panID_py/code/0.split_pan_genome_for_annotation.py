# -*- coding: utf-8 -*-
# -*- python 3 -*-
# -*- hongzhong Lu -*-

'''
The pan genome annotation could be divided into the followed three steps
1. EggNog annotation
2. KEGG annotation
3. RAVEN2 based on biocyc database

Once the annotation is obtained, a find_homolog_for_panID_py will be obtained.
'''
import os    ##for directory
import sys

# set the directory
sys.path.append(r"/Users/luho/PycharmProjects/3D_model/find_homolog_for_panID_py/code")
os.chdir('/Users/luho/PycharmProjects/3D_model/find_homolog_for_panID_py/code')
from mainFunction import *

# input the fast file
infile = '../data/representatives.fasta'
outfile = '../result/'


# divided the fasta file into three part for the function annotation from EggNOG
i = 0
j = 0
seq = []
ID = []
with open(infile) as infile0:
  for line in infile0:
      i = i +1
      print(line)
      seq.append(line)
      if line.startswith('>'):
         ID.append(line)
      else:
          continue

# estimate the separate line index
num1 = seq.index(ID[80000])-1
num2 = seq.index(ID[160000])-1


# group1 1-8000
# group2 8001-16000
# group3 16001-233478
g1 = open('../result/fasta1_yeast.fasta', 'w')
g2 = open('../result/fasta2_yeast.fasta', 'w')
g3 = open('../result/fasta3_yeast.fasta', 'w')

for i,x in enumerate(seq):
    print(i, x)
    if i <= num1:
        g1.write(x)
    elif i> num2:
        g3.write(x)
    else:
        g2.write(x)

g1.close()
g2.close()
g3.close()


''' some seq need check
>yHMPu5000035268_Wickerhamomyces_hampshirensis@Seq_3870
MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM

Aspergillus_nidulans
'''



# divided the fasta file into the reference fasta from sce-s288c and other yeast species
# this is used to conduct the blast analysis between s88c proteins and pan-proteins
from Bio import SeqIO
infile = '../data/representatives.fasta'
records = list(SeqIO.parse(infile, "fasta"))
print(records[0].id)  # first record
print(records[-1].id)  # last record

non_ref_sequences = [] # Setup an empty list
for record in records:
    print(record)
    if 'Saccharomyces_cerevisiae@' not in record.id:
        # Add this record to our list
        non_ref_sequences.append(record)
#print("Found %i short sequences" % len(short_sequences))
#save the fasta sequences
SeqIO.write(non_ref_sequences, "../result/non_ref_sequences.fasta", "fasta")



'''
# check the origial 332 yeast strain genome annotation data
from Bio import SeqIO
infile1 = '../data/saccharomyces_arboricola.max.pep'
records = list(SeqIO.parse(infile1, "fasta"))
print(records[0].id)  # first record
print(records[-1].id)  # last record

non_ref_sequences = [] # Setup an empty list
for record in records:
    print(record)
'''


# further divide the fasta file into small part
# this is used for the deepec, which can predict the EC number from protein fasta file
# due to the memory limitation, we need split it into 6 small files
from Bio import SeqIO

infile = '../data/representatives.fasta'
records = list(SeqIO.parse(infile, "fasta"))
print(records[0].id)  # first record
print(records[-1].id)  # last record
out = ['p1','p2','p3','p4','p5','p6']
start0 = [0,40000,80000,120000,160000, 200000]
end0 = [40000,80000,120000,160000, 200000, 233478]

for i in range(len(out)):
    panYeast_part1 = []  # Setup an empty list
    outfile = "../result/panYeast_" + out[i] + ".fasta"
    for id in range(start0[i], end0[i]):
        print(id)
        panYeast_part1.append(records[id])
    SeqIO.write(panYeast_part1, outfile, "fasta")


'''
# another small tasks
from Bio import SeqIO

infile = '../data/protein_sgd.fasta'
records = list(SeqIO.parse(infile, "fasta"))
print(records[0].id)  # first record
print(records[-1].id)  # last record
out = ['p1','p2','p3','p4','p5','p6']
start0 = [0,1000,2000,3000,4000, 5000]
end0 = [1000,2000,3000,4000, 5000, 6713]

for i in range(len(out)):
    panYeast_part1 = []  # Setup an empty list
    outfile = "../result/sce_sgd_" + out[i] + ".fasta"
    for id in range(start0[i], end0[i]):
        print(id)
        panYeast_part1.append(records[id])
    SeqIO.write(panYeast_part1, outfile, "fasta")
'''