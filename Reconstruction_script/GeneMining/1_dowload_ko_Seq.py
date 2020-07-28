#!/usr/bin/env python
# pathway_genes.py
# adapted from https://gist.github.com/slowkow/a2327b868f4e927ac4fb
# Feiran Li 2020.02.23

import urllib
import urllib.request
from bs4 import BeautifulSoup as bs
import sys
import random
import pandas as pd
import os


def main(checksub): # True will only extract ko sequences from species with this phenotype
    if checksub:
        sub_strain = sub_yeast_kegg() # read the data of subtrate usage info for all KEGG yeast
    unmappedko = []
    # args = sys.argv
    # if len(args) != 2:
    #     print
    #     'Usage: python 1_dowload_ko_Seq.py rxnlist.txt'
    #     sys.exit(1)
    f2 = open("rxnlist.txt","r")
    lines = f2.readlines()
    for rxn_id in lines[1:]:
        rxn_id = rxn_id.replace("\n", "")
        rxn_id = rxn_id.split(',', 1)
        sub = rxn_id[0]
        rxn_id = rxn_id[1]
        orthology_ids = get_orthology_ids(rxn_id)
        print
        'Found {} orthology ids for rxn "{}"' \
            .format(len(orthology_ids), rxn_id)

        if orthology_ids:
            for orthology_id in orthology_ids:
                gene_ids = get_gene_ids(orthology_id)

                print
                'Writing {} FASTA gene sequences to "{}.fa"' \
                    .format(len(gene_ids), orthology_id)
                # fungi ids, in order to speed it up, we choose only 20 sequences here
                fungi = []
                with open('KEGG_org_yeast.txt', 'r') as f: #can also be KEGG_ORG_FUNGI.TXT
                    for line in f:
                        fungi.append(line.strip('\n')) # all fungi ids
                x = []
                if checksub:# whether we input the sub_strain
                    x = [name for name in gene_ids if name.startswith(tuple(sub_strain[sub]))]
                if not x:
                # find wthether there is a fungi gene included in the gene_ids
                    x = [name for name in gene_ids if name.startswith(tuple(fungi))]
                    if len(gene_ids) > 5:
                        if x is not None:

                            if len(x) > 5:
                                gene_ids = random.sample(gene_ids, 5) + random.sample(x, 5) # can slao be 20
                            else:
                                gene_ids = random.sample(gene_ids, 5) + x
                            gene_ids = list(set(gene_ids))
                        else:
                            gene_ids = random.sample(gene_ids, 5)
                else:
                    gene_ids = x
                with open("ko/" + orthology_id + "_" + rxn_id + ".fa", 'w') as out:
                    for i, gene_id in enumerate(gene_ids, 1):
                        sys.stdout.write('.')
                        if not i % 5:
                            sys.stdout.write(' ')
                        if not i % 50:
                            sys.stdout.write('\n')
                        sys.stdout.flush()

                        fasta = get_fasta(gene_id)
                        out.write(fasta)
        else:
            unmappedko.append(rxn_id)
    return(unmappedko)

def sub_yeast_kegg():
    # load the yeast(in kegg) substrate usage info, in order to filter out the uncertainty
    subusage = pd.read_csv('Subusage_KEGG_yeast.tsv', sep='\t', header=0) # this one load phenotypes realted species, s that it can be used to filter the ko sequences
    sub_strain = {}
    for sub in subusage.keys()[2:]: # take out the first two whcih are not actual substrate: 0 kegg_ID 1 substrate
        idx = subusage[(subusage[sub] == '1') | (subusage[sub] == 'v')].index.tolist()
        temp = subusage['KEGG_ID'][idx]
        sub_strain[sub] = temp.values.tolist()
    return sub_strain

def get_ids(url):
    response = urllib.request.urlopen(url)
    html = response.read()
    b = bs(html)
    links = b.select("a[href*=www_bget]")
    return [link.text for link in links]


def get_orthology_ids(rxn_id):
    URL = 'http://www.genome.jp'
    FUN = '/dbget-bin/get_linkdb?-t+orthology+rn:'
    return get_ids(URL + FUN + rxn_id)


def get_gene_ids(orthology_id):
    URL = 'http://www.genome.jp'
    FUN = '/dbget-bin/get_linkdb?-t+genes+ko:'
    return get_ids(URL + FUN + orthology_id)


def get_fasta(gene_id):
    URL = 'http://www.genome.jp'
    FUN = '/dbget-bin/www_bget?-f+-n+a+' # -n+a for aa seq -n+n for nucleotide seq
    #FUN2 = 'dbget-bin/www_bget?-f+refseq+'
    response = urllib.request.urlopen(URL + FUN + gene_id)
    html = bs(response.read(),features="html.parser")
    return html.pre.text



if __name__ == '__main__':
    main(True)


