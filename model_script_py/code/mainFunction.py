# -*- coding: utf-8 -*-
# -*- python 3 -*-
# -*- hongzhong Lu -*-


# Import packages
import re
import numpy as np
import pandas as pd
import os    ##for directory
import sys
import pprint


'''general function for easy use of python'''
def splitAndCombine(gene, rxn, sep0, moveDuplicate=False):
    ## one rxn has several genes, this function was used to splite the genes
    ## used for the dataframe data

    gene = gene.fillna('NA')  # fill the NaN with 'NA'
    gene0 = gene.tolist()
    rxn0 = rxn.tolist()
    s1 = list()
    s2 = list()
    for i in range(len(gene0)):
        s1 = s1 + [rxn0[i]] * len(gene0[i].split(sep0))
        s2 = s2 + gene0[i].split(sep0)
    df0 = pd.DataFrame({'V1': s1,
                        'V2': s2}
                       )
    if moveDuplicate == True:
        df00 = df0.drop_duplicates()
    else:
        df00 = df0
    return df00


def getSimilarTarget(rxn_yeast0,rxn_newGPR0,ss):
    from fuzzywuzzy import fuzz
    from fuzzywuzzy import process
    rxn_yeast1 = np.array(rxn_yeast0)  # np.ndarray()
    rxn_yeast2 = rxn_yeast1.tolist()
    rxn_yeast3 = pd.Series((v[0] for v in rxn_yeast2))
    rxn_newGPR1 = np.array(rxn_newGPR0)  # np.ndarray()
    rxn_newGPR2 = rxn_newGPR1.tolist()
    rxn_newGPR3 = pd.Series((v[0] for v in rxn_newGPR2))
    similarTarget = [None] * ss
    for i in range(ss):
        similarTarget[i] = process.extract(rxn_newGPR3[i], rxn_yeast3, limit=2)

    return similarTarget
'''
#example
newMet = pd.read_excel('new metabolite for check.xlsx')
newMet0 = newMet[['name_unify']]
gemMet = pd.read_excel('unique metabolite in yeastGEM.xlsx')
gemMet0 = gemMet[['Description_simple']]
ss0 = len(newMet0)
similarTarget0 = getSimilarTarget(gemMet0,newMet0,ss=ss0)
'''


def singleMapping (description, item1, item2, dataframe=True):
    """get the single description of from item1 for item2 based on mapping"""
    #description = w
    #item1 = v
    #item2 = testData
    # used for the list data
    if dataframe:
        description = description.tolist()
        item1 = item1.tolist()
        item2 = item2.tolist()
    else:
        pass
    index = [None]*len(item2)
    result = [None]*len(item2)
    tt = [None]*len(item2)
    for i in range(len(item2)):
        if item2[i] in item1:
            index[i] = item1.index(item2[i])
            result[i] = description[index[i]]
        else:
            index[i] = None
            result[i] = None
    return result
'''
w=['a','b','c']
v=[1,2,3]
s=[3,1,2,4]
singleMapping(w,v,s,dataframe=False)
'''


def multiMapping (description, item1, item2, dataframe=True, sep=";", removeDuplicates=True):
    """get multiple description of from item1 for item2 based on mapping"""
    #description = w
    #item1 = v
    #item2 = testData
    #used for the list data
    if dataframe:
        description = description.tolist()
        item1 = item1.tolist()
        item2 = item2.tolist()
    else:
        pass
    result = [None]*len(item2)
    for i in range(len(item2)):
        if item2[i] in item1:
            index0 = [description[index] for index in range(len(item1)) if item1[index] == item2[i]]
            if removeDuplicates:
                index1 = pd.unique(index0).tolist()
            else:
                index1 = index0
            result[i] = sep.join(str(e) for e in index1) #string cat
        else:
            result[i] = None
    return result

'''
# example data to test all the above function
df1 = pd.DataFrame({'A' : ['one', 'one', 'two', 'three'] * 3,
                    'B' : ['A', 'B', 'C'] * 4,
                    'C' : ['foo', 'foo', 'foo', 'bar', 'bar', 'bar'] * 2}
                   )

df2 = pd.DataFrame({'A' : ['one', 'one', 'two', 'three'] * 3,
                    'B' : ['A', 'B', 'C'] * 4,
                    'D' : np.random.randn(12)})


df2['C'] = singleMapping(df1['C'], df1['A'], df2['A'])
df2['C'] = multiMapping(df1['C'], df1['A'], df2['A'])
'''


def updateOneColumn(df1, df2, key0, value0):
    """
    using dataframe df2 to update the df1

    :param df1:
    :param df2:
    :param key0: the common column name, a string, used for the mapping
    :param value0: the column in df2 used to update the df1
    :return:
    example
    df10 = pd.DataFrame({'A': ['a', 'b', 'c'],
                 'B': ['x', 'y', 'z']})

    df20 = pd.DataFrame({'A':['c','b'],
                       'B': ['e', 'd']})
    updateOneColumn(df10,df20,key0='A',value0='B')
    """
    df10 = df1.copy()
    df11 = df1.copy()
    df10[value0] = multiMapping(df2[value0], df2[key0], df10[key0])
    for i, x in df10.iterrows():
        print(x[value0])
        if x[value0] is None:
            df11[value0][i] = df11[value0][i]
        else:
            df11[value0][i] = df10[value0][i]
    return df11[value0]






def RemoveDuplicated(s1):
    """
    example:
    s1=['a // a', 'b // a', None, 'non']

    """
    s2=list()
    for x in s1:
        print(x)
        if x =='non':
            s2.append('')
        elif x is None:
            s2.append('')
        else:
            if "//" in x:
                s0= x.split(' // ')
                s0 = [x.strip() for x in s0]
                s01= list(set(s0))
                if len(s01)==1:
                    s2.append(s01[0])
                else:
                    s2.append(' // '.join(s01))
            else:
                s2.append(x)
    return s2



def nz(value):

    '''
    Convert None to string else return value.
    '''

    if value == None:
        return 'none'
    return value



def AutoUpdate(description1, para1, description2, para2):
    # using the description1 in para1 to update the description2 in para2
    description1 = description1.tolist()
    para1 = para1.tolist()
    description2 = description2.tolist()
    para2 = para2.tolist()
    ss = [None]*len(para2)
    for i in range(len(para2)):
       if para2[i] in para1:
          ss[i] = para1.index(para2[i])
       else:
          ss[i] = None

    for i in range(len(para2)):
        if ss[i] != None:
            description2[i] = description1[ss[i]]
        else:
            description2[i] = description2[i]

    return description2



'''
# example data to test the followed function
df1 = pd.DataFrame({'A' : ['one', 'one', 'two', 'three'] * 3,
                    'B' : ['A', 'B', 'C'] * 4,
                    'C' : ['foo', 'foo', 'foo', 'bar', 'bar', 'bar'] * 2}
                   )

df2 = df1.iloc[[1,2]]
df2['C'] = ['good','good']
df1['C'] = AutoUpdate(df2['C'],df2['A'],df1['C'],df1['A'])
'''


def calculateFrequency(list0, item0):
    '''
    This function is used to calculate the frequency occured in a list and turn the frequency list into a dataframe
    :param list0:  ['a','b','a']
    :param item0:
    :return: a dataframe with two columns
    '''
    summary = pd.Series(list0).value_counts()
    summary = summary.to_frame(name='number')
    summary.index.name = item0
    summary.reset_index(inplace=True)
    return summary






"""function for model part"""
from cobra.manipulation import remove_genes
def getStrainGEMrxn(s0, geneMatrix0, templateGEM, templateGene):
    '''
    This function is used to produce the strain specific model based on panYeast and gene existence matrix
    from 1011 yeast strain genome sequence project
    :param s0: strain name 'BFC'
    :param geneMatrix0: dataframe contains the gene existence matrix for each strain. geneMatrix = pd.read_csv('../data/geneMatrix0 of 1011 yeast strains.txt', sep="\t")
    :templateGEM:
    :templateGene:
    :return: the rxn list for each new reaction3

    '''

    s1 = ['geneID', s0]
    geneList = geneMatrix0.loc[:, s1]
    gene_exist = singleMapping(geneList.loc[:, s0].tolist(), geneList.loc[:, 'geneID'].tolist(), templateGene,
                               dataframe=False)
    gene_exist = [0 if v is None else v for v in gene_exist]
    gene_remove = [x for x, y in zip(templateGene, gene_exist) if y < 1]
    newModel = templateGEM.copy()
    # for i in range(len(gene_remove)):
    #    print(i)
    #    remove_genes(newModel, [gene_remove[i]], remove_reactions=True)
    remove_genes(newModel, gene_remove, remove_reactions=True)
    rxn = []
    for x in newModel.reactions:
        rxn.append(x.id)
    return rxn

def getStrainGEM(s0, geneMatrix0, templateGEM, templateGene):
    '''
    This function is used to produce the strain specific model based on panYeast and gene existence matrix
    from 1011 yeast strain genome sequence project
    :param s0: strain name 'BFC'
    :param geneMatrix0: dataframe contains the gene existence matrix for each strain. geneMatrix = pd.read_csv('../data/geneMatrix0 of 1011 yeast strains.txt', sep="\t")
    :templateGEM:
    :templateGene:
    :return: the rxn list for each new reaction3

    '''

    s1 = ['geneID', s0]
    geneList = geneMatrix0.loc[:, s1]
    gene_exist = singleMapping(geneList.loc[:, s0].tolist(), geneList.loc[:, 'geneID'].tolist(), templateGene,
                               dataframe=False)
    gene_exist = [0 if v is None else v for v in gene_exist]
    gene_remove = [x for x, y in zip(templateGene, gene_exist) if y < 1]
    newModel = templateGEM.copy()
    # for i in range(len(gene_remove)):
    #    print(i)
    #    remove_genes(newModel, [gene_remove[i]], remove_reactions=True)
    remove_genes(newModel, gene_remove, remove_reactions=True)
    return newModel







def getRemoveGeneList(s0, geneMatrix0, templateGEM, templateGene):
    '''
    This function is used to produce the strain specific model based on panYeast and gene existence matrix
    from 1011 yeast strain genome sequence project
    :param s0: strain name 'BFC'
    :param geneMatrix0: dataframe contains the gene existence matrix for each strain. geneMatrix = pd.read_csv('../data/geneMatrix0 of 1011 yeast strains.txt', sep="\t")
    :templateGEM:
    :templateGene:
    :return: the gene list removed from each strain specific model

    '''

    s1 = ['geneID', s0]
    geneList = geneMatrix0.loc[:, s1]
    gene_exist = singleMapping(geneList.loc[:, s0].tolist(), geneList.loc[:, 'geneID'].tolist(), templateGene,
                               dataframe=False)
    gene_exist = [0 if v is None else v for v in gene_exist]
    gene_remove = [x for x, y in zip(templateGene, gene_exist) if y < 1]
    newModel = templateGEM.copy()
    # for i in range(len(gene_remove)):
    #    print(i)
    #    remove_genes(newModel, [gene_remove[i]], remove_reactions=True)
    remove_genes(newModel, gene_remove, remove_reactions=True)
    gene = []
    for x in newModel.genes:
        gene.append(x.id)

    gene_remove_from_model = list(set(templateGene)-set(gene))
    return  gene_remove_from_model




def updateGPR(gpr0, nameMapping):
    '''
    This function is used to update the gpr reaction only with 'or' relation. It is used to replace the old gene name using
    the new gene name. Also it did not remove the duplicated value.
    :param: gpr0
    :nameMapping: a dataframe contains the mapping relation between the old and new gene name, has two columns-'geneID', 'panID'
    :return: gpr with the replaced new gene name
    '''
    #this function is mainly used to update the gene relation with 'or'
    s1 = gpr0
    s2 = s1.split(' ')
    s3 = singleMapping(nameMapping['panID'].tolist(),nameMapping['geneID'].tolist(),s2, dataframe=False)
    for i, x in enumerate(s3):
        if x is None:
            s3[i]=s2[i]
        else:
            s3[i] = s3[i]
    s4 = ' '.join(s3)
    return s4



def getCompartment(rxn):
    """
    This function is used to obtain the compartment information from reaction of yeastGEM
    :param rxn:  example acetyl-CoA[m] + L-glutamate[m]  -> coenzyme A[m] + H+[m] + N-acetyl-L-glutamate[m]'
    :return:
    """
    cp1 = ['[c]','[ce]','[e]','[er]','[erm]','[g]','[gm]','[lp]','[m]','[mm]','[n]','[p]','[v]','[vm]']
    cp2 = ['cytoplasm','cell envelope','extracellular','endoplasmic reticulum','endoplasmic reticulum membrane','Golgi','Golgi membrane','lipid particle',
             'mitochondrion','mitochondrial membrane','nucleus','peroxisome','vacuole','vacuolar membrane']


    cp = [None]*len(cp1)
    for i in range(len(cp1)):
       if cp1[i] in rxn:
         cp[i] = cp2[i]
       else:
          cp[i] = None
    cp1 = [x for i,x in enumerate(cp) if x is not None]
    cp0 = ';'.join(str(e) for e in cp1)
    return cp0



def getCommonCompartment(c1,c2, sep0=";"):
    '''this function could get the common part between string c1 and c2
    for example, c1="a;b", c2="a;c" '''
    if c1 is None:
        c10 = 'NONE'
    else:
        c10 = c1.split(sep0)
        c10 = [x.strip() for x in c10]
    if c2 is None:
        c20 = 'NONE'
    else:
        c20 = c2.split(sep0)
        c20 = [x.strip() for x in c20]
    c3 = list(set(c10).intersection(c20))
    c4 = sep0.join(str(e) for e in c3)
    return c4


def getRXNgeneMapping(rxn0, gpr0):
    '''this function is used to split the GPR;
    input, for example rxn0=['r1','g2']
    gpr0=['a or c','a and b']
    output, each rxn related with each gene'''
    s1 = rxn0
    s2 = gpr0
    s2 = s2.str.replace('and','@')
    s2 = s2.str.replace('or','@')
    s2 = s2.str.replace('\\( ','')
    s2 = s2.str.replace('\\(\\( ','')
    s2 = s2.str.replace('\\(', '')
    s2 = s2.str.replace('\\(\\(', '')
    s2 = s2.str.replace(' \\)','')
    s2 = s2.str.replace(' \\)\\) ','')
    s2 = s2.str.replace('\\)', '')
    s2 = s2.str.replace('\\)\\) ', '')
    s3 = splitAndCombine(s2,s1,sep0="@")
    s3['V2'] = s3['V2'].str.strip()
    s3.columns = ['rxnID', 'gene']
    return s3

def getRXNmetaboliteMapping(rxn0, met0):
    '''this function is used to split the equation of metabolites; used to produce the dataframe format of GEM using
    cobrapy
    input, for example rxn0=['r1','g2']
    gpr0=['a => c','a => b']
    output, each rxn related with each gene'''
    met_annotation = pd.read_excel('/Users/luho/PycharmProjects/model/cobrapy/result/met_yeastGEM.xlsx')
    s1 = rxn0
    s2 = met0
    s3 = splitAndCombine(s2,s1,sep0=" ")
    s3['V2'] = s3['V2'].str.strip()
    s3.columns = ['rxnID', 'met']
    s3['met_name'] = singleMapping(met_annotation['description'],met_annotation['m_name'],s3['met'])
    for i, x in s3.iterrows():
        if s3['met_name'][i] is None:
            s3['met_name'][i] = s3['met'][i]
        else:
            s3['met_name'][i] = s3['met_name'][i]
    return s3


def correctSomeWrongFormat(model0):
  """
  This function is used to correct some wrong format when read yeastGEM model from cobratoolbox
  """
  # Correct metabolite ids:
  for met in model0.metabolites:
      met.id = met.id.replace('__93__', '')
      met._id = met._id.replace('__91__', '_')
      print(met.id)
  for reaction in model0.reactions:
      reaction.gene_reaction_rule = reaction.gene_reaction_rule.replace('__45__', '-')
      print(reaction.gene_reaction_rule)
  for gene in model0.genes:
      gene.id = gene.id.replace('__45__', '-')

  return model0


def produceMetaboliteList(model0):
  #produce the dataframe for the metabolites from yeastGEM
  met_list = [None] * len(model0.metabolites)
  met_dataframe = pd.DataFrame({'m_name': met_list,
                                'description': met_list,
                                'formula': met_list,
                                'charge': met_list,
                                'chebi': met_list,
                                'kegg': met_list,
                                'MNXID': met_list})

  for i, met in enumerate(model0.metabolites):
      print(i)
      met_dataframe['m_name'][i] = met.id
      met_dataframe['description'][i] = met.name
      met_dataframe['formula'][i] = met.formula
      met_dataframe['charge'][i] = met.charge
      key = list(met.annotation.keys())
      if 'chebi' in key:
          met_dataframe['chebi'][i] = met.annotation['chebi']
      else:
          met_dataframe['chebi'][i] = None
      if 'kegg.compound' in key:
          met_dataframe['kegg'][i] = met.annotation['kegg.compound']
      else:
          met_dataframe['kegg'][i] = None
      if 'metanetx.chemical' in key:
          met_dataframe['MNXID'][i] = met.annotation['metanetx.chemical']
      else:
          met_dataframe['MNXID'][i] = None

  #s2 = met_dataframe['m_name'].str.split('_', expand=True)
  #met_dataframe['description'] = met_dataframe['description'].str.replace('\s\[', '@')
  #s3 = met_dataframe['description'].str.split('@', expand=True)
  #met_dataframe['description'] = s3.iloc[:, 0] + '[' + s2.iloc[:, 2] + ']'
  return met_dataframe


def produceGeneList(model0):
  #produce the gene list from GEM
  genelist = []
  for i in model0.genes:
      print(i)
      genelist.append(i.id)
  return genelist


def produceRxnList(model0):
  #produce the dataframe for the rxn from yeastGEM
  reaction_list =[None]*len(model0.reactions)
  gem_dataframe = pd.DataFrame({'name':reaction_list,
                              'equation':reaction_list,
                              'GPR':reaction_list,
                              'rxnID':reaction_list,
                              'formula':reaction_list
                              })

  for i, reaction in enumerate(model0.reactions):
      print(i)
      gem_dataframe['name'][i] = reaction.name
      gem_dataframe['equation'][i] = reaction.reaction
      gem_dataframe['GPR'][i] = reaction.gene_reaction_rule
      gem_dataframe['rxnID'][i] = reaction.id
  gem_dataframe['ID'] = ['R'+ str(i) for i in range(0, len(model0.reactions))]
  gem_dataframe['GPR'] = gem_dataframe['GPR'].str.replace('__45__', '-')
  #replace the metabolite name in gem_dataframe
  s0 = getRXNmetaboliteMapping(gem_dataframe['rxnID'], gem_dataframe['equation'])
  gem_dataframe['formula'] = multiMapping(s0['met_name'],s0['rxnID'],gem_dataframe['rxnID'],removeDuplicates=False)
  gem_dataframe['formula'] = gem_dataframe['formula'].str.replace(";", " ")

  return gem_dataframe



def exchange(s1, subystem):
    """
    this function is used to define the exchange reaction
    s1=['a --> b','a <=> c', 'H+ [extracellular] + L-citrulline [extracellular] <=> H+ [cytoplasm] L-citrulline [cytoplasm]', ' a--> ']
    subsystem = ['a','a','b','']

    """
    for i, x in enumerate(s1):
        print(i)
        if ' --> ' in x:
            x0 = x.split(' --> ')
            if len(x0[1]) >=1:
                #subystem.append('General')  # exchange
                subystem[i] = subystem[i]
            else:
                subystem[i] ='Exchange reaction' #exchange
                print(subystem[i])
        if ' <=> ' in x:
            x0 = x.split(' <=> ')
            if len(x0[1]) >=1:
                #subystem.append('General')  # exchange
                subystem[i] = subystem[i]
            else:
                subystem[i] ='Exchange reaction' #exchange
                print(subystem[i])
        else:
            subystem[i] = subystem[i]
    return subystem


def exchange_ecYeast(s1, subystem):
    """
    this function is used to define the exchange reaction
    s1=['a --> b','a <=> c', 'H+ [extracellular] + L-citrulline [extracellular] <=> H+ [cytoplasm] L-citrulline [cytoplasm]', ' a--> ']
    subsystem = ['a','a','b','']

    """
    for i, x in enumerate(s1):
        print(i)
        if ' --> ' in x:
            x0 = x.split(' --> ')
            if len(x0[1]) >=1 and len(x0[0]) >=1:
                #subystem.append('General')  # exchange
                subystem[i] = subystem[i]
            else:
                subystem[i] ='Exchange reaction' #exchange
                print(subystem[i])
        if ' <=> ' in x:
            x0 = x.split(' <=> ')
            if len(x0[1]) >=1 and len(x0[0]) >=1:
                #subystem.append('General')  # exchange
                subystem[i] = subystem[i]
            else:
                subystem[i] ='Exchange reaction' #exchange
                print(subystem[i])
        else:
            subystem[i] = subystem[i]
    return subystem


#SLIME rxn
def SLIME(rxnName, subsystem):
    """
    if the rxnName contains the SLIME, classify the reaction into SLIME reaction
    """
    for i,x in enumerate(rxnName):
        if 'SLIME' in x:
            subsystem[i] = 'SLIME reaction'
            print(subsystem[i])
        else:
            subsystem[i] = subsystem[i]
    return subsystem



def transport(s1, subsysem):
    """
    this function is used to define the transport reaction
    #example
     s1 =['2-methylbutyl acetate [cytoplasm] --> 2-methylbutyl acetate [extracellular]', 'H+ [extracellular] + phosphoenolpyruvate [extracellular] <=> H+ [cytoplasm] + phosphoenolpyruvate [cytoplasm]']
     subsysem = ['a','b']

    :param s1:
    :param subsysem:
    :return:
    """
    for i, x0 in enumerate(s1):
        x1 = re.findall(r"\[([A-Za-z0-9_\s]+)\]", x0)
        x0 = x0.replace('(','[')
        x0 = x0.replace(')',']')
        x2 = re.sub(r"\[([A-Za-z0-9_\s+]+)\]", '', x0)

        if "<=>" in x2:
            x3 = x2.split("<=>")
        elif "<->" in x2: #bigg database format
            x3 = x2.split("<->")
        else:
            x3 = x2.split("-->")
        x3 = [x.strip() for x in x3]
        x1=pd.unique(x1).tolist() #remove the duplicated
        if '+' in x3[0]:
            x30=x3[0].split('+')
        else:
            x30=x3[0]
        x30=[x.strip() for x in x30]
        x30 = [x for x in x30 if x != '']
        if '+' in x3[1]:
            x31 = x3[1].split('+')
        else:
            x31=x3[1]
        x31 = [x.strip() for x in x31]
        x31 = [x for x in x31 if x != '']

        if set(x30) == set(x31):
            subsysem[i] ='Transport' + '['+', '.join(x1)+']'
            print(subsysem[i])
        elif set(x30)-set(['ATP','H2O']) == set(x31) - set(['ADP','phosphate','H']):
            subsysem[i] = 'Transport' + '[' + ', '.join(x1) + ']'
            print(subsysem[i])
        else:
            subsysem[i] = subsysem[i]
    return subsysem



def findRemoveRxnBasedOnGene(rxnRemovedGene, rxnAllGene):
    '''this function is used to remove rxn based on the removed gene list
    if the all genes in a reaction were in the removed gene list, then this reaction was removed'''
    #x0 = gem_dataframe['removed_gene'].tolist()
    #y0 = gem_dataframe['all_gene'].tolist()
    x0=rxnRemovedGene.tolist()
    y0=rxnAllGene.tolist()
    removed_rxn = list()
    for x,y in zip(x0,y0):
        if x is None:
            removed_rxn.append('NO')
        else:
            if len(x) ==len(y):
                removed_rxn.append('YES')
            else:
                removed_rxn.append('NO')
    return removed_rxn



def saveExcel(infile, outfile):
    '''
    function to save the dataframe into xlsx format
    :param infile:
    :param outfile:
    :return:
    '''
    writer = pd.ExcelWriter(outfile)
    infile.to_excel(writer,'Sheet1')
    writer.save()


def find(rxn_name_list, specific_string, equal=False):
    '''
    function to find the index of element in a list contains a specific string, refer the function
    from matlab
    :param rxn_name_list:
    :param specific_string:
    :param equal: if true will return the indexes in rxn_name_list which is equal to specific_string
    :return: index of the element from the list contains the specific string
    example:
    rxn_name_list = ['a_b','e_b','c','d']
    specific_string = '_b'
    '''
    s0 = rxn_name_list
    if equal ==False:
        index =[i for i,x in enumerate(s0) if specific_string in x]
    else:
        index = [i for i, x in enumerate(s0) if specific_string == x]
    return index


'''model: function for 332 yeast species'''
def bbhFilterWithPident(a_to_b, b_to_a, pident):
    '''
    This function is used to conduct the standard BBH analysis based on the blast result from DIAMOND
    :param a_to_b: blast a to b using DIAMOND
    :param b_to_a: blast b to a using DIAMOND
    :param pident: a number
    :return:
    '''
    b1 = a_to_b.copy(deep=True)
    b2 = b_to_a.copy(deep=True)
    ss0 = ["geneID", "hitID", "pident", "length", "mismatch", "gapopen",
           "qstart", "qend", "sstart", "send", "evalues", "bitscore"]
    b1.columns = ss0
    b1['combine'] = b1["geneID"] + '@@@' + b1["hitID"]
    b2.columns = ss0
    b2['combine'] = b2["hitID"] + '@@@' + b2["geneID"]
    b1 = b1[b1['pident'] >= pident]
    b2 = b2[b2['pident'] >= pident]
    result_BBH0 = b1[b1['combine'].isin(b2['combine'])]
    return result_BBH0


def findBestHitFromBlast(panID, blast_inf):
    '''
    This function is used to find the best hit for homolog genes in non reference genomes from the reference genomes
    of s288c.
    Firstly find the hit with the highest pidentity.
    Then find the hit with the highest bitscore if the above result has multiple hits
    :param panID: A panID list
    :param blast_inf: A blast result from DIAMOND blast, which should have columns: 'geneID', 'hitID', 'pident','bitscore'
    :return:
    '''
    hit_best = []
    for x in panID:
        x_detail = blast_inf[blast_inf['geneID'] == x]
        pidentity = x_detail['pident'].tolist()
        hit = x_detail['hitID'].tolist()
        bitscore = x_detail['bitscore'].tolist()
        pi_max = max(pidentity)
        index_max_pi = [i for i, x in enumerate(pidentity) if x == pi_max]
        if (len(index_max_pi) >= 2):
            bitscore0 = [x for i, x in enumerate(bitscore) if i in index_max_pi]
            b_max = max(bitscore0)
            index_max_b = [i for i, x in enumerate(bitscore) if x == b_max]
            index = list(set(index_max_pi) & set(index_max_b))
            hitID = hit[index[0]]
        else:
            hitID = hit[index_max_pi[0]]
        hit_best.append(hitID)
    return hit_best




def singleBlastFilterWithPident(a_to_b, pident):
    '''
    This function is used to conduct the filter analysis based on a single blast result
    :param a_to_b: blast a to b using DIAMOND
    :param pident: a number
    :return:
    '''
    b1 = a_to_b.copy(deep=True)
    ss0 = ["geneID", "hitID", "pident", "length", "mismatch", "gapopen",
           "qstart", "qend", "sstart", "send", "evalues", "bitscore"]
    b1.columns = ss0
    b1['combine'] = b1["geneID"] + '@@@' + b1["hitID"]
    b1 = b1[b1['pident'] >= pident]
    return b1


def eggnogInput(dir_input):
    '''
    This function is used to input and summarize the eggnog annotation
    :param dir_input: The original data file from EggNOG database
    :return: A dataframe contains the detailed annotation of protein from EggNOG database
    '''
    f1_kegg = pd.read_fwf(dir_input, header =None)
    f1_kegg = f1_kegg[0].str.split('\t', expand=True)
    headers = f1_kegg.iloc[0, 0:17].tolist()
    other_col = ['annot lvl', 'matching OGs', 'Best OG', 'COG cat', 'description']
    headers0 = headers + other_col
    f1_kegg0 = f1_kegg[1:]
    f1_kegg0.columns = headers0
    return f1_kegg0


''' all simple plot'''
import matplotlib.pyplot as plt
import seaborn as sns

def simpleLineXY(x,y,x_title, y_title, title):
    '''
    This is a simple function to display the relation between x and y
    :param x:
    :param y:
    :param x_title:
    :param y_title:
    :param title:
    :return:
    '''
    fig=plt.figure(figsize=(4, 3))
    plt.plot(x, y, 'k', color='black')
    plt.xlabel(x_title)
    plt.ylabel(y_title)
    plt.title(title)
    plt.rc('xtick', labelsize=10)
    plt.rc('ytick', labelsize=10)
    plt.grid(linestyle='--', linewidth=0.5)
    plt.show()
    out = '../result/' + title + '.eps'
    fig.savefig(out, bbox_inches='tight', pad_inches=0)


def simpleBarPlot(x, y, x_title=None, y_title=None, title=''):
    '''
    This is a simple function to draw the bar plot
    :param x:
    :param y:
    :param x_title:
    :param y_title:
    :param title:
    :return:

    # example
    x = ['a','b']
    y = [19, 20]
    '''
    fig = plt.figure(figsize=(4, 3))
    y_pos = np.arange(len(x))
    plt.bar(y_pos, y, align='center', alpha=0.5)
    plt.xticks(y_pos, x)
    plt.xlabel(x_title)
    plt.ylabel(y_title)
    plt.title(title)
    plt.rc('xtick', labelsize=3)
    plt.rc('ytick', labelsize=3)
    plt.grid(linestyle='--', linewidth=0.5)
    plt.show()
    out = '../result/' + title + '.eps'
    fig.savefig(out, bbox_inches='tight', pad_inches=0)


def plotDensityProfile(ss, x_title, y_title):
    '''
    Function to calculate the frequency of element of ss (can be a series), then plot the density plot
    :param ss:
    :return:
    '''
    summary = ss.value_counts()
    summary0 = pd.DataFrame({'KO': summary.index, 'num': summary.values})
    # visualization
    plt.figure(figsize=(8, 6))
    ax = sns.distplot(summary0['num'])
    ax.set(xlabel=x_title, ylabel=y_title)
    ax.set(xlim=(0, 50))
    plt.show()