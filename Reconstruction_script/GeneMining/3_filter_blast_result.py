
import os
import glob
import numpy as np
import pandas as pd
import csv
from Bio import SeqIO
import re
import matplotlib.pyplot as plt
plt.rc('font',family='Helvetica')

# find the best hit
def filter_results(resultID,fungi,writefile,outfile) : # writefile:true or false, outfile:filename for the outfile

    if writefile:
        csv_writer = csv.writer(outfile)
        csv_writer.writerow(["query", "subject", "identity", "qcov",
                             "qlength", "sstart", "send",
                             "bitscore", "E-value", "strain","seq"])
    k = 0
    kkk = []
    total = pd.DataFrame(columns=('strain','qcov','identity','evalue'))

    for onefile in resultID:
        try:
            results = pd.read_csv(onefile, sep="\t", header=None)
            baseout = os.path.basename(onefile)[:-4]
            baseout = baseout.replace("blastres_","",1)
            # Add headers to forward and reverse results dataframes
            headers = ["query", "subject", "identity", "coverage",
                   "qlength", "slength", "alength","sstart","send",
                   "bitscore", "evalue"]
            results.columns = headers

            # Create a new column in both dataframes: normalised bitscore
            results['norm_bitscore'] = results.bitscore/results.qlength

            # Create query and subject coverage columns in both dataframes
            results['qcov'] = results.alength/results.qlength


            # Clip maximum coverage values at 1.0
            results['qcov'] = results['qcov'].clip(upper=1)

            idx = results[(results["evalue"] < 1e-20) & (results["qcov"] > 0.7)].index.tolist()
            if idx:
                result0 = results.loc[idx, :]
                s1 = result0.loc[:, "identity"] # find the best hit
                s1_argmax = s1[s1 == s1.max()].index
                s1_argmax = np.random.choice(s1_argmax)
                resultfinal = result0.loc[s1_argmax, :]
                # find strain
                temp = baseout.split('_', 2)
                ko = temp[0] + "_" + temp[1]
                strain = temp[2]
                # append the dataframe for final result
                total = total.append(pd.DataFrame({ 'strain': [strain], 'qcov': [resultfinal.qcov], 'identity': [resultfinal.identity], 'evalue': [resultfinal.evalue]}),ignore_index=True)
                query_species = resultfinal.query.split(':',1)[0]
                if query_species in fungi:
                    iden = 35
                else:
                    iden = 35
                if not (ko.startswith("K") & len(base) == 13): # which means not from kegg but the ID is from manually collected yeast species
                    iden = 35
                if resultfinal.identity > iden:
                    if writefile:
                        try:
                            # define location
                            location = []
                            location.append(resultfinal.subject)
                            location.append(resultfinal.sstart)
                            location.append(resultfinal.send)
                            # find sequence
                            current_path = os.getcwd()
                            datadir = current_path + "/db/"
                            temp = baseout.split("_", 2)
                            FASTA = temp[-1]
                            FASTA = datadir + FASTA + ".fas"
                            seq = seqfind(FASTA, location)
                            kkk.append(strain)
                            csv_writer.writerow([resultfinal.query,resultfinal.subject,resultfinal.identity,resultfinal.qcov,resultfinal.qlength,resultfinal.sstart,resultfinal.send,resultfinal.bitscore,resultfinal.evalue, baseout,seq])
                        except:
                            print(baseout)
                    else:
                        kkk.append(strain)
            else:
                k = k+1
                if k > 330:
                    print(k)
        except:
            print(onefile)
    outfile.close()
    return kkk, total

def seqfind(FASTA,location) :
    with open(FASTA, "rU") as handle:
        for record in SeqIO.parse(handle, "fasta") :
            geneID = str(record.id)
            searchid = location[0]
            if geneID in searchid:
                s = location[1]
                e = location[2]
                minnum = min(s, e)
                maxnum = max(s, e)
                seq = record.seq[minnum:maxnum]
                seq = str(seq)
    return seq



def find_rxn_species(resultlist,ko_rxn):
    ko_rxn_strain = {}
    j = 0
    for sub in list(ko_rxn.keys()):
        rxns = ko_rxn[sub]
        for rxn in rxns:
            if rxn in resultlist.keys():
                j = j + 1
                if j == 1:
                    temp = resultlist[rxn] # list in list
                    b = []
                    [b.extend(li) for li in temp]
                    strainOut = b
                else:
                    temp = resultlist[rxn] # list in list
                    b = []
                    [b.extend(li) for li in temp]
                    strainOut = list(set(strainOut).intersection(set(b)))
                ko_rxn_strain[sub] = strainOut
            else:
                print(sub + '_' + rxn)
    return ko_rxn_strain # dict for species for each substrate

# get index small function used in match function
def get_index3(lst=None, item=''):
    return [i for i in range(len(lst)) if lst[i] == item]

def find_sub_strain(ko_rxn):
    # link the substrate IDs with strains
    subusage = pd.read_csv('SubstrateUsage.tsv', sep='\t',header = 1)
    subusage = subusage.drop([0, 1, 2]) # describtive not actual value
    sub_strain = {}
    for sub in list(ko_rxn.keys()):
        if sub not in subusage.keys():
            sub = re.sub('.$', '', sub, 1) # for alternative pathways like lysine1 lysine2
        idx = subusage[(subusage[sub] == '1')| (subusage[sub] == 'v')].index.tolist()
        temp = subusage['name in the model'][idx]
        sub_strain[sub] = temp.values.tolist()
    return sub_strain

def match(ko_rxn_strain,sub_strain):
    accuracy = []
    for sub in list(ko_rxn_strain.keys()):
        sub_strain_exp = [s.lower() for s in sub_strain[sub]]
        sub_strain_blast = [re.sub('_16....$', '', str(i), 1) for i in ko_rxn_strain[sub]]
        sub_strain_blast = [s.lower() for s in sub_strain_blast]
        accuracy.append(len(list(set(sub_strain_exp).intersection(sub_strain_blast)))/len(sub_strain_exp))
    return accuracy

# plot result
def plot_filter(sub_strain,total,sub,rxnID):
    sub_strain_exp = [s.lower() for s in sub_strain[sub]]
    total = total.sort_values(by="identity")
    total.reset_index(drop=True, inplace=True)
    ret = [re.sub('_16....$', '', str(i), 1) for i in total["strain"]]
    total_strain = [s.lower() for s in ret]
    T = [get_index3(list(total_strain), i) for i in sub_strain_exp]
    b = []
    [b.extend(li) for li in T] # LIST IN LIST
    T_new = b
    #if len(T) == len(T_new):
    T_new.sort()
    allcov = T_new[0]
    # plot the figure
    fig = plt.figure(figsize=(5,4))

    ax1 = fig.add_subplot(111)
    x = np.arange(1, len(total["identity"])+1, 1)  # x = 1:1:total species number
    ax1.set_ylim([0, total["evalue"].max() * 1.1])
    ax1.plot(x, total["evalue"], 'b')
    ax1.set_ylabel('Y values for evalue', color='b',fontsize=20)
    ax1.set_title(sub + "_" + rxnID,fontsize=20)

    ax2 = ax1.twinx()  # this is the important function
    ax2.plot(x, total["identity"], 'r')
    ax2.set_xlim([0, len(x)])
    ax2.set_ylabel('Y values for identity', color='r',fontsize=20)
    ax2.set_xlabel('species',fontsize=20)

    ax2.bar(T_new, total["identity"][T_new], width = 0.35, color = '#000000', label = "bottomline")
    plt.show()
    #else:
       # print(list(set(sub_strain_exp).intersection(total_strain)))


if __name__ == "__main__":

    current_path = os.getcwd()
    os.chdir(current_path)
    datadir = current_path + "/ko"
    readfilenames = glob.glob(os.path.join(datadir, '*.fa'))# all ko
    resultdir = current_path + "/blast_result/" # all blast results folder

    # yeast ids, in order to speed it up, we choose only 20 sequences here, if yeast idS presented, then we prefer yeast sequences
    fungi = []
    with open('KEGG_org_yeast.txt', 'r') as f:
        for line in f:
            fungi.append(line.strip('\n'))  # all yeast ids

    # define ko_rxn_mapping
    f2 = csv.reader(open("rxnlistnew.txt","r")) # get mappings for ko sequence file with substrate phenotypes
    ko_rxn = {}
    for row in f2:
        k, v = row
        ko_rxn.setdefault(k, []).append(v)
    print(ko_rxn)
    # find phenotype related strains
    sub_strain = find_sub_strain(ko_rxn)

    resultlist = {}
    for conSeqfile in readfilenames:
        base = os.path.basename(conSeqfile)[:-3] # 'K03339_R05378'
        rxnID = base.split('_',1)[1] # R05378
        print("This KO is: %s" % base)
        resultID = glob.glob(resultdir + 'blastres_' + base + '*.txt') # all blast results for this ko
        outfile = open("../resultblast/best_hit_blast_" + base + ".csv", "w") # outfile filename
        kkk,total = filter_results(resultID,fungi,False,outfile) # filter the results above the cutoff, ref from fungi 35% ref from others 35%
        resultlist.setdefault(rxnID, []).append(kkk)
        i = 0
        for m in list(ko_rxn.values()):
            i = i + 1
            if rxnID in str(m):
                sub = list(ko_rxn.keys())[list(ko_rxn.values()).index(m)] # find the substrate matched to this reaction
                plot_filter(sub_strain, total, sub, rxnID) # plot the result to see the distribution of strians with this phenotype on blasted species
        print("-----------------------------------------")
    ko_rxn_strain = find_rxn_species(resultlist,ko_rxn)
    accurary = match(ko_rxn_strain,sub_strain) # calculate the accuray for the filtered results

    count_dict = dict()
    for i in resultlist.keys():
        a = list(resultlist[i])
        for item in a:
            if item in count_dict:
                count_dict[item] += 1
            else:
                count_dict[item] = 1