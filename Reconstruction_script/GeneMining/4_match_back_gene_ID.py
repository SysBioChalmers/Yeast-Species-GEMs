import json
import os
import glob
import numpy as np
import pandas as pd
import csv
from Bio import SeqIO
import re
# Import Biopython tools for running local BLASTX
from Bio.Blast.Applications import NcbiblastxCommandline
import matplotlib.pyplot as plt
plt.rc('font',family='Helvetica')

def create_blast_command(strainfile) :
    strain_id = os.path.basename(strainfile)[:-8]
    file1 = os.path.join("../sortedseqbyspecise/","%s.fa" % strain_id)
    if os.path.getsize(file1):
        file2 = strainfile
        fwd_out = os.path.join("../sortedseqbyspecise/", "%s_fwd.tab" % strain_id)

        # Create BLAST command-lines for forward and reverse BLAST searches
        # blastp -out ../Data/reciprocal_blast/Candida_albicans_fwd.tab -outfmt "6 qseqid sseqid pident qcovs qlen slen length bitscore evalue"
        # -query ../Data/reciprocal_blast/Candida_albicans_query.fasta -max_target_seqs 1 -subject ../Data/reciprocal_blast/Candida_albicans_ref.fasta
        fwd_blastx = NcbiblastxCommandline(query=file1, subject=file2, out=fwd_out,
                                           outfmt="6 qseqid sseqid pident qcovs qlen slen length bitscore evalue",
                                           max_target_seqs=3)

        # Inspect command-lines
        print("FORWARD: %s" % fwd_blastx)
        os.system(str(fwd_blastx))
    else:
        print("the query is empty: " + strain_id)

# https://widdowquinn.github.io/2018-03-06-ibioic/02-sequence_databases/05-blast_for_rbh.html#reciprocal
def blast_results(strain) :
    try:
        fwd_out = os.path.join("../sortedseqbyspecise/", "%s_fwd.tab" % strain)

        fwd_results = pd.read_csv(fwd_out, sep="\t", header=None)

        # Add headers to forward and reverse results dataframes
        headers = ["query", "subject", "identity", "coverage",
                   "qlength", "slength", "alength",
                   "bitscore", "E-value"]

        fwd_results.columns = headers

        # Create a new column in both dataframes: normalised bitscore
        fwd_results['norm_bitscore'] = fwd_results.bitscore/fwd_results.qlength

        # Create query and subject coverage columns in both dataframes
        fwd_results['qcov'] = fwd_results.alength/fwd_results.qlength
        fwd_results['scov'] = fwd_results.alength/fwd_results.slength

        # Clip maximum coverage values at 1.0
        fwd_results['qcov'] = fwd_results['qcov'].clip(upper=1)
        fwd_results['scov'] = fwd_results['scov'].clip(upper=1)

        fwd_results2 = fwd_results.loc[fwd_results["identity"]>95.0]
        gene_protein = fwd_results2.set_index('subject')['query'].to_dict()
        myset = set(fwd_results["query"])
        print(" blast best hits with Pidentity more than 95%%: %s" % len(gene_protein) + " total number: %s " % len(myset))
        fwd_results.to_excel("../mappedresult/" + strain + ".xlsx")
    except:
        gene_protein = {}
    return gene_protein

def sortresultbyspecies(strainlst,current_path) :
    readfilenames = glob.glob(os.path.join(current_path, '../best_hit_blast_*.csv'))
    #readfilenames = ['/Users/feiranl/Documents/Python/blast/resultblast/new_best_hit_blast_larab_R01903.csv', '/Users/feiranl/Documents/Python/blast/resultblast/new_best_hit_blast_K00849_R01092.csv', '/Users/feiranl/Documents/Python/blast/resultblast/new_best_hit_blast_K00965_R00955.csv', '/Users/feiranl/Documents/Python/blast/resultblast/new_best_hit_blast_K01784_R00291.csv', '/Users/feiranl/Documents/Python/blast/resultblast/new_best_hit_blast_lad_MNXR117337.csv', '/Users/feiranl/Documents/Python/blast/resultblast/new_best_hit_blast_K17066_R00608.csv', '/Users/feiranl/Documents/Python/blast/resultblast/new_best_hit_blast_K03781_R00009.csv', '/Users/feiranl/Documents/Python/blast/resultblast/new_best_hit_blast_dha1_R01440.csv', '/Users/feiranl/Documents/Python/blast/resultblast/new_best_hit_blast_K00121_R06983.csv', '/Users/feiranl/Documents/Python/blast/resultblast/new_best_hit_blast_K01070_R00527.csv', '/Users/feiranl/Documents/Python/blast/resultblast/new_best_hit_blast_K00122_R00519.csv', '/Users/feiranl/Documents/Python/blast/resultblast/new_best_hit_blast_K18338_R03772.csv', '/Users/feiranl/Documents/Python/blast/resultblast/new_best_hit_blast_K12661_R03774.csv', '/Users/feiranl/Documents/Python/blast/resultblast/new_best_hit_blast_inu_MNXR115392.csv', '/Users/feiranl/Documents/Python/blast/resultblast/new_best_hit_blast_YALI0F01650p_MNXR121148.csv', '/Users/feiranl/Documents/Python/blast/resultblast/new_best_hit_blast_YALI0F01606g_R11461.csv', '/Users/feiranl/Documents/Python/blast/resultblast/new_best_hit_blast_K21911_R11333.csv', '/Users/feiranl/Documents/Python/blast/resultblast/new_best_hit_blast_K02575_MNXR101992.csv', '/Users/feiranl/Documents/Python/blast/resultblast/new_best_hit_blast_K10534_R00794.csv', '/Users/feiranl/Documents/Python/blast/resultblast/new_best_hit_blast_K17877_R00787.csv', '/Users/feiranl/Documents/Python/blast/resultblast/new_best_hit_blast_K02575_MNXR101912.csv', '/Users/feiranl/Documents/Python/blast/resultblast/new_best_hit_blast_K17877_R00787.csv', '/Users/feiranl/Documents/Python/blast/resultblast/new_best_hit_blast_LYC1_R00454.csv']
    for conSeqfile in readfilenames:
        # read census.csv into a dataframe : census_df
        census_df = pd.read_csv(conSeqfile, header = 1)
        # rename the columns of the census dataframe
        census_df.columns = ["query", "subject", "identity", "qcov",
                             "qlength", "sstart", "send",
                             "bitscore", "E-value", "strain","seq"]
        a = census_df[['strain', 'seq']]
        a = a.dropna()
        a["strain"] = [re.sub('_16....$', '', str(i), 1) for i in a["strain"]]
        dict_seq = a.set_index('strain').T.to_dict('list')
        key = list(dict_seq.keys())
        for strain in strainlst:
            key_idx = [x for x in key if x.lower().endswith(strain.lower())]
            ofile = open('sortedseqbyspecise/' + strain + ".fa", "a")
            for i in range(len(key_idx)):
                ofile.write(">" + key_idx[i] + "\n" + ''.join(dict_seq[key_idx[i]]) + "\n")
            ofile.close()


if __name__ == "__main__":

    # reference protein sequence for each species
    current_path = os.getcwd()
    datadir = current_path + "/db_protein"
    dbs = glob.glob(os.path.join(datadir, '*.pep'))
    os.chdir(current_path)
    strainlst = [os.path.basename(x)[:-8] for x in dbs] # taking away .max.pep
    sortresultbyspecies(strainlst, current_path) # sort all blasted result by species, each species has a fasta file

    k = 0
    os.chdir(datadir)
    for strainfile in dbs: # Run Candida glabrata
        k = k + 1
        strain_id = os.path.basename(strainfile)[:-8]
        print("This strain is: %s" % strain_id.replace("_", " "))
        create_blast_command(strainfile) # get the linux command line for BLAST
        # # create_blast_command()
        gene_protein = blast_results(strain_id)
        if k == 1:
            allresult = gene_protein
        else:
            allresult.update(gene_protein.items())
        print("-----------------------------------------")
    # save all result into a file
    df = pd.DataFrame(data=allresult, index=[0])
    df = (df.T)
    df.to_excel('../mappedresult/allresult.xlsx', header=False)
