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
    file1 = os.path.join("sortedseqbyspecise/","%s.fa" % strain_id)
    file2 = strainfile
    fwd_out = os.path.join("sortedseqbyspecise/", "%s_fwd.tab" % strain_id)
    rev_out = os.path.join("sortedseqbyspecise/", "%s_rev.tab" % strain_id)

    # Create BLAST command-lines for forward and reverse BLAST searches
    # blastp -out ../Data/reciprocal_blast/Candida_albicans_fwd.tab -outfmt "6 qseqid sseqid pident qcovs qlen slen length bitscore evalue"
    # -query ../Data/reciprocal_blast/Candida_albicans_query.fasta -max_target_seqs 1 -subject ../Data/reciprocal_blast/Candida_albicans_ref.fasta
    fwd_blastx = NcbiblastxCommandline(query=file1, subject=file2, out=fwd_out,
                                       outfmt="6 qseqid sseqid pident qcovs qlen slen length bitscore evalue",
                                       max_target_seqs=3)

    rev_blastx = NcbiblastxCommandline(query=file2, subject=file1, out=rev_out,
                                       outfmt="6 qseqid sseqid pident qcovs qlen slen length bitscore evalue",
                                       max_target_seqs=3)
    # Inspect command-lines
    print("FORWARD: %s" % fwd_blastx)
    #print("REVERSE: %s" % rev_blastx)
    os.system(str(fwd_blastx))
    #os.system(str(rev_blastx))

    return fwd_blastx, #rev_blastx

# https://widdowquinn.github.io/2018-03-06-ibioic/02-sequence_databases/05-blast_for_rbh.html#reciprocal
def blast_results(strain) :
    fwd_out = os.path.join("sortedseqbyspecise/", "%s_fwd.tab" % strain)
    #rev_out = os.path.join("sortedseqbyspecise/", "%s_rev.tab" % strain)
    # fwd_out = os.path.join("../Data/reciprocal_blast/", "Candida_albicans_fwd.tab")
    # rev_out = os.path.join("../Data/reciprocal_blast/", "Candida_albicans_rev.tab")

    fwd_results = pd.read_csv(fwd_out, sep="\t", header=None)
    #rev_results = pd.read_csv(rev_out, sep="\t", header=None)

    # Add headers to forward and reverse results dataframes
    headers = ["query", "subject", "identity", "coverage",
               "qlength", "slength", "alength",
               "bitscore", "E-value"]

    fwd_results.columns = headers
    rev_results.columns = headers

    # Create a new column in both dataframes: normalised bitscore
    fwd_results['norm_bitscore'] = fwd_results.bitscore/fwd_results.qlength
    rev_results['norm_bitscore'] = rev_results.bitscore/rev_results.qlength

    # Create query and subject coverage columns in both dataframes
    fwd_results['qcov'] = fwd_results.alength/fwd_results.qlength
    rev_results['qcov'] = rev_results.alength/rev_results.qlength
    fwd_results['scov'] = fwd_results.alength/fwd_results.slength
    rev_results['scov'] = rev_results.alength/rev_results.slength

    # Clip maximum coverage values at 1.0
    # fwd_results['qcov'] = fwd_results['qcov'].clip_upper(1)
    # rev_results['qcov'] = rev_results['qcov'].clip_upper(1)
    # fwd_results['scov'] = fwd_results['scov'].clip_upper(1)
    # rev_results['scov'] = rev_results['scov'].clip_upper(1)
    fwd_results['qcov'] = fwd_results['qcov'].clip(upper=1)
    rev_results['qcov'] = rev_results['qcov'].clip(upper=1)
    fwd_results['scov'] = fwd_results['scov'].clip(upper=1)
    rev_results['scov'] = rev_results['scov'].clip(upper=1)

    return fwd_results, rev_results

def plot_rbbh(fwd_results, rev_results, strain) : # reciprocal blast best hits
    # Merge forward and reverse results
    rbbh = pd.merge(fwd_results, rev_results[['query', 'subject']],
                    left_on='subject', right_on='query',
                    how='outer')

    # Discard rows that are not RBH
    rbbh = rbbh.loc[rbbh.query_x == rbbh.subject_y]

    # Group duplicate RBH rows, taking the maximum value in each column
    rbbh = rbbh.groupby(['query_x', 'subject_x']).max()

    # Plot 2D density histograms

    # Calculate 2D density histograms for counts of matches at several coverage levels
    # (H, xedges, yedges) = np.histogram2d(rbbh.qcov, rbbh.scov, bins=20)

    # # Create a 1x2 figure array
    # fig, ax = plt.subplots(1, 1, figsize=(6, 6), sharex=True, sharey=True)

    # # Plot histogram for RBBH
    # im = ax.imshow(H, cmap=plt.cm.Blues, norm=LogNorm(),
    #                  extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
    #                  origin='lower', aspect=1)
    # ax.set_title("RBBH")
    # ax.set_xlabel("query")
    # ax.set_ylabel("subject")

    # # Add colourbar
    # fig.colorbar(im, ax=ax)

    plt.savefig("../figure/rbbh/%s_rbbh.png" % strain, dpi=400)

    return rbbh

def id_mapping(rbbh, strain) :
    rbbh2 = rbbh.loc[rbbh["identity"]>95.0]
    # for col in rbbh2.columns :
    #     print(col)
    gene_protein = rbbh2.set_index('subject_y')['query_y'].to_dict()
    print("Reciprocal blast best hits with Pidentity more than 95%%: %s" % len(gene_protein))

    with open("/json/%s.json" % strain) as f :
        proteinAll = json.load(f)

    print("The number of protein complex data: %s" % len(proteinAll))

    print("Reciprocal blast best hits: %s" % len(rbbh.index))  # including all the reciprocal blast best hits results

    gene_mapping = list()
    gene_unmapping = list()
    for gene in proteinAll :
        try :
            gene["id"] = gene_protein[list(gene.keys())[0]]
            gene_mapping.append(gene)
        except :
            gene_unmapping.append(list(gene.keys())[0])

    print("The number of gene mapping: %s" % len(gene_mapping))  # including all the reciprocal blast best hits results pidentity > 95%
    # print(gene_unmapping)


    with open("../processed_data/%s_include_id.json" % strain, "w") as outfile :
        json.dump(gene_mapping, outfile, indent=4)


    rbbh.to_excel("/mappedresult/rbbh" + strain + ".xlsx")
    # fwd_results.to_csv("../Data/reciprocal_blast/fwd_Candida_albicans.csv")
    # rev_results.to_csv("../Data/reciprocal_blast/rev_Candida_albicans.csv")
    # len(rbbh.index)
def sortresultbyspecies(strainlst,filtedresults) :
    for conSeqfile in filtedresults:
        print(conSeqfile)
        # read census.csv into a dataframe : census_df
        census_df = pd.read_csv(conSeqfile, header = 1)
        # rename the columns of the census dataframe
        census_df.columns = ["query", "subject", "identity", "qcov",
                             "qlength", "sstart", "send",
                             "bitscore", "E-value", "strain","seq"]
        a = census_df[['strain', 'seq']]
        dict_seq = a.set_index('strain').T.to_dict('list')
        key = list(dict_seq.keys())
        for strain in strainlst:
            key_idx = [x for x in key if strain.lower() in x.lower()]
            if len(key_idx) > 0:
                ofile = open('../sortedseqbyspecise/' + strain + ".fa", "a")
                for i in range(len(key_idx)):
                    ofile.write(">" + key_idx[i] + "\n" + ''.join(dict_seq[key_idx[i]]) + "\n")
                ofile.close()

if __name__ == "__main__":

    # reference protein sequence for each species
    current_path = os.getcwd()
    datadir = current_path + "/db_protein" # protein sequence folder
    dbs = glob.glob(os.path.join(datadir, '*.pep'))  # get all species protein sequence file names
    os.chdir(datadir)
    strainlst = [os.path.basename(x)[:-8] for x in dbs] # taking away .max.pep
    filtedresults = glob.glob(os.path.join(current_path, '../resultblast/best_hit_blast_*.csv'))
    sortresultbyspecies(strainlst, filtedresults) # sort all blasted result by species, each species has a fasta file to fasten the blast

    for strainfile in dbs[1:2] : # Run Candida glabrata
        strain_id = os.path.basename(strainfile)[:-8]
        print("This strain is: %s" % strain_id.replace("_", " "))
        create_blast_command(strainfile) # get the linux command line for BLAST
        # # create_blast_command()
        fwd_results, rev_results = blast_results(strain_id)
        # plot_oneway(fwd_results, rev_results, strain)
        rbbh = plot_rbbh(fwd_results, rev_results, strain_id)
        id_mapping(rbbh, strain_id)
        print("-----------------------------------------")