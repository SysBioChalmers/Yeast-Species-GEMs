import os
from os import path
import glob
from Bio.Blast.Applications import NcbitblastnCommandline
import pandas as pd

if __name__ == "__main__" :

    # generate the db for each strain for further use of tblastn
    current_path = os.getcwd()
    datadir = current_path + "/db"
    dbs = glob.glob(os.path.join(datadir, '*.fas'))
    os.chdir(datadir)
    for db in dbs[:332]:
        db_id = os.path.basename(db)[:-4]
        print(db)     #printed here to make it is printing the file
        cmd = "makeblastdb -in {} -dbtype nucl -input_type fasta -title {} -parse_seqids -out {} -logfile".format(db, db_id, db_id)
        os.system(cmd)

    # find all ko to balst
    os.chdir(current_path)
    datadir = current_path + "/ko"
    readfilenames = glob.glob(os.path.join(datadir, '*.fa'))

    #A directory for the output data
    outdir = '../blast_result_new'
    datadir = current_path + "/db" # stores reference genome sequences for all species
    os.chdir(datadir)
    for conSeqfile in readfilenames: # all ko
        for db in dbs: # all genome sequences
            db_id = os.path.basename(db)[:-4]
            print(conSeqfile)     #printed here the seached ko
            out_blast = path.join(outdir, "blastres_" + os.path.basename(conSeqfile)[:-3] + "_" + os.path.basename(db)[:-4] + ".txt")
            tblastn = NcbitblastnCommandline(query=conSeqfile, db=db_id, out=out_blast, evalue = 1e-10,
                     outfmt="6 qseqid sseqid pident qcovs qlen slen length sstart send bitscore evalue",
                     max_target_seqs=5)
            os.system(str(tblastn))
