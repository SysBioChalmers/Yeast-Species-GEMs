# this code is used to find the new GPR from panID based on kegg and eggnog web services
# if the new KO id occured in kegg database and eggnog, then the annotation from kegg will be used

#library(tidyverse)
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(readr)
library(hongR)
source('main_function.R')

# panGenome annotation based on : https://www.genome.jp/tools/kaas/  
# single-directional best hit to get the annotation
# It seems that two_directial best hit takes a very long time
panGenome_kegg1 <- read_excel("data/fasta1_kegg_SBH.xlsx", col_names = FALSE)
colnames(panGenome_kegg1) <- c('query','ko')
panGenome_kegg2 <- read_excel("data/fasta2_kegg_SBH.xlsx", col_names = FALSE)
colnames(panGenome_kegg2) <- c('query','ko')
panGenome_kegg3 <- read_excel("data/fasta3_kegg_SBH.xlsx", col_names = FALSE)
colnames(panGenome_kegg3) <- c('query','ko')
panGenome_kegg_annotation <- rbind.data.frame(panGenome_kegg1, panGenome_kegg2,panGenome_kegg3)
panGenome_kegg_annotation <- filter(panGenome_kegg_annotation, ko != "")
# input the s288c annotation
sce_kegg <- read.delim2("data/sce_ko_2019_8.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(sce_kegg) <- c('query','ko')
sce_kegg <- filter(sce_kegg, ko != "")


# get the rxn based on KO
panGenome_Rxn_kegg1 <- getRXNfromKO(panGenome_kegg_annotation)
sce_Rxn_kegg1 <- getRXNfromKO(sce_kegg)
test <- sce_Rxn_kegg1[str_detect(sce_Rxn_kegg1$query,'YMR217W'), ]

# find the new Rxn ID from panGenome
newRxnID <- setdiff(panGenome_Rxn_kegg1$rxn, sce_Rxn_kegg1$rxn)
index1 <- which(panGenome_Rxn_kegg1$rxn %in% newRxnID ==TRUE)
newRxn_panGenome_kegg <- panGenome_Rxn_kegg1[index1,]

rxn_newGene <- splitAndCombine(newRxn_panGenome_kegg$query, newRxn_panGenome_kegg$rxn, sep0 = ";")


# panGenome annotation based on : http://eggnogdb.embl.de/#/app/home
# it is quite detailed for the panGenome annotation from eggnog web service
# the results could the direct evidences for the article
# with the mapping onto BiGG, no new reactions were found for these new panID compared with mapping onto kegg database
# it may be due to the fact that the metabolic coverage in the bigg database is quite small
egg10 <- read.delim2("data/query_seqs.fa.emapper.annotations1", header = FALSE, sep = "\t",stringsAsFactors = FALSE)
egg20 <- read.delim2("data/query_seqs.fa.emapper.annotations2", header = FALSE, sep = "\t",stringsAsFactors = FALSE)
egg30 <- read.delim2("data/query_seqs.fa.emapper.annotations3", header = FALSE, sep = "\t",stringsAsFactors = FALSE)
header0 <- egg10[1,]
header0 <- unlist(header0)[1:17]
other_col <- c('annot lvl', 'matching OGs', 'Best OG', 'COG cat', 'description')
header_new <- c(header0, other_col)
egg10 <- egg10[-1,]
egg20 <- egg20[-1,]
egg30 <- egg30[-1,]
panGenome_eggnog <- rbind.data.frame(egg20, egg10, egg30)
colnames(panGenome_eggnog) <- header_new

# input the eggnog annotation for sce s288c
sce_eggnog <- read.delim2("data/sce288c.fa.emapper.annotations", header = FALSE, sep = "\t",stringsAsFactors = FALSE)
sce_eggnog <- sce_eggnog[-1,]
colnames(sce_eggnog) <- header_new
write.table(panGenome_eggnog, "result/panGenome_eggnog.txt", row.names = FALSE, sep = "\t")

#analysis the results mapping onto kegg database based on eggnog annotation
panGenome_eggnog1 <- select(panGenome_eggnog, query_name, KEGG_ko)
panGenome_eggnog1 <- filter(panGenome_eggnog1, KEGG_ko!="")
panGenome_eggnog1_rxn <- splitAndCombine(panGenome_eggnog1$KEGG_ko,panGenome_eggnog1$query,sep0 = ",")
colnames(panGenome_eggnog1_rxn) <- c('ko','query')
panGenome_eggnog1_rxn$ko <- str_replace_all(panGenome_eggnog1_rxn$ko, "ko:", "")

#Obtain the GPR from kegg database 
panGenome_Rxn_eggnog1<- getRXNfromKO(panGenome_eggnog1_rxn)
#Obtain rxn connected with each panID
#panGenome_eggnog1<- getRXNfromKO(panGenome_eggnog1_rxn, outputRxn = FALSE)
#get the new rxn

#analysis the results mapping onto kegg database based on eggnog annotation
sce_eggnog1 <- select(sce_eggnog, query_name, KEGG_ko)
sce_eggnog1 <- filter(sce_eggnog1, KEGG_ko!="")
sce_eggnog1_rxn <- splitAndCombine(sce_eggnog1$KEGG_ko,sce_eggnog1$query,sep0 = ",")
colnames(sce_eggnog1_rxn) <- c('ko','query')
sce_eggnog1_rxn$ko <- str_replace_all(sce_eggnog1_rxn$ko, "ko:", "")

#Obtain the GPR from kegg database 
sce_Rxn_eggnog1 <- getRXNfromKO(sce_eggnog1_rxn)
#Obtain rxn connected with each panID
#panGenome_eggnog1<- getRXNfromKO(panGenome_eggnog1_rxn, outputRxn = FALSE)
#get the new rxn
newRxnID <- setdiff(panGenome_Rxn_eggnog1$rxn, sce_Rxn_eggnog1$rxn)
index1 <- which(panGenome_Rxn_eggnog1$rxn %in% newRxnID ==TRUE)
newRxn_panGenome_eggnog <- panGenome_Rxn_eggnog1[index1,]





#comparsion between the new RXN obtained by kegg web service and eggnog gene annotation
#all reactions could be further divided into three types
kegg_only <- setdiff(newRxn_panGenome_kegg$rxn, newRxn_panGenome_eggnog$rxn)
eggnog_only <- setdiff(newRxn_panGenome_eggnog$rxn,newRxn_panGenome_kegg$rxn)
kegg_and_eggnog <- intersect(newRxn_panGenome_kegg$rxn, newRxn_panGenome_eggnog$rxn)
kegg_or_eggnog <- union(newRxn_panGenome_kegg$rxn, newRxn_panGenome_eggnog$rxn)

newRxn_panGenome_kegg$source <- "kegg_web"
newRxn_panGenome_eggnog$source <- "eggnog_web"

newRxn_all <-data.frame(rxnID=kegg_or_eggnog, stringsAsFactors = FALSE)
newRxn_all$type <- "NA"
kegg_only_index <- which(newRxn_all$rxn %in% kegg_only==TRUE)
eggnog_only_index <- which(newRxn_all$rxn %in% eggnog_only==TRUE)
both_index <-  which(newRxn_all$rxn %in% kegg_and_eggnog==TRUE)
newRxn_all$type[kegg_only_index] <- "kegg_only"
newRxn_all$type[eggnog_only_index] <- "eggnog_only"
newRxn_all$type[both_index] <- "eggnog_and_kegg"

newRxn_all$panID_kegg <- getMultipleReactionFormula(newRxn_panGenome_kegg$query,newRxn_panGenome_kegg$rxn,newRxn_all$rxnID)
newRxn_all$panID_eggnog <- getMultipleReactionFormula(newRxn_panGenome_eggnog$query,newRxn_panGenome_eggnog$rxn,newRxn_all$rxnID)

#save the result
write.table(sce_Rxn_kegg1, "result/all s288c reactions from kegg.txt", row.names = FALSE, sep = "\t")
write.table(newRxn_all, "result/newRxn_all based on kegg and eggnog annotation.txt", row.names = FALSE, sep = "\t")

