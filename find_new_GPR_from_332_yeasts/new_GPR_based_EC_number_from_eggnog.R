# this code is used to find the new GPR from panID based on EC number from eggnog web services annotation
# Revised by hongzhong 2019-8-12

#library(tidyverse)
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(readr)
library(hongR)
source('main_function.R')

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

#Establish relation between ec and gene
panGenome_eggnog1 <- select(panGenome_eggnog, query_name, EC)
panGenome_eggnog1 <- filter(panGenome_eggnog1, EC!="")
panGenome_eggnog1_EC <- splitAndCombine(panGenome_eggnog1$EC,panGenome_eggnog1$query,sep0 = ",")
colnames(panGenome_eggnog1_EC) <- c('EC','query')

sce_eggnog1 <- select(sce_eggnog, query_name, EC)
sce_eggnog1 <- filter(sce_eggnog1, EC!="")
sce_eggnog1_EC <- splitAndCombine(sce_eggnog1$EC,sce_eggnog1$query,sep0 = ",")
colnames(sce_eggnog1_EC) <- c('EC','query')

#new ec
newEC <- setdiff(panGenome_eggnog1_EC$EC, sce_eggnog1_EC$EC)
newEC_eggnog <- panGenome_eggnog1_EC[panGenome_eggnog1_EC$EC %in% newEC, ]

#save the result
write.table(newEC_eggnog, "result/new EC based eggnog annotation.txt", row.names = FALSE, sep = "\t")













# The followed will be removed!!
#-----------------------------------------

# Also we can get the new kegg reactions directly from eggnog annotation
# Establish relation between kegg RXN and gene
# This method is not used
# Revised by Hongzhong Lu
panGenome_eggnog1 <- select(panGenome_eggnog, query_name, KEGG_Reaction)
panGenome_eggnog1 <- filter(panGenome_eggnog1, KEGG_Reaction!="")
panGenome_eggnog1_KEGG_Reaction <- splitAndCombine(panGenome_eggnog1$KEGG_Reaction,panGenome_eggnog1$query_name,sep0 = ",")
colnames(panGenome_eggnog1_KEGG_Reaction) <- c('KEGG_Reaction','query')

sce_eggnog1 <- select(sce_eggnog, query_name, KEGG_Reaction)
sce_eggnog1 <- filter(sce_eggnog1, KEGG_Reaction!="")
sce_eggnog1_KEGG_Reaction <- splitAndCombine(sce_eggnog1$KEGG_Reaction,sce_eggnog1$query_name,sep0 = ",")
colnames(sce_eggnog1_KEGG_Reaction) <- c('KEGG_Reaction','query')
#new kegg rxn
newKEGG_Reaction <- setdiff(panGenome_eggnog1_KEGG_Reaction$KEGG_Reaction, sce_eggnog1_EC$KEGG_Reaction)
newKEGG_Reaction_eggnog <- panGenome_eggnog1_KEGG_Reaction[panGenome_eggnog1_KEGG_Reaction$KEGG_Reaction %in% newKEGG_Reaction, ]
# further establish the map between rxn and genes
length(unique(newKEGG_Reaction_eggnog$query))
# Note :
# based on eggnog web service, there are 7436 panIDs connect with 1873 new reactions. So many!!
