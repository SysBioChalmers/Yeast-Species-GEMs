# Compare the reaction from different sources: RAVEN, kegg and eggnog
# Revised by Hongzhong 2019-8-5

# load library
library(readxl)
library(stringr)
library(tidyverse)
library(hongR)
source('function_general.R')

#-----------------------------------------------------------
# initially compare the new RXN from different sources
#-----------------------------------------------------------

# RAVEN biocyc
# newRxn_biocyc <- read.table("data/newRxn_biocyc_RAVEN.txt", header= TRUE, stringsAsFactors = FALSE)
newRxn_biocyc <- read.table("data/newRxn_biocyc_RAVEN_55_110.txt", header= TRUE, stringsAsFactors = FALSE)
newRxn_biocyc$MNXID <- findRxnMNXid(rxnID = newRxn_biocyc$ID, id_type = 'metacyc')
newRxn_biocyc <- getRxnInfFromMNX(newRxn_biocyc,newRxn_biocyc$MNXID)

# RAVEN KEGG
newRxn_kegg <- read.table("data/newRxn_kegg_RAVEN.txt", header= TRUE, stringsAsFactors = FALSE)
newRxn_kegg$MNXID <- findRxnMNXid(rxnID = newRxn_kegg$ID, id_type = 'kegg')
newRxn_kegg <- getRxnInfFromMNX(newRxn_kegg,newRxn_kegg$MNXID)

# KEGG and eggnog web services
newRxn_kegg_eggnog <- read.table("data/newRxn_all based on kegg and eggnog annotation.txt", header= TRUE, stringsAsFactors = FALSE)
newRxn_kegg_eggnog$rxnID <- str_replace_all(newRxn_kegg_eggnog$rxnID, "rn:", "")
newRxn_kegg_eggnog$MNXID <- findRxnMNXid(rxnID = newRxn_kegg_eggnog$rxnID, id_type = 'kegg')
newRxn_kegg_eggnog <- getRxnInfFromMNX(newRxn_kegg_eggnog,newRxn_kegg_eggnog$MNXID)

rxn_kegg_web <- newRxn_kegg_eggnog[str_detect(newRxn_kegg_eggnog$type, 'kegg'),]
rxn_eggnog_web <- newRxn_kegg_eggnog[str_detect(newRxn_kegg_eggnog$type, 'eggnog'),]


# compare the common reaction from raven and from kegg and eggnog directly
# plot the vnn graph
kegg_web <- unique(rxn_kegg_web$MNXID)
eggnog_web <- unique(rxn_eggnog_web$MNXID)
RAVEN_kegg <- unique(newRxn_kegg$MNXID)
RAVEN_biocyc <- unique(newRxn_biocyc$MNXID)

new_rxn_all <- unique(c(kegg_web, eggnog_web,RAVEN_kegg,RAVEN_biocyc))


#plot the graph
VennDiagram::venn.diagram(x= list(kegg_web = kegg_web, eggnog_web = eggnog_web, RAVEN_kegg = RAVEN_kegg, RAVEN_biocyc = RAVEN_biocyc), 
             filename = "result/new reactions for 332 yeast species from different sources.png", height = 1000, width = 1000,resolution =300, imagetype="png", col="transparent",
             fill=c("blue","green","red", "grey"),alpha = 0.50, cex=0.45, cat.cex=0.45)






#---------------------------------------------
# if we only choose the balanced reactions
#---------------------------------------------
newRxn_biocyc_b <- newRxn_biocyc[newRxn_biocyc$balance_MNX=='true', ]
newRxn_kegg_b <- newRxn_kegg[newRxn_kegg$balance_MNX=='true', ]
newRxn_kegg_eggnog_b <- newRxn_kegg_eggnog[newRxn_kegg_eggnog$balance_MNX=='true', ]

rxn_kegg_web <- newRxn_kegg_eggnog_b[str_detect(newRxn_kegg_eggnog_b$type, 'kegg'),]
rxn_eggnog_web <- newRxn_kegg_eggnog_b[str_detect(newRxn_kegg_eggnog_b$type, 'eggnog'),]

kegg_web <- unique(rxn_kegg_web$MNXID)
eggnog_web <- unique(rxn_eggnog_web$MNXID)
RAVEN_kegg <- unique(newRxn_kegg_b$MNXID)
RAVEN_biocyc <- unique(newRxn_biocyc_b$MNXID)


#plot the graph
venn.diagram(x= list(kegg_web = kegg_web, eggnog_web = eggnog_web, RAVEN_kegg = RAVEN_kegg, RAVEN_biocyc = RAVEN_biocyc), 
             filename = "result/new balanced reactions for 332 yeast species from different sources.png", height = 1000, width = 1000,resolution =300, imagetype="png", col="transparent",
             fill=c("blue","green","red", "grey"),alpha = 0.50, cex=0.45, cat.cex=0.45)







#----------------------------------------------------------------------------------------------
# specially, here we found MNXID for the panID through the EC number based on eggnog annotation
#----------------------------------------------------------------------------------------------
newEC_eggnog <- read.table("data/new EC based eggnog annotation.txt", header= TRUE, stringsAsFactors = FALSE)
newEC_eggnog$MNXID <- findRxnMNXidFromEC(newEC_eggnog$EC)

newEC_eggnog$rxn_num <- str_count(newEC_eggnog$MNXID, ";")
# as a EC number could be connected with so many reactions, we choose ec with no more than 5 rxns
newEC_eggnog_filter <- newEC_eggnog[newEC_eggnog$rxn_num <= 5,]
newEC_eggnog0 <- splitAndCombine(newEC_eggnog_filter$MNXID, newEC_eggnog_filter$query, sep0 = ";")
rxn_ec1 <- unique(newEC_eggnog0$v1)


# Also we found MNXID for the panID through the EC number based on deepec
newEC_deepec <- read.table("data/newEC_predicted_by_deep_ec_for_pan_genome.txt", header= TRUE, stringsAsFactors = FALSE)
newEC_deepec$MNXID <- findRxnMNXidFromEC(newEC_deepec$Predicted.EC.number)
newEC_deepec$rxn_num <- str_count(newEC_deepec$MNXID, ";")
# as a EC number could be connected with so many reactions, we choose ec with no more than 5 rxns
newEC_deepec_filter <- newEC_deepec[newEC_deepec$rxn_num <= 5,]
newEC_deepec0 <- splitAndCombine(newEC_deepec_filter$MNXID, newEC_deepec_filter$Query.ID, sep0 = ";")
rxn_ec2 <- unique(newEC_deepec0$v1)
# combine new EC from different sources
rxn_ec_all <- union(rxn_ec1, rxn_ec2)
rxn_ec_all0 <- data.frame(MNXID =rxn_ec_all, stringsAsFactors = FALSE)
rxn_ec_all0 <- getRxnInfFromMNX(rxn_ec_all0, rxn_ec_all0$MNXID)
rxn_ec_all0_b <- rxn_ec_all0[rxn_ec_all0$balance_MNX=='true', ]
rxn_ec_combine <- unique(rxn_ec_all0_b$MNXID)

# plot the graph
venn.diagram(x= list(kegg_web = kegg_web,rxn_ec = rxn_ec_combine, RAVEN_kegg = RAVEN_kegg, RAVEN_biocyc = RAVEN_biocyc), 
             filename = "result/new balanced reactions for 332 yeast species with rxn found by EC number.png", height = 1000, width = 1000,resolution =300, imagetype="png", col="transparent",
             fill=c("blue","green","red", "grey"),alpha = 0.50, cex=0.45, cat.cex=0.45)
