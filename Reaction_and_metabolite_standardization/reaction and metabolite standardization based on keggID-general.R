# Give the detailed annotation for reaction based on its KEGG ID.
# Also give the name of metabolite from the latest yeast GEM
# Revised by Hongzhong 2019-8-5


# load library
library(readxl)
library(stringr)
library(tidyverse)
library(hongR)
source('function_general.R')

#input the rection information with RxnID from kegg databas
#input data should have a column named 'rxnID', in this column, there are reactionID from kegg database,like R06513
#one example
rxn <- read.delim2("data/newRxn_all based on kegg and eggnog annotation.txt", stringsAsFactors = FALSE)
rxn$rxnID <- str_replace_all(rxn$rxnID, "rn:", "")

#---------------------------------------------------------------
# Find other reaction ID, like Rhea ID, MNX ID based on KEGG ID
#---------------------------------------------------------------

#find the rhea id based on keggid
#input the rhea id and keggid
rhea2kegg <- read_tsv('data/rhea2kegg_reaction.tsv')
rhea2kegg$RHEA_ID <- paste("RHEA:",rhea2kegg$RHEA_ID, sep = "")
rhea2kegg$MASTER_ID <- paste("RHEA:",rhea2kegg$MASTER_ID, sep = "")
rxn$RHEA_ID <- getSingleReactionFormula(rhea2kegg$MASTER_ID,rhea2kegg$ID,rxn$rxnID)


#find the rhea id based on biocycid
rhea2biocyc <- read_tsv('data/rhea2metacyc.tsv')
rhea2biocyc$RHEA_ID <- paste("RHEA:",rhea2biocyc$RHEA_ID, sep = "")
rhea2biocyc$MASTER_ID <- paste("RHEA:",rhea2biocyc$MASTER_ID, sep = "")
rxn$RHEA_ID2 <- getSingleReactionFormula(rhea2biocyc$MASTER_ID,rhea2biocyc$ID,rxn$rxnID)

#merge the two rhea ID
for (i in seq(length(rxn$rxnID))){
  if(rxn$RHEA_ID[i]=='NA'){
    rxn$RHEA_ID[i] <- rxn$RHEA_ID2[i]
  } else{
    rxn$RHEA_ID[i] <- rxn$RHEA_ID[i]
  }
  
}


#find the metnetid based on keggid
reac_metnet <- read_tsv('data/reac_prop_metanetx.tsv')
#reac_xref <- read_tsv('data/reac_xref.tsv')
reac_xref <- read_tsv('data/reac_xref_2019_8.tsv') # using the latest version from MNX database
reac_metnet$Source <- str_replace_all(reac_metnet$Source, "kegg:","")
reac_xref$XREF <- str_replace_all(reac_xref$XREF, "kegg:","")
rxn$metnetID <- getMultipleReactionFormula(reac_xref$MNX_ID,reac_xref$XREF,rxn$rxnID)

#find the metnetid based on rheaid
rxn$RHEA_ID <-str_replace_all(rxn$RHEA_ID, "RHEA","rhea")
rxn$metnetID2 <- getMultipleReactionFormula(reac_metnet$MNX_ID,reac_metnet$Source,rxn$RHEA_ID) # why not use reac_xref

#merge the metnetID
#metnetID 2 is the final version of ID for each reaction
for (i in seq(length(rxn$metnetID2))){
  if(is.na(rxn$metnetID2[i])){
    rxn$metnetID2[i] <- rxn$metnetID[i]
  } else{
    rxn$metnetID2[i] <- rxn$metnetID2[i]
  }

  if(str_detect(rxn$rxnID[i],'MNXR' )){
    rxn$metnetID2[i] <- rxn$rxnID[i]
  } else{
    rxn$metnetID2[i] <- rxn$metnetID2[i]
  }

}


#-----------------------------------------------------------------------------------------------
# now some rxnID from kegg could be mapped onto Rhea or metnetx database, while some can't
# in the above first case, the metabolite formula could be obtained based on Rhea or metanetx
# while for the second case, keggid can be obtained for each reactions
#-----------------------------------------------------------------------------------------------

#rxn_final <- filter(rxn, !is.na(rxn$metnetID2) | rxn$RHEA_ID != 'NA' )
rxn_final <- rxn #sometimes need some filters

# find reaction based on metnet id and rhea id
# input the metnet database
reac_metnet <- read_tsv('data/reac_prop_metanetx.tsv')
rxn_final$formula_metnet <- getSingleReactionFormula(reac_metnet$Description,reac_metnet$MNX_ID,rxn_final$metnetID2)
rxn_final$formula_metnet_balance <- getSingleReactionFormula(reac_metnet$Balance,reac_metnet$MNX_ID,rxn_final$metnetID2)
rxn_final$Equation <- getMultipleReactionFormula(reac_metnet$Equation,reac_metnet$MNX_ID,rxn_final$metnetID2)
rxn_final$EC <- getMultipleReactionFormula(reac_metnet$EC,reac_metnet$MNX_ID,rxn_final$metnetID2)

# input the rhea database
rhea_reaction_summary <- read_csv("data/rhea reaction summary.csv")
rhea_reaction_summary$masterID <- paste("rhea:",rhea_reaction_summary$masterID, sep = "")
rxn_final$formula_rhea <- getMultipleReactionFormula(rhea_reaction_summary$formula,rhea_reaction_summary$masterID,rxn_final$RHEA_ID)

for (i in seq(length(rxn_final$rxnID))){
  if(rxn_final$formula_metnet[i] =='NA'){
    rxn_final$formula_metnet[i] <- rxn_final$formula_rhea[i]
    #rxn_final$formula_metnet_balance[i] <- 'true'
    } else{
    rxn_final$formula_metnet[i] <- rxn_final$formula_metnet[i] 
    rxn_final$formula_metnet_balance[i] <- rxn_final$formula_metnet_balance[i]
    
    }
  
}


# increase the reveserbility check based on Rhea database
# for example based on rhea database
reverse <- str_count(rxn_final$formula_rhea, ";")
for (i in seq(length(rxn_final$rxnID))){
  if(!is.na(rxn_final$formula_rhea[i]) & reverse[i] >= 1 ){
    rxn_final$reverse_rhea[i] <- 'yes'
  } 
  else if (!is.na(rxn_final$formula_rhea[i]) & reverse[i] != 1) {
    rxn_final$reverse_rhea[i] <- 'no'
    
  } else{
    rxn_final$reverse_rhea[i] <- 'not sure'
  }
  
}


# check the reveserbility based on metacyc database
# get the metacyc id based on the MNXID
# it should be noted that a MNXID could have multiple metacyc ids
mnx_biocyc <- reac_xref[str_detect(reac_xref$XREF,'metacyc'),]
mnx_biocyc$XREF <- str_replace_all(mnx_biocyc$XREF, "metacyc:", "")
biocyc <- read.delim2("data/All_instances_of_Reactions_in_MetaCyc.txt", stringsAsFactors = FALSE)
biocyc_reaction <- read_excel("data/metcyc_rxn_RAVEN.xlsx")
rxn_final$biocycID <- getSingleReactionFormula(mnx_biocyc$XREF, mnx_biocyc$MNX_ID, rxn_final$metnetID2)
# find the reaction reversiblity from biocyc database
rxn_final$reverse_biocyc <- getSingleReactionFormula(biocyc$Reaction.Direction, biocyc$Reactions, rxn_final$biocycID)
rxn_final$formula_biocyc <- getSingleReactionFormula(biocyc_reaction$formula, biocyc_reaction$Reactions, rxn_final$biocycID)


# check the reveserbility based on modelSeed database
# get the modelSEEDid based on the MNXID
mnx_seed <- reac_xref[str_detect(reac_xref$XREF,'seed'),]
mnx_seed$XREF <- str_replace_all(mnx_seed$XREF, "seed:", "")
seed <- read_excel("data/ModelSEED-reactions-db.xls")
rxn_final$seedID <- getSingleReactionFormula(mnx_seed$XREF, mnx_seed$MNX_ID, rxn_final$metnetID2)
# find the reaction reversiblity from seed database
rxn_final$reverse_seed <- getSingleReactionFormula(seed$`THERMODYNAMIC FEASIBILTY`, seed$DATABASE, rxn_final$seedID)
rxn_final$formula_seed <- getSingleReactionFormula(seed$`NAME EQ`, seed$DATABASE, rxn_final$seedID)






#-------------------------------------------------------------------------------------------
#find the standard formula of metabolite based on reaction compostions from metnetx database
#-------------------------------------------------------------------------------------------
rxn_final0 <- select(rxn_final, rxnID, Equation)
colnames(rxn_final0) <- c('ID0', 'Equation')

#remove the compartment information from metnet database
rxn_met <- splitRxnToMetabolite(rxn_final0, sep0 = " = ")

chem_prop <- read_tsv('data/chem_prop_metanetx.tsv')
rxn_met$Formula <- getMultipleReactionFormula(chem_prop$Formula,chem_prop$MNX_ID,rxn_met$MetID)
rxn_met$Charge <- getMultipleReactionFormula(chem_prop$Charge,chem_prop$MNX_ID,rxn_met$MetID)
rxn_met$Description <- getMultipleReactionFormula(chem_prop$Description,chem_prop$MNX_ID,rxn_met$MetID)

#find the keggid and chebiID based on metnetID
chem_xref <- read_tsv('data/chem_xref.tsv')
chem_xref_kegg <- chem_xref[str_detect(chem_xref$XREF,'kegg:'),]
chem_xref_chebi <- chem_xref[str_detect(chem_xref$XREF,'chebi:'),]
rxn_met$chebi_multiple <- getMultipleReactionFormula(chem_xref_chebi$XREF,chem_xref_chebi$MNX_ID,rxn_met$MetID)
rxn_met$kegg_multiple <- getMultipleReactionFormula(chem_xref_kegg$XREF,chem_xref_kegg$MNX_ID,rxn_met$MetID)

chem_prop_chebi <- chem_prop[str_detect(chem_prop$Source, "chebi"), ]
chem_prop_kegg <- chem_prop[str_detect(chem_prop$Source, "kegg"), ]

rxn_met$chebi <- getMultipleReactionFormula(chem_prop_chebi$Source,chem_prop_chebi$MNX_ID,rxn_met$MetID)
rxn_met$kegg <- getMultipleReactionFormula(chem_prop_kegg$Source,chem_prop_kegg$MNX_ID,rxn_met$MetID)
# find the keggid based on chebiID
chebi <- read_excel("data/chebi_compound_comprehensive.xlsx",  sheet = "Sheet1")
chebi$chebiID <-str_replace_all(chebi$chebiID, "CHEBI","chebi")
rxn_met$kegg2 <- getMultipleReactionFormula(chebi$keggID,chebi$chebiID,rxn_met$chebi)
rxn_met$kegg2[is.na(rxn_met$kegg2)] <- "NA"

for (i in seq_along(rxn_met$kegg2)){
  if(rxn_met$kegg2[i]=="NA"){
    rxn_met$kegg2[i] <- rxn_met$kegg_multiple[i]
  } else{
    rxn_met$kegg2[i] <-  rxn_met$kegg2[i]
  }
}
rxn_met$kegg2 <- str_replace_all(rxn_met$kegg2, "kegg:", "")
rxn_met <- subset(rxn_met, select = -c(kegg, kegg_multiple, chebi_multiple))



#--------------------------------------------------
# find the standard metabolites name from yeast GEM
#--------------------------------------------------
model_8_13 <- read_excel("data/model_8.13.xlsx", sheet = "MNX", col_names = FALSE)
colnames(model_8_13) <- c('MNXID','MetName')
model_8_13 <- model_8_13 %>% separate(., MetName, into = c('met','compartment'), sep =" \\[")
model_8_13 <- model_8_13[!duplicated(model_8_13$MNXID),]
rxn_met$name_GEM <- getSingleReactionFormula(model_8_13$met,model_8_13$MNXID,rxn_met$MetID)
write.table(rxn_met, "result/met for new reactions based on panGenome.txt", row.names = FALSE, sep = "\t")
