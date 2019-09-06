# One general function to find the reaction ID from kegg database according to the KO ID
# Revised by Hongzhong 2019-8-12

#library(tidyverse)
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(readr)
library(hongR)

# This function is used to the get the detailed reaction information based on the KO ID
getRXNfromKO <- function(ss, outputRxn=TRUE) {
  
  # input: ss the KO annotation from kegg database, a dataframe should have the column 'ko'
  # KO_rxn mapping
  KO_rxn <- read.delim2("data/reaction-KO from kegg .txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(KO_rxn) <- c("rxnID", "ko")
  KO_rxn$ko <- str_replace_all(KO_rxn$ko, "ko:", "")
  # rxn formual in kegg
  rxn_kegg <- read.delim2("data/reaction_kegg summary.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  # KO_pathway mapping
  KO_pathway <- read.delim2("data/ko_pathway.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(KO_pathway) <- c("pathway", "ko")
  KO_pathway$ko <- str_replace_all(KO_pathway$ko, "ko:", "")
  KO_pathway <- filter(KO_pathway, str_detect(KO_pathway$pathway, "ko") == FALSE)
  # pathway annotation
  pathway <- read.delim2("data/pathway_list_kegg.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(pathway) <- c("pathway", "pathway_name")
  
  # obtain the rxn information based on ko mapping
  ss$rxns <- getMultipleReactionFormula(KO_rxn$rxnID, KO_rxn$ko, ss$ko)
  
  print('Find  rxn based on EC number mapping!!')
  # here we add another method: find the rxn based on the EC number while the EC number can be found from KO
  EC_ko_mapping <- read.delim2("data/EC_KO_mapping_kegg.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  EC_rxn_mapping <- read.delim2("data/EC_rxn_mapping_kegg.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  EC_ko_mapping$V2 <- str_replace_all(EC_ko_mapping$V2, "ko:", "")
  ss$EC <- getMultipleReactionFormula(EC_ko_mapping$V1, EC_ko_mapping$V2, ss$ko)
  # divided the ss into two types
  ss1 <- ss[!is.na(ss$rxns),]
  ss2 <- ss[is.na(ss$rxns),]
  # remove the row without ec number in ss2
  ss2 <- ss2[!is.na(ss2$EC),]
  # then we find the rxnID based on the ec number for each row
  # it seems that from ec number we can find more rections connect with each KO
  # YAR071W as an example
  for (i in 1:nrow(ss2)){
    ec0 <- ss2$EC[i]
    ec0 <- unlist(str_split(ec0, ";"))
    rxn_with_ec0 <- EC_rxn_mapping[EC_rxn_mapping$V1 == ec0,]
    if(nrow(rxn_with_ec0) >=1){
      rxn_choose <- rxn_with_ec0$V2
      rxn_choose <- paste0(rxn_choose,collapse = ";")
    } else(
      rxn_choose <- NA
    )
    print(rxn_choose)
    ss2$rxns[i] <- rxn_choose
  }
  # then further combine the ss1 and ss2
  ss_new <- rbind.data.frame(ss1, ss2)
  all_newRxn <- filter(ss_new, !is.na(rxns))
  
  
  # find other reaction information through KO and rxnID mapping
  all_kegg <- splitAndCombine(all_newRxn$rxns, all_newRxn$query, sep0 = ";")
  colnames(all_kegg) <- c("rxn", "query")
  
  all_kegg$formula <- getMultipleReactionFormula(rxn_kegg$reaction, rxn_kegg$rxnID, all_kegg$rxn)
  all_kegg$description <- getMultipleReactionFormula(rxn_kegg$name, rxn_kegg$rxnID, all_kegg$rxn)
  all_kegg$ko <- getMultipleReactionFormula(ss$ko, ss$query, all_kegg$query)
  all_kegg$pathway <- getMultipleReactionFormula(KO_pathway$pathway, KO_pathway$ko, all_kegg$ko)
  
  gene_pathway <- splitAndCombine(all_kegg$pathway, all_kegg$query, sep0 = ";")
  colnames(gene_pathway) <- c("pathID", "gene")
  gene_pathway$pathway_name <- getMultipleReactionFormula(pathway$pathway_name, pathway$pathway, gene_pathway$pathID)
  
  # change the pathway into sce
  gene_pathway$pathID <- str_replace_all(gene_pathway$pathID, "path:map", "sce")
  gene_pathway$pathway_name0 <- paste(gene_pathway$pathID, gene_pathway$pathway_name, sep = "  ")
  
  # merge the pathways
  all_kegg$pathway_name0 <- getMultipleReactionFormula(gene_pathway$pathway_name0, gene_pathway$gene, all_kegg$query)
  
  # obtain the possible GPRs
  unique_rxn <- unique(all_kegg$formula)
  index2 <- which(duplicated(all_kegg$formula) == FALSE)
  all_Rxn_kegg <- all_kegg[index2, ]
  all_Rxn_kegg$query <- getMultipleReactionFormula(all_kegg$query, all_kegg$formula, all_Rxn_kegg$formula)
  if(outputRxn){
    return(all_Rxn_kegg) #output the information based on the rxn, that is one rxn could be connected with several panIDs
  } else{
    return(all_kegg) #output the information based on the gene
  }
  
}


# an example
#colnames(sce_kegg) <- c('query','ko')
#sce_kegg <- filter(sce_kegg, ko != "")
# get the rxn based on KO
#sce_Rxn_kegg1 <- getRXNfromKO(ss = sce_kegg)




# this is previous version
# This function is used to the get the detailed reaction information based on the KO ID
getRXNfromKO_old_version <- function(ss, outputRxn=TRUE) {
  
  # input: ss the KO annotation from kegg database, a dataframe should have the column 'ko'
  # KO_rxn mapping
  KO_rxn <- read.delim2("data/reaction-KO from kegg .txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(KO_rxn) <- c("rxnID", "ko")
  KO_rxn$ko <- str_replace_all(KO_rxn$ko, "ko:", "")
  # rxn formual in kegg
  rxn_kegg <- read.delim2("data/reaction_kegg summary.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  # KO_pathway mapping
  KO_pathway <- read.delim2("data/ko_pathway.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(KO_pathway) <- c("pathway", "ko")
  KO_pathway$ko <- str_replace_all(KO_pathway$ko, "ko:", "")
  KO_pathway <- filter(KO_pathway, str_detect(KO_pathway$pathway, "ko") == FALSE)
  # pathway annotation
  pathway <- read.delim2("data/pathway_list_kegg.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(pathway) <- c("pathway", "pathway_name")
  
  # obtain the rxn information based on ko mapping
  ss$rxns <- getMultipleReactionFormula(KO_rxn$rxnID, KO_rxn$ko, ss$ko)
  panGenome_newRxn <- filter(ss, !is.na(rxns))
  
  panGenome_kegg <- splitAndCombine(panGenome_newRxn$rxns, panGenome_newRxn$query, sep0 = ";")
  colnames(panGenome_kegg) <- c("rxn", "query")
  
  panGenome_kegg$formula <- getMultipleReactionFormula(rxn_kegg$reaction, rxn_kegg$rxnID, panGenome_kegg$rxn)
  panGenome_kegg$description <- getMultipleReactionFormula(rxn_kegg$name, rxn_kegg$rxnID, panGenome_kegg$rxn)
  panGenome_kegg$ko <- getMultipleReactionFormula(ss$ko, ss$query, panGenome_kegg$query)
  panGenome_kegg$pathway <- getMultipleReactionFormula(KO_pathway$pathway, KO_pathway$ko, panGenome_kegg$ko)
  
  gene_pathway <- splitAndCombine(panGenome_kegg$pathway, panGenome_kegg$query, sep0 = ";")
  colnames(gene_pathway) <- c("pathID", "gene")
  gene_pathway$pathway_name <- getMultipleReactionFormula(pathway$pathway_name, pathway$pathway, gene_pathway$pathID)
  
  # change the pathway into sce
  gene_pathway$pathID <- str_replace_all(gene_pathway$pathID, "path:map", "sce")
  gene_pathway$pathway_name0 <- paste(gene_pathway$pathID, gene_pathway$pathway_name, sep = "  ")
  
  # merge the pathways
  panGenome_kegg$pathway_name0 <- getMultipleReactionFormula(gene_pathway$pathway_name0, gene_pathway$gene, panGenome_kegg$query)
  
  # obtain the possible GPRs
  unique_rxn <- unique(panGenome_kegg$formula)
  index2 <- which(duplicated(panGenome_kegg$formula) == FALSE)
  panGenome_Rxn_kegg <- panGenome_kegg[index2, ]
  panGenome_Rxn_kegg$query <- getMultipleReactionFormula(panGenome_kegg$query, panGenome_kegg$formula, panGenome_Rxn_kegg$formula)
  if(outputRxn){
    return(panGenome_Rxn_kegg) #output the information based on the rxn, that is one rxn could be connected with several panIDs
  } else{
    return(panGenome_kegg) #output the information based on the gene
  }
  
}

