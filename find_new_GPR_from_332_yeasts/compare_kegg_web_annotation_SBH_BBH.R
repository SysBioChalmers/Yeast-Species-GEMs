# this code is used to compare the different annotation-SBH and BBH using kegg web service.

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
# panGenome annotation based on : https://www.genome.jp/tools/kaas/  
# single-directional best hit to get the annotation
panGenome_kegg1 <- read_excel("data/fasta1_kegg_SBH.xlsx", col_names = FALSE)
colnames(panGenome_kegg1) <- c('query','ko')
panGenome_kegg2 <- read_excel("data/fasta2_kegg_SBH.xlsx", col_names = FALSE)
colnames(panGenome_kegg2) <- c('query','ko')
panGenome_kegg3 <- read_excel("data/fasta3_kegg_SBH.xlsx", col_names = FALSE)
colnames(panGenome_kegg3) <- c('query','ko')
panGenome_kegg_annotation <- rbind.data.frame(panGenome_kegg1, panGenome_kegg2,panGenome_kegg3)
panGenome_kegg_annotation <- filter(panGenome_kegg_annotation, ko != "")

sce_kegg <- read.delim2("data/sce_ko_2019_8.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(sce_kegg) <- c('query','ko')
sce_kegg <- filter(sce_kegg, ko != "")

# get the rxn based on KO
panGenome_Rxn_kegg1 <- getRXNfromKO(panGenome_kegg_annotation)
# panGenome_kegg1 <- getRXNfromKO(panGenome_kegg_annotation, outputRxn = FALSE)
sce_Rxn_kegg1 <- getRXNfromKO(sce_kegg)

# find the new Rxn ID from panGenome
newRxnID <- setdiff(panGenome_Rxn_kegg1$rxn, sce_Rxn_kegg1$rxn)
index1 <- which(panGenome_Rxn_kegg1$rxn %in% newRxnID ==TRUE)
newRxn_panGenome_kegg <- panGenome_Rxn_kegg1[index1,]








#bi-directional annotation by kegg web service
#It seems that this is not a good way as it get much less KO IDs for each protein.
#There are total 7610 KO IDs much less than the above analysis.
#Thus only the SBH result is used in the further analysis
panGenome_BBH <- read_excel("data/fasta_non_reference_kegg_BBH.xlsx", col_names = FALSE)
colnames(panGenome_BBH) <- c('query','ko')
panGenome_BBH <- filter(panGenome_BBH, ko != "")







