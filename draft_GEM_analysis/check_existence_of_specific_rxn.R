# Compare the original s288c reactions from different sources: RAVEN, kegg
# Revised by Hongzhong 2020-7-27

# load library
library(readxl)
library(stringr)
library(tidyverse)


#---------------------------------------------------
## Small task- to check whether the specific reactions existing in all yeast species
# RAVEN kegg
# RAVEN biocyc information should be also considered in future version

gene_all <- vector()
rxn_all <- vector()
exist_R01867 <- vector()
ec <- c('ec:1.3.98.1','ec:1.3.5.2')#,'ec:1.3.1.14')
ec_all <- vector()
strain <- list.files('../ComplementaryData/draft_yeast_GEMs/strain_specific_model_from_RAVEN_kegg')
ec_rxn <- read.table('data/EC_rxn_mapping_kegg.txt', stringsAsFactors = FALSE) # ec rxn mapping from kegg
ec_rxn$V2 <- str_replace_all(ec_rxn$V2, "rn:", "")


for (i in strain) {
  print(i)
  #i <- "biocyc_panYeast2_45_100"
  inputfile <- paste('../ComplementaryData/draft_yeast_GEMs/strain_specific_model_from_RAVEN_kegg/',i, '/excelGenes.txt', sep = "")
  inputfile2 <- paste('../ComplementaryData/draft_yeast_GEMs/strain_specific_model_from_RAVEN_kegg/',i, '/excelRxns.txt', sep = "")
  
  gene_one_species <- read.table(inputfile, header =FALSE, sep = "\t", stringsAsFactors = FALSE)
  gene <- length(unique(gene_one_species$V2))
  
  rxn_one_species <- read_table2(inputfile2)
  rxn <- length(unique(rxn_one_species$`#`))
  
  rxn_detail <- data_frame(rxn=rxn_one_species$`#`)
  rxn_detail <- merge(rxn_detail, ec_rxn, by.x = 'rxn', by.y = 'V2', all.x = TRUE)
  existence <- length(which(rxn_detail$rxn %in% 'R01867'))
  

  ec_existence <- length(which(rxn_detail$V1 %in% ec))
  #save the result
  exist_R01867 <- c(exist_R01867, existence)
  ec_all <- c(ec_all, ec_existence)
  gene_all <- c(gene_all, gene) # just calculate total gene number from draft model
  rxn_all <- c(rxn_all, rxn) # just calculate total rxn number from draft model
}

