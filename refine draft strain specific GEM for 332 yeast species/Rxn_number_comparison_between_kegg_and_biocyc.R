# It can be found for the raven-biocyc method, if we use pi=55%, bitscore=110, the reaction number from
# each yeast species decreases a lot comparing with the original parameter pi=45% and bitscore=100. 
# Revised by Hongzhong 2019-11-03

# load library
library(readxl)
library(stringr)
library(tidyverse)
library(hongR)
source('function_general.R')

genome_yeasts <- read_excel("data/genome_summary_332_yeasts.xlsx") 
strain_index <- read_excel("data/332taxa_index.xlsx")
genome_yeasts$genomeID <- getSingleReactionFormula(strain_index$original_genome_id, strain_index$old_speceis_names,genome_yeasts$old_species_id)

#---------------------------------------------------
# input the model from biocyc and kegg
strain <- list.files('strain specific model from RAVEN_biocyc_55_110')
rxn_biocyc0 <- vector()
rxn_kegg0 <- vector()
all_rxn <- list()
for (i in strain) {
  print(i)
  inputfile1 <- paste('strain specific model from RAVEN_biocyc_55_110/',i, '/excelRxns.txt', sep = "")
  rxn_biocyc <- read_table2(inputfile1)
  inputfile2 <- paste('strain specific model from RAVEN_kegg/',i, '/excelRxns.txt', sep = "")
  rxn_kegg <- read_table2(inputfile2)
  rxn_biocyc0 <- c(rxn_biocyc0, length(unique(rxn_biocyc$`#`)))
  rxn_kegg0 <- c(rxn_kegg0, length(unique(rxn_kegg$`#`)))
  
}

reaction_num_summary <- data.frame(strain=strain, kegg=rxn_kegg0, biocyc=rxn_biocyc0, stringsAsFactors = FALSE)
plot(reaction_num_summary$kegg, reaction_num_summary$biocyc)








# when the original parameter is used
# input the model from biocyc and kegg
strain <- list.files('strain specific model from RAVEN_biocyc_45_100')
rxn_biocyc0 <- vector()
rxn_kegg0 <- vector()

for (i in strain) {
  print(i)
  inputfile1 <- paste('strain specific model from RAVEN_biocyc_45_100/',i, '/excelRxns.txt', sep = "")
  rxn_biocyc <- read_table2(inputfile1)
  inputfile2 <- paste('strain specific model from RAVEN_kegg/',i, '/excelRxns.txt', sep = "")
  rxn_kegg <- read_table2(inputfile2)
  rxn_biocyc0 <- c(rxn_biocyc0, length(unique(rxn_biocyc$`#`)))
  rxn_kegg0 <- c(rxn_kegg0, length(unique(rxn_kegg$`#`)))
  
}

reaction_num_summary <- data.frame(strain=strain, kegg=rxn_kegg0, biocyc=rxn_biocyc0, stringsAsFactors = FALSE)
plot(reaction_num_summary$kegg, reaction_num_summary$biocyc)


