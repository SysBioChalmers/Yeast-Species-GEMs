# Compare the original s288c reactions from different sources: RAVEN, kegg
# Revised by Hongzhong 2020-7-27

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
# RAVEN biocyc
gene_all <- vector()
rxn_all <- vector()
strain <- list.files('../ComplementaryData/draft_yeast_GEMs/strain_specific_model_from_RAVEN_biocyc_55_110')

for (i in strain) {
  print(i)
  inputfile <- paste('../ComplementaryData/draft_yeast_GEMs/strain_specific_model_from_RAVEN_biocyc_55_110/',i, '/excelGenes.txt', sep = "")
  inputfile2 <- paste('../ComplementaryData/draft_yeast_GEMs/strain_specific_model_from_RAVEN_biocyc_55_110/',i, '/excelRxns.txt', sep = "")
  gene_biocyc <- read.table(inputfile, header =FALSE, sep = "\t", stringsAsFactors = FALSE)

}



#---------------------------------------------------
# RAVEN kegg
gene_all <- vector()
rxn_all <- vector()
strain <- list.files('../ComplementaryData/draft_yeast_GEMs/strain_specific_model_from_RAVEN_kegg')

all_m_kegg <- vector()
for (i in strain) {
  print(i)
  inputfile <- paste('../ComplementaryData/draft_yeast_GEMs/strain_specific_model_from_RAVEN_kegg/',i, '/excelGenes.txt', sep = "")
  inputfile2 <- paste('../ComplementaryData/draft_yeast_GEMs/strain_specific_model_from_RAVEN_kegg/',i, '/excelRxns.txt', sep = "")
  gene_biocyc <- read.table(inputfile, header =FALSE, sep = "\t", stringsAsFactors = FALSE)
  gene_list_biocyc <- gene_biocyc$V2
  all_m_kegg <- c(all_m_kegg, gene_list_biocyc)
}

# how to change gene id into OG id










