# Compare the original s288c reactions from different sources: RAVEN, kegg
# Revised by Hongzhong 2019-8-10

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
strain <- list.files('strain specific model from RAVEN_biocyc_55_110')

for (i in strain) {
  print(i)
  #i <- "biocyc_panYeast2_45_100"
  inputfile <- paste('strain specific model from RAVEN_biocyc_55_110/',i, '/excelGenes.txt', sep = "")
  inputfile2 <- paste('strain specific model from RAVEN_biocyc_55_110/',i, '/excelRxns.txt', sep = "")
  
  gene_biocyc <- read.table(inputfile, header =FALSE, sep = "\t", stringsAsFactors = FALSE)
  gene <- length(unique(gene_biocyc$V2))
  
  rxn_biocyc <- read_table2(inputfile2)
  rxn <- length(unique(rxn_biocyc$`#`))
  
  #save the result

  gene_all <- c(gene_all, gene)
  
  rxn_all <- c(rxn_all, rxn)
}

summary_paramter <- data.frame(strain = strain, gene=gene_all, rxn= rxn_all, stringsAsFactors = FALSE)
summary_paramter$gene <- as.numeric(summary_paramter$gene)
summary_paramter$rxn <- as.numeric(summary_paramter$rxn)
summary_paramter$gene_set <- getSingleReactionFormula(genome_yeasts$`No. genes`, genome_yeasts$genomeID,summary_paramter$strain)
summary_paramter$gene_set <- as.numeric(summary_paramter$gene_set)


# plot
plot(genome_yeasts$`Assembly size`, genome_yeasts$`No. genes`, xlab = 'Genome size (bp)', ylab = "Gene set")
genome_yeasts0 <- data.frame(x=genome_yeasts$`Assembly size`, y=genome_yeasts$`No. genes`, stringsAsFactors = FALSE)
plotXYdotGraph(data_frame=genome_yeasts0, paraX = 'x', paraY = 'y', xlab_name='Genome size (bp)', ylab_name='Gene set')

plot(density(summary_paramter$gene))
plot(density(summary_paramter$rxn))
plot(summary_paramter$gene, summary_paramter$rxn)
plot(summary_paramter$gene_set, summary_paramter$rxn, xlab = 'Gene set', ylab = "Rxn number_biocyc")
plot(summary_paramter$gene_set, summary_paramter$gene, xlab = 'Gene set', ylab = "Gene number_biocyc")


# plot of high quality
plotDensityGraph(data_frame=summary_paramter, para = "gene", xlab_name='Metabolic gene number_biocyc')
plotDensityGraph(data_frame=summary_paramter, para = "rxn", xlab_name='Metabolic rxn number_biocyc')
plotXYdotGraph(data_frame=summary_paramter, paraX = 'gene', paraY = 'rxn', xlab_name='Metabolic gene number_biocyc', ylab_name='rxn number_biocyc')
plotXYdotGraph(data_frame=summary_paramter, paraX = 'gene_set', paraY = 'rxn', xlab_name='Gene set', ylab_name='rxn number_biocyc')
plotXYdotGraph(data_frame=summary_paramter, paraX = 'gene_set', paraY = 'gene', xlab_name='Gene set', ylab_name='Metabolic gene_biocyc')

plotXYdotGraph2(data_frame=summary_paramter, paraX = 'gene', paraY = 'rxn', xlab_name='Metabolic gene number_biocyc', ylab_name='rxn number_biocyc')
plotXYdotGraph2(data_frame=summary_paramter, paraX = 'gene_set', paraY = 'rxn', xlab_name='Gene set', ylab_name='rxn number_biocyc')
plotXYdotGraph2(data_frame=summary_paramter, paraX = 'gene_set', paraY = 'gene', xlab_name='Gene set', ylab_name='Metabolic gene_biocyc')







#---------------------------------------------------
# RAVEN kegg
gene_all <- vector()
rxn_all <- vector()
strain <- list.files('strain specific model from RAVEN_kegg')

for (i in strain) {
  print(i)
  #i <- "biocyc_panYeast2_45_100"
  inputfile <- paste('strain specific model from RAVEN_kegg/',i, '/excelGenes.txt', sep = "")
  inputfile2 <- paste('strain specific model from RAVEN_kegg/',i, '/excelRxns.txt', sep = "")
  
  gene_biocyc <- read.table(inputfile, header =FALSE, sep = "\t", stringsAsFactors = FALSE)
  gene <- length(unique(gene_biocyc$V2))
  
  rxn_biocyc <- read_table2(inputfile2)
  rxn <- length(unique(rxn_biocyc$`#`))
  
  #save the result
  
  gene_all <- c(gene_all, gene)
  
  rxn_all <- c(rxn_all, rxn)
}

summary_paramter2 <- data.frame(strain = strain, gene=gene_all, rxn= rxn_all, stringsAsFactors = FALSE)
summary_paramter2$gene <- as.numeric(summary_paramter2$gene)
summary_paramter2$rxn <- as.numeric(summary_paramter2$rxn)
summary_paramter2$gene_set <- getSingleReactionFormula(genome_yeasts$`No. genes`,genome_yeasts$genomeID,summary_paramter2$strain)
summary_paramter2$gene_set <- as.numeric(summary_paramter2$gene_set)

# plot
plot(density(summary_paramter2$gene))
plot(density(summary_paramter2$rxn))
plot(summary_paramter2$gene, summary_paramter2$rxn)
plot(summary_paramter2$gene_set, summary_paramter2$rxn, xlab = 'Gene set', ylab = "Rxn number_kegg")
plot(summary_paramter2$gene_set, summary_paramter2$gene, xlab = 'Gene set', ylab = "Gene number_kegg")

# plot of high quality
plotDensityGraph(data_frame=summary_paramter2, para = "gene", xlab_name='Metabolic gene number_kegg')
plotDensityGraph(data_frame=summary_paramter2, para = "rxn", xlab_name='Metabolic rxn number_kegg')
plotXYdotGraph(data_frame=summary_paramter2, paraX = 'gene', paraY = 'rxn', xlab_name='Metabolic gene number_kegg', ylab_name='rxn number_kegg')
plotXYdotGraph(data_frame=summary_paramter2, paraX = 'gene_set', paraY = 'rxn', xlab_name='Gene set', ylab_name='rxn number_kegg')
plotXYdotGraph(data_frame=summary_paramter2, paraX = 'gene_set', paraY = 'gene', xlab_name='Gene set', ylab_name='Metabolic gene_kegg')

plotXYdotGraph2(data_frame=summary_paramter2, paraX = 'gene', paraY = 'rxn', xlab_name='Metabolic gene number_kegg', ylab_name='rxn number_kegg')
plotXYdotGraph2(data_frame=summary_paramter2, paraX = 'gene_set', paraY = 'rxn', xlab_name='Gene set', ylab_name='rxn number_kegg')
plotXYdotGraph2(data_frame=summary_paramter2, paraX = 'gene_set', paraY = 'gene', xlab_name='Gene set', ylab_name='Metabolic gene_kegg')

