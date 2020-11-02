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
summary_paramter$gene_set <- getSingleReactionFormula(genome_yeasts$`No. genes`, genome_yeasts$old_species_id,summary_paramter$strain)
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
summary_paramter2$gene_set <- getSingleReactionFormula(genome_yeasts$`No. genes`,genome_yeasts$old_species_id,summary_paramter2$strain)
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





#---------------------------------------------------
# connect the EC number and domain number with rxn number, gene set
# EC number is annotated using deepEC
# domain number is annotated using pfam

# files to store the EC number and domain number
EC_file <- "/Users/luho/Documents/332_yeast_EC/"
Domain_file <- "/Users/luho/Documents/332_yeast_domain/output/"

# parse the ec and domain for each yeast species
EC_all <- vector()
strain <- list.files(EC_file)
for (i in strain) {
  print(i)
  #i <- "biocyc_panYeast2_45_100"
  inputfile <- paste(EC_file,i, '/DeepEC_Result.txt', sep = "")
  items <- read.table(inputfile, header =TRUE, sep = "\t", stringsAsFactors = FALSE)
  total_num <- length(unique(items$Predicted.EC.number))
  #save the result
  EC_all <- c(EC_all, total_num)
}

EC_summary <- data.frame(strain = strain, EC=EC_all, stringsAsFactors = FALSE)
# strain id mapping
EC_summary$strain_unify <- getSingleReactionFormula(strain_index$old_speceis_names, strain_index$original_genome_id, EC_summary$strain)



Domain_all <- vector()
strain <- list.files(Domain_file)
for (i in strain) {
  print(i)
  #i <- "biocyc_panYeast2_45_100"
  inputfile <- paste(Domain_file,i, sep = "")
  items <- read.table(inputfile, stringsAsFactors = FALSE)
  total_num <- length(unique(items$V6))
  #save the result
  Domain_all <- c(Domain_all, total_num)
}

Domain_summary <- data.frame(strain = strain, Domain=Domain_all, stringsAsFactors = FALSE)
# strain id mapping
Domain_summary$strain <- str_replace_all(Domain_summary$strain,".txt", "")
Domain_summary$strain_unify <- getSingleReactionFormula(strain_index$old_speceis_names, strain_index$original_genome_id, Domain_summary$strain)

# connec the EC number and domain number with rxn obtained by biocyc
summary_paramter$EC <- getSingleReactionFormula(EC_summary$EC,EC_summary$strain_unify,summary_paramter$strain)
summary_paramter$Domain <- getSingleReactionFormula(Domain_summary$Domain,Domain_summary$strain_unify,summary_paramter$strain)
summary_paramter$EC <- as.numeric(summary_paramter$EC)
summary_paramter$Domain <- as.numeric(summary_paramter$Domain)

plot(summary_paramter$EC, summary_paramter$rxn)
plot(summary_paramter$Domain, summary_paramter$rxn)
plot(summary_paramter$gene, summary_paramter$rxn)

plotXYdotGraph2(data_frame=summary_paramter, paraX = 'EC', paraY = 'rxn', xlab_name='Unique EC number', ylab_name='rxn number_biocyc')
plotXYdotGraph2(data_frame=summary_paramter, paraX = 'Domain', paraY = 'rxn', xlab_name='Unique domain number', ylab_name='rxn number_biocyc')
plotXYdotGraph2(data_frame=summary_paramter, paraX = 'gene', paraY = 'rxn', xlab_name='Metabolic gene number_biocyc', ylab_name='rxn number_biocyc')

cor.test(summary_paramter$EC, summary_paramter$rxn)
cor.test(summary_paramter$Domain, summary_paramter$rxn)
cor.test(summary_paramter$gene, summary_paramter$rxn)


plotXYdotGraph2(data_frame=summary_paramter, paraX = 'EC', paraY = 'gene', xlab_name='Unique EC number', ylab_name='Metabolic gene number_biocyc')
plotXYdotGraph2(data_frame=summary_paramter, paraX = 'Domain', paraY = 'gene', xlab_name='Unique domain number', ylab_name='Metabolic gene number_biocyc')
plotXYdotGraph2(data_frame=summary_paramter, paraX = 'Domain', paraY = 'EC', xlab_name='Unique domain number', ylab_name='Unique EC number')


cor.test(summary_paramter$EC, summary_paramter$gene)
cor.test(summary_paramter$Domain, summary_paramter$gene)
cor.test(summary_paramter$Domain, summary_paramter$EC)


#---------------------------------------------------
# check the effect of expanded gene number on rxn number, gene set, EC number is annotated using deepEC
# input the gene expanded data
# it seems there are good correlation between the expanded gene number and gene size, rxn size.
gene_expansion <- read.table("data/gene_family_expansion_extraction_for_332_species.txt", header =TRUE, sep = "\t", stringsAsFactors = FALSE)
summary_paramter$expanded_gene <- getSingleReactionFormula(gene_expansion$expanded,gene_expansion$species,summary_paramter$strain)
summary_paramter$contracted_gene <- getSingleReactionFormula(gene_expansion$extracted,gene_expansion$species,summary_paramter$strain)

# remove the 11 out-group fungal species
summary_paramter <- summary_paramter[summary_paramter$expanded_gene!="NA",]
summary_paramter$expanded_gene <- as.numeric(summary_paramter$expanded_gene)
summary_paramter$contracted_gene <- as.numeric(summary_paramter$contracted_gene)
plotXYdotGraph2(data_frame=summary_paramter, paraX = 'expanded_gene', paraY = 'gene_set', xlab_name='Expanded gene number', ylab_name='gene_set')
cor.test(summary_paramter$gene_set, summary_paramter$expanded_gene)
plotXYdotGraph2(data_frame=summary_paramter, paraX = 'expanded_gene', paraY = 'rxn', xlab_name='Expanded gene number', ylab_name='rxn')
cor.test(summary_paramter$rxn, summary_paramter$expanded_gene)



#---------------------------------------------------
## Small task- to check whether the specific reactions existing in all yeast species
# RAVEN kegg
gene_all <- vector()
rxn_all <- vector()
exist_R01867 <- vector()
ec <- c('ec:1.3.98.1','ec:1.3.5.2')#,'ec:1.3.1.14')
ec_all <- vector()
strain <- list.files('strain specific model from RAVEN_kegg')
ec_rxn <- read.table('data/EC_rxn_mapping_kegg.txt', stringsAsFactors = FALSE) # ec rxn mapping from kegg
ec_rxn$V2 <- str_replace_all(ec_rxn$V2, "rn:", "")


for (i in strain) {
  print(i)
  #i <- "biocyc_panYeast2_45_100"
  inputfile <- paste('strain specific model from RAVEN_kegg/',i, '/excelGenes.txt', sep = "")
  inputfile2 <- paste('strain specific model from RAVEN_kegg/',i, '/excelRxns.txt', sep = "")
  
  gene_biocyc <- read.table(inputfile, header =FALSE, sep = "\t", stringsAsFactors = FALSE)
  gene <- length(unique(gene_biocyc$V2))
  
  rxn_biocyc <- read_table2(inputfile2)
  rxn <- length(unique(rxn_biocyc$`#`))
  
  rxn_detail <- data_frame(rxn=rxn_biocyc$`#`)
  rxn_detail <- merge(rxn_detail, ec_rxn, by.x = 'rxn', by.y = 'V2', all.x = TRUE)
  existence <- length(which(rxn_detail$rxn %in% 'R01867'))
  

  ec_existence <- length(which(rxn_detail$V1 %in% ec))
  #save the result
  exist_R01867 <- c(exist_R01867, existence)
  ec_all <- c(ec_all, ec_existence)
  gene_all <- c(gene_all, gene)
  rxn_all <- c(rxn_all, rxn)
}



