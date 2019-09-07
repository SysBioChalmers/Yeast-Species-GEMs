# Compare the original s288c reactions from different sources: RAVEN_kegg and RAVEN_biocyc
# Explore the effect of parameters on the rxn number in RAVEN_biocyc
# Revised by Hongzhong 2019-8-10

# load library
library(readxl)
library(stringr)
library(tidyverse)
library(hongR)
source('function_general.R')

# input the rection information with RxnID from biocycID
s288c_Rxn_biocyc <- read_excel("data/biocyc_rxn_S288c.xlsx")
s288c_Rxn_biocyc <- s288c_Rxn_biocyc[!is.na(s288c_Rxn_biocyc$geneID),]
s288c_Rxn_biocyc$MNXID <- findRxnMNXid(rxnID = s288c_Rxn_biocyc$Reactions, id_type = 'metacyc')
s288c_Rxn_biocyc0 <- getRxnInfFromMNX(s288c_Rxn_biocyc, s288c_Rxn_biocyc$MNXID)


s288c_Rxn_kegg <- read.table("data/all s288c reactions from kegg.txt", header= TRUE, stringsAsFactors = FALSE)
s288c_Rxn_kegg$rxn <- str_replace_all(s288c_Rxn_kegg$rxn, "rn:", "")
s288c_Rxn_kegg$MNXID <- findRxnMNXid(rxnID = s288c_Rxn_kegg$rxn, id_type = 'kegg')
s288c_Rxn_kegg0 <- getRxnInfFromMNX(s288c_Rxn_kegg, s288c_Rxn_kegg$MNXID)
#test <- s288c_Rxn_kegg0[str_detect(s288c_Rxn_kegg0$query,'YMR217W'),]

# compare the common reaction from raven and from kegg and eggnog directly
# plot the vnn graph
library(VennDiagram)
s288c_biocyc <- unique(s288c_Rxn_biocyc$MNXID)
s288c_kegg <- unique(s288c_Rxn_kegg$MNXID)

#plot the graph
venn.diagram(x= list(s288c_biocyc = s288c_biocyc, s288c_kegg= s288c_kegg), 
             filename = "result/original_rxns_of_s288c_reactions in biocyc and kegg.png", height = 1000, width = 1000,resolution =300, imagetype="png", col="transparent",
             fill=c("blue", "red"),alpha = 0.50, cex=0.45, cat.cex=0.45)





# analysis the reactions which only exist in biocyc or kegg
biocyc_only <- setdiff(s288c_biocyc, s288c_kegg)
s288c_Rxn_biocyc_only <- s288c_Rxn_biocyc0[s288c_Rxn_biocyc0$MNXID %in%biocyc_only, ]


kegg_only <- setdiff(s288c_kegg, s288c_biocyc)
s288c_Rxn_kegg_only <- s288c_Rxn_kegg0[s288c_Rxn_kegg0$MNXID %in%kegg_only, ]



# systematic study how the pidenty and bitscore affect the result in RAVEN2-get model from biocyc for yeast
origi_m_gene_biocyc <- splitAndCombine(s288c_Rxn_biocyc$geneID, s288c_Rxn_biocyc$Reactions, sep0 = "//")
origi_m_gene_biocyc$v1 <- str_trim(origi_m_gene_biocyc$v1, side = "both")

id_mapping_biocyc <- read.table('data/s288c_geneID_mapping_biocyc.txt', header = TRUE, sep = "\t", stringsAsFactors = FALSE)
origi_m_gene_biocyc$gene <- getSingleReactionFormula(id_mapping_biocyc$Accession.1, id_mapping_biocyc$Gene.Name, origi_m_gene_biocyc$v1)
origi_m_gene <- unique(origi_m_gene_biocyc$gene) # metabolic gene in biocyc
origi_non_metabolic_gene <- setdiff(id_mapping_biocyc$Accession.1, origi_m_gene) # non metaboic gene in biocyc


# calculate the True positive, True negtive, False positive, False negative
pi <- vector()
bt <- vector()
score_all <- vector()
file_name <- list.files('data/biocyc_parameters_influence')

for (i in file_name) {
  print(i)
  n1 <- unlist(str_split(i, "_"))
  inputfile <- paste('data/biocyc_parameters_influence/',i, '/excelGenes.txt', sep = "")
  sce_biocyc <- read.table(inputfile, header =FALSE, sep = "\t", stringsAsFactors = FALSE)
  gene <- unique(sce_biocyc$V2)
  #TP
  tp <- intersect(gene, origi_m_gene)
  fn <- intersect(gene, origi_non_metabolic_gene)
  fp <- intersect(setdiff(id_mapping_biocyc$Accession.1, gene), origi_m_gene)
  tn <- intersect(setdiff(id_mapping_biocyc$Accession.1, gene), origi_non_metabolic_gene)
  score <- (length(tp)+length(tn))/(length(tp) + length(tn) + length(fp) + length(fn))
  #save the result
  pi <- c(pi, n1[3])
  bt <- c(bt, n1[4])
  score_all <- c(score_all, score)
}

summary_paramter <- data.frame(pident = pi, score=score_all, bitscore=bt, stringsAsFactors = FALSE)
summary_paramter$pident <- as.numeric(summary_paramter$pident)
summary_paramter$score <- as.numeric(summary_paramter$score)
summary_paramter$bitscore <- factor(summary_paramter$bitscore)
summary_paramter$bitscore <- factor(summary_paramter$bitscore, levels = c("70", "80", "90", "100", "110", "120", "130", "140", "150"))

# plot
ggplot(summary_paramter, aes(x = pident, y = score, shape = bitscore, colour = bitscore)) + geom_line(size = 1) +
  geom_point(size = 4) +
  xlab("pidenty") +
  ylab("Score") +
  theme_bw() +
  theme(legend.position = c(0.85, 0.6)) +
  theme(
    axis.text = element_text(size = 20, family = "Arial"),
    axis.title = element_text(size = 24, family = "Arial"),
    legend.text = element_text(size = 13, family = "Arial")
  ) +
  ggtitle("") +
  theme(panel.background = element_rect(fill = "white", color = "black", size = 1)) +
  ggsave(out <- paste("result/", "Effect of parameters on model reconstruction in RAVEN2-biocyc", ".eps", sep = ""), width = 8, height = 6, dpi = 300)



#--------------------------------
# systematic study how the pidenty and bitscore affect panYeast in RAVEN2-get model from biocyc for yeast
#--------------------------------

pi <- vector()
bt <- vector()
gene_all <- vector()
rxn_all <- vector()
file_name <- list.files('data/biocyc_pan_Yeast_different_parameter')

for (i in file_name) {
  print(i)
  #i <- "biocyc_panYeast2_45_100"
  n1 <- unlist(str_split(i, "_"))
  inputfile <- paste('data/biocyc_pan_Yeast_different_parameter/',i, '/excelGenes.txt', sep = "")
  inputfile2 <- paste('data/biocyc_pan_Yeast_different_parameter/',i, '/excelRxns.txt', sep = "")
  
  gene_biocyc <- read.table(inputfile, header =FALSE, sep = "\t", stringsAsFactors = FALSE)
  gene <- length(unique(gene_biocyc$V2))
  
  rxn_biocyc <- read_table2(inputfile2)
  rxn <- length(unique(rxn_biocyc$`#`))
  
  #save the result
  pi <- c(pi, n1[3])
  bt <- c(bt, n1[4])
  gene_all <- c(gene_all, gene)
  
  rxn_all <- c(rxn_all, rxn)
}

summary_paramter <- data.frame(pident = pi, gene=gene_all, rxn= rxn_all, bitscore=bt, stringsAsFactors = FALSE)
summary_paramter$pident <- as.numeric(summary_paramter$pident)
summary_paramter$gene <- as.numeric(summary_paramter$gene)
summary_paramter$rxn <- as.numeric(summary_paramter$rxn)
summary_paramter$bitscore <- factor(summary_paramter$bitscore)
summary_paramter$bitscore <- factor(summary_paramter$bitscore, levels = c("70", "80", "90", "100", "110", "120", "130", "140", "150"))

# plot
ggplot(summary_paramter, aes(x = pident, y = gene, shape = bitscore, colour = bitscore)) + geom_line(size = 1) +
  geom_point(size = 4) +
  xlab("pidenty") +
  ylab("Gene number") +
  theme_bw() +
  theme(legend.position = c(0.85, 0.6)) +
  theme(
    axis.text = element_text(size = 20, family = "Arial"),
    axis.title = element_text(size = 24, family = "Arial"),
    legend.text = element_text(size = 13, family = "Arial")
  ) +
  ggtitle("") +
  theme(panel.background = element_rect(fill = "white", color = "black", size = 1)) +
  ggsave(out <- paste("result/", "Effect of parameters on pan_model model reconstruction in RAVEN2-biocyc-gene number", ".eps", sep = ""), width = 8, height = 6, dpi = 300)


# plot
ggplot(summary_paramter, aes(x = pident, y = rxn, shape = bitscore, colour = bitscore)) + geom_line(size = 1) +
  geom_point(size = 4) +
  xlab("pidenty") +
  ylab("Reaction number") +
  theme_bw() +
  theme(legend.position = c(0.85, 0.6)) +
  theme(
    axis.text = element_text(size = 20, family = "Arial"),
    axis.title = element_text(size = 24, family = "Arial"),
    legend.text = element_text(size = 13, family = "Arial")
  ) +
  ggtitle("") +
  theme(panel.background = element_rect(fill = "white", color = "black", size = 1)) +
  ggsave(out <- paste("result/", "Effect of parameters on pan_model model reconstruction in RAVEN2-biocyc-reaction number", ".eps", sep = ""), width = 8, height = 6, dpi = 300)