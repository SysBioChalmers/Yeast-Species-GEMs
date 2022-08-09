# Merge the model from Biocyc and KEGG
# It can be found for the raven-biocyc method, if we use pi=55%, bitscore=110, the reaction number from
# each yeast species decreases a lot comparing with the original parameter pi=45% and bitscore=100. 
# Revised by Hongzhong 2020-07-01

# load library
library(readxl)
library(stringr)
library(tidyverse)
library(hongR)
source('function_general.R')

# input the strain information
genome_yeasts <- read_excel("data/genome_summary_332_yeasts.xlsx") 
strain_index <- read_excel("data/332taxa_index.xlsx")
genome_yeasts$genomeID <- getSingleReactionFormula(strain_index$original_genome_id, strain_index$old_speceis_names,genome_yeasts$old_species_id)
# input the model from biocyc and kegg
strain <- list.files('../ComplementaryData/draft_yeast_GEMs/strain_specific_model_from_RAVEN_biocyc_55_110')
rxn_biocyc0 <- vector()
rxn_kegg0 <- vector()

for (n in strain) {
  print(n)
  inputfile1 <- paste('../ComplementaryData/draft_yeast_GEMs/strain_specific_model_from_RAVEN_biocyc_55_110/',n, '/draft_GEM.tsv', sep = "")
  rxn_biocyc <- read_tsv(inputfile1)
  inputfile2 <- paste('../ComplementaryData/draft_yeast_GEMs/strain_specific_model_from_RAVEN_kegg/',n, '/draft_GEM.tsv', sep = "")
  rxn_kegg <- read_tsv(inputfile2)
  rxn_biocyc0 <- rxn_biocyc[, c('ID','EQUATION','GENE ASSOCIATION')]
  rxn_kegg0 <- rxn_kegg[, c('ID','EQUATION','GENE ASSOCIATION')]
  colnames(rxn_biocyc0) <- c('rxnID','equation','geneID')
  colnames(rxn_kegg0) <- c('rxnID','equation','geneID')
  rxn_biocyc0$MNXID <- findRxnMNXid(rxnID = rxn_biocyc0$rxnID, id_type = "metacyc")
  rxn_kegg0$MNXID <- findRxnMNXid(rxnID = rxn_kegg0$rxnID, id_type = "kegg")
  
  
  # get the unique rxnID
  RAVEN_kegg <- unique(rxn_kegg0$MNXID)
  RAVEN_biocyc <- unique(rxn_biocyc0$MNXID)
  # merge the rxn from different source
  union_rxn <- unique(c(RAVEN_kegg, RAVEN_biocyc))
  union_rxn_ann <- data.frame(rxnID = union_rxn, stringsAsFactors = FALSE)
  union_rxn_ann$source <- NA
  union_rxn_ann$panID_RAVEN_kegg <- NA
  union_rxn_ann$panID_RAVEN_biocyc <- NA
  union_rxn_ann$panID_inter <- NA # the interaction panID sets from at least two different sources
  union_rxn_ann$panID_union <- NA # the union panID sets from different sources
  
  # summarize the panID
  # integrate the source
  E3 <- union_rxn_ann$rxnID %in% RAVEN_kegg
  E4 <- union_rxn_ann$rxnID %in% RAVEN_biocyc
  source_merge <- vector()
  for (i in 1:length(E3)) {
    if (E3[i]) {
      s3 <- "RAVEN_kegg"
    } else {
      s3 <- NA
    }
    if (E4[i]) {
      s4 <- "RAVEN_biocyc"
    } else {
      s4 <- NA
    }
    s <- paste(s3, s4, sep = ";")
    source_merge <- c(source_merge, s)
  }
  
  source_merge <- str_replace_all(source_merge, "NA;", "")
  source_merge <- str_replace_all(source_merge, ";NA", "")
  union_rxn_ann$source <- source_merge
  
  # union_rxn_ann0 <- filter(union_rxn_ann, balance_MNX=='true') # 1275 balanced reactions
  # union_rxn_ann0$source_num <- str_count(union_rxn_ann0$source, ";") # 605 reactions from at least two evidences
  union_rxn_ann$panID_RAVEN_kegg <- getMultipleReactionFormula(rxn_kegg0$geneID, rxn_kegg0$MNXID, union_rxn_ann$rxnID)
  union_rxn_ann$panID_RAVEN_biocyc <- getMultipleReactionFormula(rxn_biocyc0$geneID, rxn_biocyc0$MNXID, union_rxn_ann$rxnID)
  
  # unify the panID for each reactions
  g_union_all <- vector()
  g_inter_all <- vector()
  for (i in 1:nrow(union_rxn_ann)) {
    print(i)
    g3 <- union_rxn_ann$panID_RAVEN_kegg[i]
    g3 <- unlist(str_split(g3, " or "))
    g3 <- str_trim(g3, side = "both")
    
    g4 <- union_rxn_ann$panID_RAVEN_biocyc[i]
    g4 <- unlist(str_split(g4, " or "))
    g4 <- str_trim(g4, side = "both")
    g_union_string <- paste(union(g3,g4), collapse = ";") %>%
      str_replace_all('NA;',"") %>%
      str_replace_all(';NA', "")
    
    g_inter_string <- paste(intersect(g3,g4), collapse = ";")

    g_union_all <- c(g_union_all, g_union_string)
    g_inter_all <- c(g_inter_all, g_inter_string)
  }
  
  union_rxn_ann$panID_inter <- g_inter_all
  union_rxn_ann$panID_union <- g_union_all
  write.table(union_rxn_ann, paste("result/Merged_strain_specific_model_from_RAVEN2/",n,'.txt', sep=""), row.names = FALSE, sep = "\t")
}






# Conduce the pan and core rxn analysis for 332 yeast species, not including the 11 outgroup
strain_332 <- list.files("result/Merged_strain_specific_model_from_RAVEN2/")
strain_332 <- str_replace_all(strain_332, ".txt", "")
all_rxn <- vector()
for (n in strain_332) {
  print(n)
  inputfile1 <- paste("result/Merged_strain_specific_model_from_RAVEN2/",n,'.txt', sep="")
  model <- read.table(inputfile1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  all_rxn <- c(all_rxn, model$rxnID)
}

all_rxn_unique <- unique(all_rxn)

rxn_existence <- data_frame(rxnID= all_rxn_unique)

for (n in strain_332) {
  print(n)
  inputfile1 <- paste("result/Merged_strain_specific_model_from_RAVEN2/",n,'.txt', sep="")
  model <- read.table(inputfile1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  rxn <- model$rxnID
  existence <- all_rxn_unique %in% rxn
  rxn_existence[, n] <- existence
}

rxn_existence0 <- select(rxn_existence, -c(1))
sum0 <-apply(rxn_existence0,1,sum)
rxn_result <- data_frame(rxnID=all_rxn_unique, Occur_num=sum0)


# plot the result
ggplot(rxn_result, aes(x=Occur_num)) + geom_histogram(binwidth=5) +
  ylab('Reaction number') + 
  theme(panel.background = element_rect(fill = "white", color="black", size = 1),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=24, family="Arial"),
        legend.text = element_text(size=20, family="Arial"))

ggplot(rxn_result, aes(x=Occur_num)) + stat_ecdf() +
  ylab('Percentage of reactions') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  theme(axis.text=element_text(size=10, family="Arial"),
        axis.title=element_text(size=24, family="Arial"),
        legend.text = element_text(size=20, family="Arial"))

write.table(rxn_result, "result/pan_reactome.txt", row.names = FALSE, sep = "\t")


# calculate the pan and core reactome
strain_num <- 1:332
core_all <- vector()
pan_all <- vector()
accessory_all <- vector()
for (i in strain_num){
  print(i)
  strain_choose <- sample(strain_332, i)
  print(strain_choose)
  rxn_set <- select(rxn_existence0, strain_choose)
  sum_set <-apply(rxn_set,1,sum)
  rxn_set_new <- data_frame(rxnID=rxn_existence$rxnID, num= sum_set)
  core <-rxn_set_new[rxn_set_new$num==i,]
  pan <- rxn_set_new[rxn_set_new$num >=1,]
  accessory <- setdiff(pan$rxnID, core$rxnID)
  core_all <- c(core_all, length(core$rxnID))
  pan_all <- c(pan_all, length(pan$rxnID))
  accessory_all <- c(accessory_all, length(accessory))
}
reactome <- data_frame(num= strain_num, pan=pan_all, core=core_all, accessory=accessory_all)


# repeat the above steps ten times
repeat0 <- 1:10
sample_name <- paste('s', repeat0, sep = "")
sample_result <- list()
for (j in sample_name){
  print(j)
  core_all <- vector()
  pan_all <- vector()
  accessory_all <- vector()
  for (i in strain_num){
    strain_choose <- sample(strain_332, i)
    rxn_set <- select(rxn_existence0, strain_choose)
    sum_set <-apply(rxn_set,1,sum)
    rxn_set_new <- data_frame(rxnID=rxn_existence$rxnID, num= sum_set)
    core <-rxn_set_new[rxn_set_new$num==i,]
    pan <- rxn_set_new[rxn_set_new$num >=1,]
    accessory <- setdiff(pan$rxnID, core$rxnID)
    core_all <- c(core_all, length(core$rxnID))
    pan_all <- c(pan_all, length(pan$rxnID))
    accessory_all <- c(accessory_all, length(accessory))
  }
  reactome <- data_frame(num= strain_num, pan=pan_all, core=core_all, accessory=accessory_all)
  sample_result[[j]] <- reactome
}


# calculate the average value
strain_num <- 0
pan_av <- 0
core_av <- 0
acce_av <-0
for (j in sample_name){
  print(j)
  strain_num <- strain_num +  sample_result[[j]]$num
  pan_av <- pan_av + sample_result[[j]]$pan
  core_av <- core_av + sample_result[[j]]$core
  acce_av <- acce_av + sample_result[[j]]$accessory
}


strain_num <- strain_num/10
pan_av <- pan_av/10
core_av <- core_av/10
acce_av <- acce_av/10


reactome <- data_frame(num= strain_num, pan=pan_av, core=core_av, accessory=acce_av)

number_ticks <- function(n) {
  function(limits) pretty(limits, n)
}
ggplot(reactome, aes(num)) + 
  geom_line(aes(y = core, colour = "Core gene")) + 
  geom_line(aes(y = pan, colour = "Pan gene")) +
  geom_line(aes(y = accessory, colour = "Accessory gene")) +
  
  scale_x_continuous(breaks = number_ticks(10)) +
  scale_y_continuous(limits = c(0, 4000), breaks = number_ticks(4)) +
  labs(x="Sampled strain number",y="Reaction number") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text = element_text(size=15)) +
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=24, family="Arial") ) +
  ggtitle('') +
  theme(legend.position="none") #+   theme(panel.background = element_rect(fill = "white", color="black", size = 1)) 










#--------------------------------------
# other tasks
# Build the id conversion between biocyc and keeg reaction ID
# Rxn_biocyc <- read_excel("data/All_reactions_of_MetaCyc.xlsx")
# Rxn_biocyc$KEGGID <- NA

#for (i in seq_along(Rxn_biocyc$Reaction)) {
#  s0 <- Rxn_biocyc$KEGG[i]
#  s01 <- unlist(str_split(s0, " // "))
#  s02 <- str_extract_all(s01, "(?<=rn:)[^.]+(?=\\'>)")
#  s02 <- unlist(s02)
#  s02 <- paste0(s02,collapse = ";")
#  print(paste(i, s02, sep = " ==>"))
#  Rxn_biocyc$KEGGID[i] <- s02
#}