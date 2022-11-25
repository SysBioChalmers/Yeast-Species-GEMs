# Annote the pan-reactome
# Revised by Hongzhong 2019-11-03

# load library
library(readxl)
library(stringr)
library(tidyverse)
library(hongR)
source('function_general.R')



# note: chem_prop_metanetx.tsv file is missing, so part of code report errors.
# currently, this missing datasets could be downloaded from https://www.metanetx.org/ftp/3.2/






# read the data
union_rxn_ann0 <- read.table("result/pan_reactome.txt", header=TRUE, sep = "\t", stringsAsFactors = FALSE)

# mnx
union_rxn_mnx <- union_rxn_ann0[str_detect(union_rxn_ann0$rxnID, 'MNXR'),]
union_rxn_mnx <- getRxnInfFromMNX(union_rxn_mnx, union_rxn_mnx$rxnID)
union_rxn_mnx <- getRxnReversibilty(rxn_frame=union_rxn_mnx, MNXID_list=union_rxn_mnx$rxnID)
# establish mapping between met and each rxnID
# metabolites unify based on MNXID and MNX equation
rxn_final0 <- data.frame(ID0=union_rxn_mnx$rxnID, Equation=union_rxn_mnx$Equation_MNX)
rxn_met <- splitRxnToMetabolite(rxn_final0, sep0 = " = ")
met_mnx <- metStandard.mnx(rxn_met, rxn_met$MetID)
# further refine the rxn reversiblity information for the reactions with MNXID
result_refine <- refineRxnReversiblity(rxn_ini = union_rxn_mnx, met_ini = met_mnx)
union_rxn_mnx <- result_refine[['rxn']]
met_mnx <- result_refine[['met']]
# further refine the rxn and metabolites
rxn_unbalance <- union_rxn_mnx[union_rxn_mnx$balance_MNX =='false' | union_rxn_mnx$balance_MNX =='ambiguous' | union_rxn_mnx$balance_MNX =='NA', ]
met_m_remove <- met_mnx[str_detect(met_mnx$Description, "^a ") | str_detect(met_mnx$Description, "protein") |str_detect(met_mnx$Description, "^an "),]
rxn_rm <- unique(c(rxn_unbalance$rxnID, met_m_remove$ID))
# get the refind the rxn with MNXID
union_rxn_mnx_refine <- union_rxn_mnx[!(union_rxn_mnx$rxnID %in%rxn_rm), ]
rxn_remove_all <- setdiff(union_rxn_mnx$rxnID, union_rxn_mnx_refine$rxnID)
met_mnx <- met_mnx[!(met_mnx$ID %in% rxn_remove_all),]





# biocyc
# id_mapping_metacyc <- read_excel("data/All_compounds_of_MetaCyc.xlsx")
union_rxn_biocyc <- union_rxn_ann0[str_detect(union_rxn_ann0$rxnID, 'metacyc'),]
union_rxn_biocyc <- getBiocycRxnInfFromBiocycID(union_rxn_biocyc, union_rxn_biocyc$rxnID)
# establish mapping between met and each rxnID
met_biocyc <- metStandard.biocyc (rxn_inf = union_rxn_biocyc)
# here we further remove rxns which contains met without mnxid or contains a general metabolite formula
met_b_remove <- met_biocyc[is.na(met_biocyc$MetID) | str_detect(met_biocyc$Description, "^a "),]
rxn_b_remove <- unique(met_b_remove$ID)
union_rxn_biocyc_refine <- union_rxn_biocyc[!(union_rxn_biocyc$rxnID %in%rxn_b_remove), ]
met_biocyc <- met_biocyc[!(met_biocyc$ID %in%rxn_b_remove), ]



# kegg
# for the rxn with keegID we firstly find rheaID to check the reversibility
# for the keggID which can't found MNXID, it difficult to find the biocycID
# so we only choose the keggid with rheaID
# then we can standardize the metabolites based on the reaction annotation from Rhea
union_rxn_kegg <- union_rxn_ann0[str_detect(union_rxn_ann0$rxnID, 'kegg'),]
union_rxn_kegg <- getRxnInfFromKEGG(union_rxn_kegg, union_rxn_kegg$rxnID)
union_rxn_kegg <- getRheaRxnInfFromKEGGID(union_rxn_kegg, union_rxn_kegg$rxnID)
union_rxn_kegg_refine <- union_rxn_kegg[!is.na(union_rxn_kegg$formula_rhea), ]
# establish mapping between met and each rxnID
rxn_kegg <- union_rxn_kegg_refine[, c("rxnID", "equation_rhea")]
colnames(rxn_kegg) <- c("ID0", "Equation")
rxn_kegg1 <- splitRxnToMetabolite(reationFrame = rxn_kegg, sep0 = "=", source0 = "kegg")
# find MNXID of metabolite based on keggid
chem_xref <- read_tsv("data/chem_xref.tsv")
chem_xref_chebi <- chem_xref[str_detect(chem_xref$XREF,'chebi'),]
chem_xref_chebi$XREF <- str_replace_all(chem_xref_chebi$XREF, "chebi:", "CHEBI:")

rxn_kegg1$MetID <- getMultipleReactionFormula(chem_xref_chebi$MNX_ID, chem_xref_chebi$XREF, rxn_kegg1$MetID)
met_kegg <- rxn_kegg1



# combine all the metabolite information
union_rxn_mnx0 <- select(union_rxn_mnx, rxnID, Occur_num, formula_MNX, EC_MNX, keggID, biocycID)
colnames(union_rxn_mnx0) <- c('rxnID','Occur_num','formula','EC','keggID','biocycID')
union_rxn_mnx0$annote_method <- 'MNX'
union_rxn_biocyc0 <- select(union_rxn_biocyc, rxnID, Occur_num, formula_biocyc, EC_biocyc)
colnames(union_rxn_biocyc0) <- c('rxnID','Occur_num','formula','EC')
union_rxn_biocyc0$keggID <- NA
union_rxn_biocyc0$biocycID <- union_rxn_biocyc0$rxnID
union_rxn_biocyc0$annote_method <- 'biocyc'
union_rxn_kegg0 <- select(union_rxn_kegg, rxnID, Occur_num, formula_kegg, EC_kegg)
colnames(union_rxn_kegg0) <- c('rxnID','Occur_num','formula','EC')
union_rxn_kegg0$keggID <- union_rxn_kegg$rxnID
union_rxn_kegg0$biocycID <- NA
union_rxn_kegg0$annote_method <- 'kegg'

pan_reactome <- rbind.data.frame(union_rxn_mnx0, union_rxn_biocyc0, union_rxn_kegg0)

write.table(pan_reactome, "result/draft_pan_reactome_for_yeast_species.txt", row.names = FALSE, sep = "\t")
write.table(met_mnx, "result/met_mnx_pan_reactome_for_yeast_species.txt", row.names = FALSE, sep = "\t")
write.table(met_kegg, "result/met_kegg_pan_reactome_for_yeast_species.txt", row.names = FALSE, sep = "\t")
write.table(met_biocyc, "result/met_biocyc_pan_reactome_for_yeast_species.txt", row.names = FALSE, sep = "\t")
















# part 2
# pathway mapping
pan_reactome <- read.table('result/draft_pan_reactome_for_yeast_species.txt', header = TRUE, stringsAsFactors = FALSE)
rxn_pathway <- read.table('data/kegg/reactionID_pathway.txt', header = FALSE, stringsAsFactors = FALSE)
pathway_list_kegg <- read_excel("data/kegg/pathway_list_kegg.xlsx", col_names = FALSE)
colnames(pathway_list_kegg) <- c('ID','pathway')
rxn_pathway <- merge(rxn_pathway, pathway_list_kegg, by.x = 'V2', by.y = 'ID', all.x = TRUE)
rxn_pathway$V1 <- str_replace_all(rxn_pathway$V1, 'rn:', 'kegg:')
rxn_pathway <- filter(rxn_pathway, str_detect(rxn_pathway$V2,'map'))

rxn_pathway <- select(rxn_pathway, V1, pathway)
pan_reactome0 <- merge(pan_reactome, rxn_pathway, by.x = 'keggID', by.y = 'V1', all.x = TRUE)
pan_reactome0 <- filter(pan_reactome0, !is.na(pathway)) %>% filter(., pathway !='Metabolic pathways')

# plot box plot
pan_reactome1 <- select(pan_reactome0,pathway,Occur_num)
pathway_with_reaction_num <- as.data.frame(table(pan_reactome1$pathway))
colnames(pathway_with_reaction_num) <- c('pathway','rxn_num')
pathway_with_reaction_num0 <- filter(pathway_with_reaction_num , rxn_num >=10 & rxn_num < 100) # remove general subsystems or subsystems with too few reactions
Result0 <- pan_reactome1[pan_reactome1$pathway %in% pathway_with_reaction_num0$pathway,]
Result0$pathway <- as.factor(Result0$pathway)


# set the factor level by the mean value of occur number
Factor <-Result0 %>% group_by(pathway) %>% summarise(median=mean(Occur_num))
Factor <- Factor[order(Factor$median),]
Result0$pathway <-factor(Result0$pathway, levels=Factor$pathway)

ggplot(data=Result0, aes(x=pathway, y=Occur_num)) + geom_boxplot(fill='grey93') +
  xlab('') + ylab('') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
  theme(legend.position = c(0.85, 0.2)) +
  theme(axis.text=element_text(size=10, family="Arial"),
        axis.title=element_text(size=12,family="Arial"),
        legend.text = element_text(size=10, family="Arial")) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) #+


# calculate the number of fraction of accessor and core
for (i in seq_along(Factor$pathway)){
  print(i)
  s1 <- Result0[Result0$pathway==Factor$pathway[i],]
  core_num <- filter(s1, Occur_num > 332*0.95)
  acc_num <-  filter(s1, Occur_num < 332*0.95)
  core_len <- length(core_num$Occur_num)
  acc_len <- length(acc_num$Occur_num)
  Factor$Core[i] <- core_len/(core_len + acc_len)
  Factor$Acc[i] <-  acc_len/(core_len + acc_len)
}

Factor$pathway <-factor(Factor$pathway, levels=Factor$pathway)
ggplot(data=Factor, aes(x=pathway, y=Core)) + geom_point(colour = "red") +
  xlab('') + ylab('') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
  theme(legend.position = c(0.85, 0.2)) +
  theme(axis.text=element_text(size=10, family="Arial"),
        axis.title=element_text(size=12,family="Arial"),
        legend.text = element_text(size=10, family="Arial")) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) +
  geom_point(data = Factor,
          aes(x=pathway, y=Acc), size = 2, colour = "blue")
