# Classification and standardization of the new reactions from different sources: RAVEN, kegg and eggnog
# Revised by Hongzhong 2019-8-15

# load library
source('function_general.R')


# merge the panID and reactionID from RAVEN_kegg, RAVEN_biocyc, kegg_web and eggnog_web
union_rxn_ann0 <- mergePanID()

# we can split reactions into three types: with MNXID, with keggID, with biocycID
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




# step 3 
# possible we can remove some reactions based on the concept of dead end metabolites
# merge the new rxns with yeast8
# check the number of reactions connected with each each metabolite
# then we can regard the metabolites with only one reaction as the dead end metabolites if no exchange reactions is connected further
rxn_met_yeast8 <- rxnMetMappingYeast8()
# combine the rxn_met_new with rxn_met_yeast8
rxn_met_new <- met_mnx[, c('ID','MetID')]
colnames(rxn_met_new) <- c('ID','MNXID')
rxn_met_new2 <- met_biocyc[, c('ID','MetID')]
colnames(rxn_met_new2) <- c('ID','MNXID')
rxn_met_new3 <- met_kegg[, c('ID','MetID')]
colnames(rxn_met_new3) <- c('ID','MNXID')

rxn_met_new_all <- rbind.data.frame(rxn_met_new, rxn_met_new2, rxn_met_new3)
rxn_met_mix <- rbind.data.frame(rxn_met_yeast8, rxn_met_new, rxn_met_new2, rxn_met_new3)
all_new_metID <- setdiff(rxn_met_new_all$MNXID, rxn_met_yeast8$MNXID)
all_metID <- union(rxn_met_new_all$MNXID, rxn_met_yeast8$MNXID)
# calculate how many reactions connect with each new metabolites
rxn_num_of_met <- data.frame(MNXID=all_metID, stringsAsFactors = FALSE)
for (i in 1:nrow(rxn_num_of_met)){
  print(i)
  m0 <- rxn_num_of_met$MNXID[i]
  rxn0 <- rxn_met_mix[rxn_met_mix$MNXID==m0, ]
  num0 <- length(unique(rxn0$ID))
  rxn_num_of_met$rxn_num[i] <- num0
}

rxn_num_of_met$type[rxn_num_of_met$MNXID %in% all_new_metID] <- 'metabolite_new'
rxn_num_of_met$type[!(rxn_num_of_met$MNXID %in% all_new_metID)] <- 'metabolite_Yeast8'

#here for the metabolite only occured in one new reaction is defined as potential dead-end met
dead_end_met <-filter(rxn_num_of_met, rxn_num == 1 & type=='metabolite_new')

# calculate the number of dead_end_met from each new reactions
new_rxn_all <- data.frame(rxnID = unique(rxn_met_new_all$ID), stringsAsFactors = FALSE)
for (i in 1: nrow(new_rxn_all)){
  r0 <- new_rxn_all$rxnID[i]
  rxn_met0 <- filter(rxn_met_new_all, ID==r0)
  met0 <- unique(rxn_met0$MNXID)
  dead_met_num <- length(which(met0 %in% dead_end_met$MNXID ==TRUE))
  new_rxn_all$dead_met_num[i] <- dead_met_num
}

# then we can integrate the dead end metabolite number informatin with the refined new rxn informatiion
union_rxn_mnx_refine$dead_met_num <- getSingleReactionFormula(new_rxn_all$dead_met_num,new_rxn_all$rxnID,union_rxn_mnx_refine$rxnID)
union_rxn_kegg_refine$dead_met_num <- getSingleReactionFormula(new_rxn_all$dead_met_num,new_rxn_all$rxnID,union_rxn_kegg_refine$rxnID)
union_rxn_biocyc_refine$dead_met_num <- getSingleReactionFormula(new_rxn_all$dead_met_num,new_rxn_all$rxnID,union_rxn_biocyc_refine$rxnID)


# step 4
# summarize the current panYeast reactions and choose the common reactions sets
reportRxnPanYeast <- parseReportedPanYeast()
rxn_common <- intersect(union_rxn_mnx$rxnID, reportRxnPanYeast)
union_rxn_mnx_refine$exist_panYeast <- union_rxn_mnx_refine$rxnID %in% rxn_common


# step 5 
# remove some reactions with dead end metabolite or just choose the reactions existing in history panYeast
union_rxn_mnx_refine <- filter(union_rxn_mnx_refine, dead_met_num <=0 | exist_panYeast == TRUE)
union_rxn_kegg_refine <- filter(union_rxn_kegg_refine, dead_met_num <=0 )
union_rxn_biocyc_refine <- filter(union_rxn_biocyc_refine, dead_met_num <=0 )

met_mnx_refine <- met_mnx[met_mnx$ID %in% union_rxn_mnx_refine$rxnID, ]
met_kegg_refine <- met_kegg[met_kegg$ID %in% union_rxn_kegg_refine$rxnID, ]
met_biocyc_refine <- met_biocyc[met_biocyc$ID %in% union_rxn_biocyc_refine$rxnID, ]


write.table(union_rxn_mnx_refine, "result/new rxn information from MNX database.txt", row.names = FALSE, sep = "\t")
write.table(met_mnx_refine, "result/new met information from MNX database.txt", row.names = FALSE, sep = "\t")




# potential dead end metabolite analysis
# ----------------------------
ggplot(filter(rxn_num_of_met, rxn_num < 25), aes(rxn_num, fill = type, colour = type)) +
  geom_density(alpha = 0.1) +
  labs(x='rxn number connected with each metabolite') +
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=15, family="Arial") ) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 1)) +
  theme(legend.position = c(0.8, 0.2))


# dead end met annotation
colnames(dead_end_met) <- c('MetID','rxn_num','type')
dead_end_met  <- metStandard.mnx(dead_end_met, dead_end_met$MetID)
