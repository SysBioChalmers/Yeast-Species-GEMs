# give the standard reaction and metabolite annotation information based on reactionID from MNX database

library(readxl)
library(stringr)
library(tidyverse)
library(hongR)
library(ggExtra)


splitAndCombine <- function(gene, rxn,sep0) { 
  # one rxn has several genes, this function was used to splite the genes
   gene <- str_split(gene, sep0)
   tt<- length(gene)
   gene0 <- list()
   for (i in 1:tt){
    gene0[[i]] <- paste(rxn[i], gene[[i]], sep = "@@@")
  
   }
  
  gene1 <- unique(unlist(gene0))
  gene2 <- str_split(gene1, "@@@" )
  rxnGene <- data.frame(v1=character(length(gene2)),stringsAsFactors = FALSE)
  tt1 <- length(gene2)
  for (j in 1:tt1){
  rxnGene$v1[j] <- gene2[[j]][2]
  rxnGene$v2[j] <- gene2[[j]][1]
  }
 
return(rxnGene)
}


splitRxnToMetabolite <- function(reationFrame, sep0, source0 ='MNX'){
  #reationFrame <- rxn_final0
  #input:
  #dataframe: column 1--ID0; column 2--Equation, which is the reaction formula
  rxn_list <- str_split(reationFrame$Equation, sep0)
  for (i in seq(length(rxn_list))){
    rxn_list[[i]][1] <- paste("reactant",rxn_list[[i]][1],sep = "@@" )
    rxn_list[[i]][1] <- paste(reationFrame$ID0[i],rxn_list[[i]][1],sep = "@" )
    
    rxn_list[[i]][2] <- paste("product",rxn_list[[i]][2],sep = "@@" )
    rxn_list[[i]][2] <- paste(reationFrame$ID0[i],rxn_list[[i]][2],sep = "@" )
  }
  
  rxn_unlist <- unlist(rxn_list)
  rxn0_list <- str_split(rxn_unlist, "@@")
  ss1 <- vector()
  ss2 <- vector()
  for (i in seq(length(rxn0_list))){
    ss1[i] <- rxn0_list[[i]][1]
    ss2[i] <- rxn0_list[[i]][2]
  }
  ss2_list <- str_split(ss2, " \\+ ")
  
  for (i in seq(length(rxn0_list))){
    ss2_list[[i]] <- paste(ss1[[i]], ss2_list[[i]], sep = "@@")
  }
  
  ss2_unlist <- unlist(ss2_list)
  ss3_list <- str_split(ss2_unlist, "@@")
  ss4 <- vector()
  ss5 <- vector()
  for (i in seq(length(ss3_list))){
    ss4[i] <- ss3_list[[i]][1]
    ss5[i] <- ss3_list[[i]][2]
  }
  rxn_met <- data.frame(reaction = ss4, MetID = ss5, stringsAsFactors = FALSE)
  rxn_met <- rxn_met %>% separate(.,reaction, into = c('ID','compostion'), sep = "@")
  #remove the coefficient for each metabolites, which will be easy to standardize this metabolite
  # here we should make this step more reasonable
  
  if(source0=='MNX'){
    rxn_met <- rxn_met %>% separate(.,MetID, into = c('coefficient','MetID'), sep = " ")}
  
  else if(source0 =='Yeast8'){
    for(i in 1:nrow(rxn_met)){
      rxn_met$MetID[i] <- str_trim(rxn_met$MetID[i], side = "both")
      if(str_detect(rxn_met$MetID[i], '^\\d+\\.*\\d* ')){
        rxn_met$coefficient[i] <- str_extract(rxn_met$MetID[i], '^\\d+\\.*\\d* ') %>% str_trim(.,side = "both")
        rxn_met$MetID[i] <- str_replace_all(rxn_met$MetID[i], '^\\d+\\.*\\d* ','')
      } else{
        rxn_met$coefficient[i] <- 1
        rxn_met$MetID[i] <- rxn_met$MetID[i]
      }
      
    }
    rxn_met$coefficient <- as.numeric(rxn_met$coefficient)
  }
  
  else{
    for(i in 1:nrow(rxn_met)){
      rxn_met$MetID[i] <- str_trim(rxn_met$MetID[i], side = "both")
      if(str_detect(rxn_met$MetID[i], '^[:digit:] ')){
        rxn_met$coefficient[i] <- str_extract(rxn_met$MetID[i], '^[:digit:] ') %>% str_trim(.,side = "both")
        rxn_met$MetID[i] <- str_replace_all(rxn_met$MetID[i], '^[:digit:] ','')
      } else{
        rxn_met$coefficient[i] <- 1
        rxn_met$MetID[i] <- rxn_met$MetID[i]
      }
      
    }
    rxn_met$coefficient <- as.numeric(rxn_met$coefficient)
    
  }
  
  #rxn_met$MetID <- str_replace_all(rxn_met$MetID, "[:digit:] ","") %>%
  #   str_replace_all(.,"\\(n\\)","") %>%
  #   str_replace_all(.,"\\(n\\+1\\)","") %>%
  #   str_replace_all(.,"\\(n\\-2\\)","")
  rxn_met$MetID <- str_replace_all(rxn_met$MetID, "@MNXD[:alnum:]","")
  rxn_met$MetID <- str_replace_all(rxn_met$MetID, "@BOUNDARY","")
  return(rxn_met)
}


mergePanID <- function() {
  # using this function, we can:
  # merge the panID and reactionID from RAVEN_kegg, RAVEN_biocyc, kegg_web and eggnog_web
  
  # RAVEN biocyc
  # newRxn_biocyc <- read.table("data/newRxn_biocyc_RAVEN.txt", header= TRUE, stringsAsFactors = FALSE)
  newRxn_biocyc <- read.table("data/newRxn_biocyc_RAVEN_55_110.txt", header = TRUE, stringsAsFactors = FALSE)
  
  newRxn_biocyc$MNXID <- findRxnMNXid(rxnID = newRxn_biocyc$ID, id_type = "metacyc")
  newRxn_biocyc <- getRxnInfFromMNX(newRxn_biocyc, newRxn_biocyc$MNXID)
  
  # RAVEN KEGG
  newRxn_kegg <- read.table("data/newRxn_kegg_RAVEN.txt", header = TRUE, stringsAsFactors = FALSE)
  newRxn_kegg$MNXID <- findRxnMNXid(rxnID = newRxn_kegg$ID, id_type = "kegg")
  newRxn_kegg <- getRxnInfFromMNX(newRxn_kegg, newRxn_kegg$MNXID)
  
  # KEGG and eggnog web services
  newRxn_kegg_eggnog <- read.table("data/newRxn_all based on kegg and eggnog annotation.txt", header = TRUE, stringsAsFactors = FALSE)
  newRxn_kegg_eggnog$rxnID <- str_replace_all(newRxn_kegg_eggnog$rxnID, "rn:", "")
  newRxn_kegg_eggnog$MNXID <- findRxnMNXid(rxnID = newRxn_kegg_eggnog$rxnID, id_type = "kegg")
  newRxn_kegg_eggnog <- getRxnInfFromMNX(newRxn_kegg_eggnog, newRxn_kegg_eggnog$MNXID)
  
  rxn_kegg_web <- newRxn_kegg_eggnog[str_detect(newRxn_kegg_eggnog$type, "kegg"), ]
  rxn_eggnog_web <- newRxn_kegg_eggnog[str_detect(newRxn_kegg_eggnog$type, "eggnog"), ]
  
  kegg_web <- unique(rxn_kegg_web$MNXID)
  eggnog_web <- unique(rxn_eggnog_web$MNXID)
  RAVEN_kegg <- unique(newRxn_kegg$MNXID)
  RAVEN_biocyc <- unique(newRxn_biocyc$MNXID)
  
  
  # merge the rxn from different source
  union_rxn <- unique(c(kegg_web, eggnog_web, RAVEN_kegg, RAVEN_biocyc))
  union_rxn_ann <- data.frame(rxnID = union_rxn, stringsAsFactors = FALSE)
  union_rxn_ann$source <- NA
  union_rxn_ann$panID_kegg_web <- NA
  union_rxn_ann$panID_eggnog_web <- NA
  union_rxn_ann$panID_RAVEN_kegg <- NA
  union_rxn_ann$panID_RAVEN_biocyc <- NA
  union_rxn_ann$panID_inter <- NA # the interaction panID sets from at least two different sources
  union_rxn_ann$panID_union <- NA # the union panID sets from different sources
  
  # summarize the panID
  # integrate the source
  E1 <- union_rxn_ann$rxnID %in% kegg_web
  E2 <- union_rxn_ann$rxnID %in% eggnog_web
  E3 <- union_rxn_ann$rxnID %in% RAVEN_kegg
  E4 <- union_rxn_ann$rxnID %in% RAVEN_biocyc
  source_merge <- vector()
  for (i in 1:length(E1)) {
    if (E1[i]) {
      s1 <- "kegg_web"
    } else {
      s1 <- NA
    }
    if (E2[i]) {
      s2 <- "eggnog_web"
    } else {
      s2 <- NA
    }
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
    s <- paste(s1, s2, s3, s4, sep = ";")
    source_merge <- c(source_merge, s)
  }
  
  source_merge <- str_replace_all(source_merge, "NA;", "")
  source_merge <- str_replace_all(source_merge, ";NA", "")
  union_rxn_ann$source <- source_merge
  # union_rxn_ann0 <- filter(union_rxn_ann, balance_MNX=='true') # 1275 balanced reactions
  # union_rxn_ann0$source_num <- str_count(union_rxn_ann0$source, ";") # 605 reactions from at least two evidences
  union_rxn_ann$panID_kegg_web <- getMultipleReactionFormula(rxn_kegg_web$panID_kegg, rxn_kegg_web$MNXID, union_rxn_ann$rxnID)
  union_rxn_ann$panID_eggnog_web <- getMultipleReactionFormula(rxn_eggnog_web$panID_eggnog, rxn_eggnog_web$MNXID, union_rxn_ann$rxnID)
  union_rxn_ann$panID_RAVEN_kegg <- getMultipleReactionFormula(newRxn_kegg$GENE.ASSOCIATION, newRxn_kegg$MNXID, union_rxn_ann$rxnID)
  union_rxn_ann$panID_RAVEN_biocyc <- getMultipleReactionFormula(newRxn_biocyc$GENE.ASSOCIATION, newRxn_biocyc$MNXID, union_rxn_ann$rxnID)
  
  # unify the panID for each reactions
  g_union_all <- vector()
  g_inter_all <- vector()
  for (i in 1:nrow(union_rxn_ann)) {
    print(i)
    g1 <- union_rxn_ann$panID_kegg_web[i]
    g1 <- unlist(str_split(g1, ";"))
    g1 <- str_trim(g1, side = "both")
    
    g2 <- union_rxn_ann$panID_eggnog_web[i]
    g2 <- unlist(str_split(g2, ";"))
    g2 <- str_trim(g2, side = "both")
    
    g3 <- union_rxn_ann$panID_RAVEN_kegg[i]
    g3 <- unlist(str_split(g3, " or "))
    g3 <- str_trim(g3, side = "both")
    
    g4 <- union_rxn_ann$panID_RAVEN_biocyc[i]
    g4 <- unlist(str_split(g4, " or "))
    g4 <- str_trim(g4, side = "both")
    
    
    g_union <- Reduce(union, list(g1, g2, g3, g4)) # the union panID sets from different sources
    g_union <- g_union[!is.na(g_union)]
    g_union_string <- paste(g_union, collapse = ";")
    e1 <- g_union %in% g1
    e2 <- g_union %in% g2
    e3 <- g_union %in% g3
    e4 <- g_union %in% g4
    e_add <- e1 + e2 + e3 + e4
    g_inter <- g_union[which(e_add >= 2)] # here we choose panID from at least two sources
    if (length(g_inter)) {
      g_inter_string <- paste(g_inter, collapse = ";")
    } else {
      g_inter_string <- NA
    }
    g_union_all <- c(g_union_all, g_union_string)
    g_inter_all <- c(g_inter_all, g_inter_string)
  }
  
  union_rxn_ann$panID_inter <- g_inter_all
  union_rxn_ann$panID_union <- g_union_all
  
  return(union_rxn_ann)
}



parseReportedPanYeast <- function() {
  # This function is used to parse the reaction from the reported panYeast doi: http://dx.doi.org/10.1101/412593.
  # input a early panYeast
  aybraham <- read_excel("data/aybraham.xlsx", sheet = "reactions")
  aybraham2 <- aybraham[2:nrow(aybraham), ]
  colnames(aybraham2) <- aybraham[1, ]
  
  # find mnxid based on biggid
  rxnID_mapping <- data.frame(rxnID = aybraham2$"Rxn name", stringsAsFactors = FALSE)
  ID1 <- findRxnMNXid(rxnID = rxnID_mapping$rxnID, id_type = "bigg")
  ID1 <- str_replace_all(ID1, "bigg:", "")
  rxnID_mapping$MNXID1 <- ID1
  rxnID_mapping$MNXID1[!str_detect(rxnID_mapping$MNXID1, "MNXR")] <- NA
  
  
  rxnID_mapping$keggID <- aybraham2$kegg %>%
    str_replace_all(., ":proton_less", "") %>%
    str_replace_all(., ":ammonium", "") %>%
    str_replace_all(., ":cofactor", "")
  ID2 <- findRxnMNXid(rxnID = rxnID_mapping$keggID, id_type = "kegg")
  ID2 <- str_replace_all(ID2, "kegg:", "")
  rxnID_mapping$MNXID2 <- ID2
  rxnID_mapping$MNXID2[!str_detect(rxnID_mapping$MNXID2, "MNXR")] <- NA
  
  
  rxnID_mapping$rheaID <- aybraham2$rhea %>%
    str_replace_all(., ":cofactor", "") %>%
    str_replace_all(., ":proton_less", "") %>%
    str_replace_all(., "rhea:", "")
  ID3 <- findRxnMNXid(rxnID = rxnID_mapping$rheaID, id_type = "rhea")
  rxnID_mapping$MNXID3 <- ID3
  rxnID_mapping$MNXID3[!str_detect(rxnID_mapping$MNXID3, "MNXR")] <- NA
  
  
  
  rxnID_mapping$biocycID <- aybraham2$metacyc %>%
    str_replace_all(., ":cofactor", "") %>%
    str_replace_all(., ":coeff", "") %>%
    str_replace_all(., ":proton_less", "")
  ID4 <- findRxnMNXid(rxnID = rxnID_mapping$biocycID, id_type = "metacyc")
  rxnID_mapping$MNXID4 <- ID4
  rxnID_mapping$MNXID4[!str_detect(rxnID_mapping$MNXID4, "MNXR")] <- NA
  
  # some rxn have multiple external ID
  s1 <- rxnID_mapping[, c("rxnID", "keggID")]
  s10 <- splitAndCombine(s1$keggID, s1$rxnID, sep0 = "\\|")
  s11 <- findRxnMNXid(rxnID = s10$v1, id_type = "kegg")
  s11[!str_detect(s11, "MNXR")] <- NA
  
  
  s2 <- rxnID_mapping[, c("rxnID", "rheaID")]
  s20 <- splitAndCombine(s2$rheaID, s2$rxnID, sep0 = ";")
  s21 <- findRxnMNXid(rxnID = s20$v1, id_type = "rhea")
  s21[!str_detect(s21, "MNXR")] <- NA
  
  
  s3 <- rxnID_mapping[, c("rxnID", "biocycID")]
  s30 <- splitAndCombine(s3$biocycID, s3$rxnID, sep0 = ";")
  s31 <- findRxnMNXid(rxnID = s30$v1, id_type = "metacyc")
  s31[!str_detect(s31, "MNXR")] <- NA
  
  mnx_all <- unique(c(rxnID_mapping$MNXID1, s11, s21, s31))
  mnx_all <- mnx_all[!is.na(mnx_all)]
  
  return(mnx_all)
}



rxnMetMappingYeast8 <- function() {
  # This function is used to establish the rxnID in yeast8 with its metabolite MNXID
  rxn_yeast8 <- read_excel("data/rxn_yeast8_2019_8.xlsx")
  met_yeast8 <- read_excel("data/met_yeast8_2019_8.xlsx")
  rxn_yeast8 <- rxn_yeast8[, c("equation", "rxnID")]
  rxn_yeast8$equation <- str_replace_all(rxn_yeast8$equation, "-->", "<=>")
  rxn_yeast80 <- data.frame(ID0 = rxn_yeast8$rxnID, Equation = rxn_yeast8$equation)
  rxn_yeast8_split <- splitRxnToMetabolite(rxn_yeast80, sep0 = "<=>", source0 = "Yeast8")
  # correct the coefficient
  rxn_yeast8_split1 <- rxn_yeast8_split[str_detect(rxn_yeast8_split$MetID, " s_"), ]
  rxn_yeast8_split2 <- rxn_yeast8_split[!str_detect(rxn_yeast8_split$MetID, " s_"), ]
  rxn_yeast8_split1 <- rxn_yeast8_split1 %>% separate(., MetID, into = c("coefficient", "MetID"), sep = " ")
  rxn_yeast8_split2 <- rxn_yeast8_split2[rxn_yeast8_split2$MetID != "", ]
  rxn_yeast8_split_c <- rbind.data.frame(rxn_yeast8_split2, rxn_yeast8_split1)
  # get the MNXID for each metabolite in the Yeast8
  rxn_yeast8_split_c$MNXID <- getSingleReactionFormula(met_yeast8$MNXID_new, met_yeast8$m_name, rxn_yeast8_split_c$MetID)
  rxn_met <- rxn_yeast8_split_c[, c("ID", "MNXID")]
  return(rxn_met)
}



findRxnMNXid <- function(rxnID, id_type = "metacyc") {
  # This function is used to find the MNXid based on one of other reactionID, like kegg id or biocyc id
  # Input
  # rxnID: A vector of rxnID
  # id_type: The kind of ID source, like seed, metacyc, kegg, rhea
  
  # here we found that the latest id mapping from mnx is not right for some KO id from kegg database
  # so we used the data downloaded in 2018
  # Revised by Hongzhong 2019-8-12
  
  # reac_xref <- read_tsv("data/reac_xref_2019_8.tsv") # using the latest version from MNX database
  reac_xref <- read_tsv("data/reac_xref.tsv") # using the latest version from MNX database
  mnx_special <- reac_xref[str_detect(reac_xref$XREF, id_type), ]
  id_type0 <- paste(id_type, ":", sep = "")
  mnx_special$XREF <- str_replace_all(mnx_special$XREF, id_type0, "")
  MNXID <- getMultipleReactionFormula(mnx_special$MNX_ID, mnx_special$XREF, rxnID)
  newID <- c()
  for (x in 1:length(MNXID)) {
    print(x)
    if (!is.na(MNXID[x])) {
      newID[x] <- MNXID[x]
    } else {
      newID[x] <- paste(id_type0,rxnID[x],sep = "")
    }
  }
  return(newID)
}




findRxnMNXidFromEC <- function(EC_list){
  # This function is used to get the detailed reaction annotation from MNX database based on the EC number
  # Input
  # EC_list: A vector contains EC number
  reac_metnet <- read_tsv('data/reac_prop_metanetx.tsv')
  #EC_list <- newEC_eggnog$EC
  ec_rxn <- reac_metnet[,c('MNX_ID','EC')]
  ec_rxn <- ec_rxn[!is.na(ec_rxn$EC),]
  ec_rxn0 <- splitAndCombine(ec_rxn$EC, ec_rxn$MNX_ID, sep0 = ";")
  ec_rxn0$v1 <- str_trim(ec_rxn0$v1)
  newID <- getMultipleReactionFormula(ec_rxn0$v2 ,ec_rxn0$v1, EC_list)
  return(newID)
}



getRxnInfFromMNX <- function(rxn_frame, MNXID_list) {
  # This function is used to get the detailed reaction annotation from MNX database based on the mnxID of the reaction
  # Input
  # rxn_frame: A dataframe contains a column of MNXID
  # MNID_list: A vector contains the above MNXID in the rxn_frame

  # rxn_frame <- union_rxn_mnx # for test
  # MNXID_list <- union_rxn_mnx$rxnID # for test
  reac_metnet <- read_tsv("data/reac_prop_metanetx.tsv")
  # find the the name from kegg and seed data base
  seed_id_mapping <- read_table2("data/model_seed/Aliases/Unique_ModelSEED_Reaction_Aliases.txt")
  seed_rxn <- read_tsv("data/model_seed/reactions.tsv")
  seed_id_mapping <- filter(seed_id_mapping, str_detect(seed_id_mapping$External, "metanetx.reaction"))
  rxn_frame$seedID <- getSingleReactionFormula(seed_id_mapping$ModelSEED, seed_id_mapping$ID, MNXID_list)

  # reaction summary in kegg
  rxn_formula_kegg <- read.delim2("data/reaction_kegg summary_2019_8.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  f1 <- rxn_formula_kegg$V2
  rxn_name <- vector()
  rxn_formula0 <- vector()
  for (i in 1:length(f1)) {
    s0 <- f1[i]
    s0 <- unlist(str_split(s0, "; "))
    n1 <- s0[!str_detect(s0, "=")]
    n2 <- s0[str_detect(s0, "=")]
    n1 <- paste(n1, collapse = ";")
    n2 <- paste(n2, collapse = ";")
    rxn_name <- c(rxn_name, n1)
    rxn_formula0 <- c(rxn_formula0, n2)
  }
  rxn_formula_kegg$formula <- rxn_formula0
  rxn_formula_kegg$rxn_name <- rxn_name
  colnames(rxn_formula_kegg) <- c("rxnID", "inf", "reaction", "rxn_name")
  rxn_formula_kegg$rxnID <- str_replace_all(rxn_formula_kegg$rxnID, "rn:", "kegg:")
  rxn_formula_kegg$reaction <- str_trim(rxn_formula_kegg$reaction, side = "both")
  # find the keggID based on MNXID, then find the reaction name
  reac_xref <- read_tsv("data/reac_xref.tsv")
  mnx_kegg <- reac_xref[str_detect(reac_xref$XREF, "kegg"), ]
  rxn_frame$keggID <- getSingleReactionFormula(mnx_kegg$XREF, mnx_kegg$MNX_ID, MNXID_list)
  # find the reaction name
  rxn_frame$rxn_name_seed <- getSingleReactionFormula(seed_rxn$name, seed_rxn$id, rxn_frame$seedID)
  rxn_frame$rxn_name_kegg <- getSingleReactionFormula(rxn_formula_kegg$rxn_name, rxn_formula_kegg$rxnID, rxn_frame$keggID)
  rxn_frame$formula_MNX <- getSingleReactionFormula(reac_metnet$Description, reac_metnet$MNX_ID, MNXID_list)
  rxn_frame$balance_MNX <- getSingleReactionFormula(reac_metnet$Balance, reac_metnet$MNX_ID, MNXID_list)
  rxn_frame$Equation_MNX <- getMultipleReactionFormula(reac_metnet$Equation, reac_metnet$MNX_ID, MNXID_list)
  rxn_frame$EC_MNX <- getMultipleReactionFormula(reac_metnet$EC, reac_metnet$MNX_ID, MNXID_list)
  return(rxn_frame)
}




getRxnReversibilty <- function(rxn_frame, MNXID_list) {
  # This function is used to get the detailed reaction reversibility information from model seed and biocyc database
  # rxn_frame: A dataframe contains a column of MNXID
  # MNID_list: A vector contains the above MNXID in the rxn_frame
  rxn_inf <- rxn_frame
  MNXid <- MNXID_list
  # it should be noted that a MNXID could have multiple metacyc ids
  reac_xref <- read_tsv('data/reac_xref.tsv')
  #reac_xref <- read_tsv('data/reac_xref_2019_8.tsv') # using the latest version from MNX database; seems not good of the new version
  
  mnx_biocyc <- reac_xref[str_detect(reac_xref$XREF, "metacyc"), ]
  
  mnx_biocyc$XREF <- str_replace_all(mnx_biocyc$XREF, "metacyc:", "")
  biocyc <- read.delim2("data/All_instances_of_Reactions_in_MetaCyc.txt", stringsAsFactors = FALSE)
  biocyc_reaction <- read_excel("data/metcyc_rxn_RAVEN.xlsx")
  
  rxn_inf$biocycID <- getSingleReactionFormula(mnx_biocyc$XREF, mnx_biocyc$MNX_ID, MNXid)
  # find the reaction reversiblity from biocyc database
  rxn_inf$reverse_biocyc <- getSingleReactionFormula(biocyc$Reaction.Direction, biocyc$Reactions, rxn_inf$biocycID)
  rxn_inf$formula_biocyc <- getSingleReactionFormula(biocyc_reaction$formula, biocyc_reaction$Reactions, rxn_inf$biocycID)
  
  # check the reveserbility based on modelSeed database
  # get the modelSEEDid based on the MNXID
  mnx_seed <- reac_xref[str_detect(reac_xref$XREF, "seed"), ]
  mnx_seed$XREF <- str_replace_all(mnx_seed$XREF, "seed:", "")
  seed_id_mapping <- read_table2('data/model_seed/Aliases/Unique_ModelSEED_Reaction_Aliases.txt')
  seed_rxn <- read_tsv('data/model_seed/reactions.tsv')
  seed_id_mapping <- filter(seed_id_mapping, str_detect(seed_id_mapping$External, "metanetx.reaction"))
  #seed <- read_excel("data/ModelSEED-reactions-db.xls")
  #rxn_inf$seedID <- getSingleReactionFormula(mnx_seed$XREF, mnx_seed$MNX_ID, MNXid)
  rxn_inf$seedID <- getSingleReactionFormula(seed_id_mapping$ModelSEED, seed_id_mapping$ID, MNXid)
  # find the reaction reversiblity from seed database
  rxn_inf$seed_equation <- getSingleReactionFormula(seed_rxn$equation, seed_rxn$id, rxn_inf$seedID)
  rxn_inf$seed_definition <- getSingleReactionFormula(seed_rxn$definition, seed_rxn$id, rxn_inf$seedID)
  rxn_inf$seed_direction <- getSingleReactionFormula(seed_rxn$direction, seed_rxn$id, rxn_inf$seedID)
  rxn_inf$seed_reversibility <- getSingleReactionFormula(seed_rxn$reversibility, seed_rxn$id, rxn_inf$seedID)
  return(rxn_inf)
}


refineRxnReversiblity <- function(rxn_ini = union_rxn_mnx, met_ini = met_mnx) {
  # This function is merge the biocyc and seed information into the mnx reaction equation
  # Using the MNXID of reaction, we can get the the biocycID and seedID and the related reaction reversiblity information
  # So here the biocyc will be the main reference, if the biocyc id cann't be found, then the seed information will be used
  
  # note:
  # some general rules to estimate reaction direction
  # Gibbs free energy change of a reaction
  # ATP-consuming reactions were not allowed to proceed in the reverse direction;
  # reactions present in thermodynamically curated models
  
  # divided it into three types
  rxn_with_biocyc <- rxn_ini[rxn_ini$formula_biocyc != "NA", ]
  rxn_with_only_seed <- rxn_ini[rxn_ini$formula_biocyc == "NA" & rxn_ini$seed_equation != "NA", ]
  rxn_no_biocyc_seed <- rxn_ini[rxn_ini$formula_biocyc == "NA" & rxn_ini$seed_equation == "NA", ]
  
  # step1
  # get the right the bound of reaction with biocycID
  rxn_with_biocyc$reverse_biocyc0 <- str_extract_all(rxn_with_biocyc$formula_biocyc, "<=>")
  rxn_with_biocyc$reverse_correction[rxn_with_biocyc$reverse_biocyc0 == "<=>"] <- 1
  rxn_with_biocyc$reverse_correction[rxn_with_biocyc$reverse_biocyc0 != "<=>"] <- 0
  rxn_with_biocyc <- subset(rxn_with_biocyc, select = -c(reverse_biocyc0))
  
  # get the right definition of metabolite in a reaction
  # establish mapping between met and each rxnID
  met_biocyc0 <- metStandard.biocyc(rxn_inf = rxn_with_biocyc)
  # rewrite reaction using MNX metabolite ID
  rxn_unqiue <- unique(met_biocyc0$ID)
  
  # prepare original met-rxn
  met_mnx1 <- met_ini[met_ini$ID %in% rxn_unqiue, ]
  
  met_mnx1_correction <- met_mnx1[FALSE, ]
  
  for (i in 1:length(rxn_unqiue)) {
    j <- rxn_unqiue[i]
    r1 <- met_mnx1[met_mnx1$ID == j, ]
    r1_biocyc <- met_biocyc0[met_biocyc0$ID == j, ]
    r1$compostion_new <- getSingleReactionFormula(r1_biocyc$compostion, r1_biocyc$MetID, r1$MetID)
    composition_type <- r1$compostion_new == r1$compostion
    print(composition_type)
    if (!all(composition_type)) {
      r1$compostion <- str_replace_all(r1$compostion, "reactant", "m")
      r1$compostion <- str_replace_all(r1$compostion, "product", "reactant")
    }
    
    met_mnx1_correction <- rbind.data.frame(r1, met_mnx1_correction)
    met_mnx1_correction$compostion <- str_replace_all(met_mnx1_correction$compostion, "m", "product")
  }
  
  
  
  # step2
  rxn_with_only_seed <- rxn_ini[rxn_ini$formula_biocyc == "NA" & rxn_ini$seed_equation != "NA", ]
  # first esimate the reaction reversibility
  rxn_with_only_seed$reverse_correction[rxn_with_only_seed$seed_reversibility == "?"] <- "need_manual_check"
  rxn_with_only_seed$reverse_correction[rxn_with_only_seed$seed_reversibility == "="] <- 1
  rxn_with_only_seed$reverse_correction[rxn_with_only_seed$seed_reversibility == "<" | rxn_with_only_seed$seed_reversibility == ">"] <- 0
  rxn_with_only_seed$seed_equation_new <- NA
  # adjust the reaction formula based on the reaction direction and reaction reversibility of seed database
  for (i in 1:nrow(rxn_with_only_seed)) {
    print(i)
    if (rxn_with_only_seed$seed_direction[i] == "=" & rxn_with_only_seed$seed_reversibility[i] == "<") {
      rxn <- rxn_with_only_seed$seed_equation[i] %>% str_replace_all(., "<", "") %>% str_replace_all(., ">", "")
      rxn0 <- str_split(rxn, "=")
      newReact <- str_trim(rxn0[[1]][2], side = "both")
      newProduct <- str_trim(rxn0[[1]][1], side = "both")
      new_rxn <- paste(newReact, newProduct, sep = " = ")
      rxn_with_only_seed$seed_equation_new[i] <- new_rxn
    } else if (rxn_with_only_seed$seed_direction[i] == ">" & rxn_with_only_seed$seed_reversibility[i] == "<") {
      rxn <- rxn_with_only_seed$seed_equation[i] %>% str_replace_all(., "<", "") %>% str_replace_all(., ">", "")
      rxn0 <- str_split(rxn, "=")
      newReact <- str_trim(rxn0[[1]][2], side = "both")
      newProduct <- str_trim(rxn0[[1]][1], side = "both")
      new_rxn <- paste(newReact, newProduct, sep = " = ")
      rxn_with_only_seed$seed_equation_new[i] <- new_rxn
    } else if (rxn_with_only_seed$seed_direction[i] == "<" & rxn_with_only_seed$seed_reversibility[i] != ">") {
      rxn <- rxn_with_only_seed$seed_equation[i] %>% str_replace_all(., "<", "") %>% str_replace_all(., ">", "")
      rxn0 <- str_split(rxn, "=")
      newReact <- str_trim(rxn0[[1]][2], side = "both")
      newProduct <- str_trim(rxn0[[1]][1], side = "both")
      new_rxn <- paste(newReact, newProduct, sep = " = ")
      rxn_with_only_seed$seed_equation_new[i] <- new_rxn
    } else {
      rxn_with_only_seed$seed_equation_new[i] <- rxn_with_only_seed$seed_equation[i] %>% str_replace_all(., ">", "") %>% str_replace_all(., "<", "")
    }
  }
  
  
  
  rxn_inf <- rxn_with_only_seed
  rxn_seed <- rxn_inf[, c("rxnID", "seed_equation_new")]
  colnames(rxn_seed) <- c("ID0", "Equation")
  rxn_seed$Equation <- str_replace_all(rxn_seed$Equation, "<=>", "=") %>%
    str_replace_all(., "\\(", "") %>%
    str_replace_all(., "\\)", "") %>%
    str_replace_all(., "\\[0\\]", "")
  
  rxn_seed1 <- splitRxnToMetabolite(reationFrame = rxn_seed, sep0 = "=", source0 = "biocyc")
  # find MNXID of metabolite based on seedid
  chem_xref <- read_tsv("data/chem_xref.tsv")
  chem_xref_seed <- chem_xref[str_detect(chem_xref$XREF, "seed"), ]
  chem_xref_seed$XREF <- str_replace_all(chem_xref_seed$XREF, "seed:", "")
  rxn_seed1$MetID <- getMultipleReactionFormula(chem_xref_seed$MNX_ID, chem_xref_seed$XREF, rxn_seed1$MetID)
  rxn_seed1 <- metStandard.mnx(met_frame = rxn_seed1, met_MNXID_list = rxn_seed1$MetID)
  
  met_seed0 <- rxn_seed1
  rxn_unique_seed <- unique(rxn_seed1$ID)
  
  # prepare original met-rxn
  met_mnx2 <- met_ini[met_ini$ID %in% rxn_unique_seed, ]
  
  met_mnx2_correction <- met_mnx2[FALSE, ]
  
  for (i in 1:length(rxn_unique_seed)) {
    j <- rxn_unique_seed[i]
    r1 <- met_mnx2[met_mnx2$ID == j, ]
    r1_seed <- met_seed0[met_seed0$ID == j, ]
    r1$compostion_new <- getSingleReactionFormula(r1_seed$compostion, r1_seed$MetID, r1$MetID)
    composition_type <- r1$compostion_new == r1$compostion
    print(composition_type)
    if (!all(composition_type)) {
      r1$compostion <- str_replace_all(r1$compostion, "reactant", "m")
      r1$compostion <- str_replace_all(r1$compostion, "product", "reactant")
    }
    
    met_mnx2_correction <- rbind.data.frame(r1, met_mnx2_correction)
    met_mnx2_correction$compostion <- str_replace_all(met_mnx2_correction$compostion, "m", "product")
  }
  
  
  # step3
  rxn_no_biocyc_seed$reverse_correction <- "need_manual_check"
  met_mnx_no_correction <- met_ini[met_ini$ID %in% rxn_no_biocyc_seed$rxnID, ]
  met_mnx_no_correction$compostion_new <- NA
  
  
  # combine all the information together
  met_mnx_new <- rbind.data.frame(met_mnx1_correction, met_mnx2_correction, met_mnx_no_correction)
  rxn_with_only_seed <- subset(rxn_with_only_seed, select = -c(seed_equation_new))
  union_rxn_mnx_new <- rbind.data.frame(rxn_with_biocyc, rxn_with_only_seed, rxn_no_biocyc_seed)
  
  result <- list()
  result[["met"]] <- met_mnx_new
  result[["rxn"]] <- union_rxn_mnx_new
  
  return(result)
}








getRxnInfFromKEGG <- function(rxn_frame, kegg_rxnID_list) {
  # This function is used to get the detailed reaction annotation from KEGG database based on the KEGG ID of the reaction
  # Input
  # rxn_frame: A dataframe contains a column of KEGG ID
  # kegg_rxnID_list: A vector contains the above KEGG ID in the rxn_frame
  ec_rxn_kegg <- read.table("data/EC_rxn_mapping_kegg.txt", sep = "\t", stringsAsFactors = FALSE)
  ec_rxn_kegg$V2 <- str_replace_all(ec_rxn_kegg$V2, "rn:", "kegg:")
  rxn_formula_kegg <- read.delim2("data/reaction_kegg summary_2019_8.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  f1 <- rxn_formula_kegg$V2
  rxn_name <- vector()
  rxn_formula0 <- vector()
  for (i in 1:length(f1)) {
    s0 <- f1[i]
    s0 <- unlist(str_split(s0, "; "))
    n1 <- s0[!str_detect(s0, "=")]
    n2 <- s0[str_detect(s0, "=")]
    n1 <- paste(n1, collapse = ";")
    n2 <- paste(n2, collapse = ";")
    rxn_name <- c(rxn_name, n1)
    rxn_formula0 <- c(rxn_formula0, n2)
  }
  rxn_formula_kegg$formula <- rxn_formula0
  rxn_formula_kegg$rxn_name <- rxn_name
  colnames(rxn_formula_kegg) <- c("rxnID", "inf", "reaction", "rxn_name")
  rxn_formula_kegg$rxnID <- str_replace_all(rxn_formula_kegg$rxnID, "rn:", "kegg:")
  rxn_formula_kegg$reaction <- str_trim(rxn_formula_kegg$reaction, side = "both")
  
  rxn_frame$rxn_name <- getSingleReactionFormula(rxn_formula_kegg$rxn_name, rxn_formula_kegg$rxnID, kegg_rxnID_list)
  rxn_frame$formula_kegg <- getSingleReactionFormula(rxn_formula_kegg$reaction, rxn_formula_kegg$rxnID, kegg_rxnID_list)
  rxn_frame$EC_kegg <- getSingleReactionFormula(ec_rxn_kegg$V1, ec_rxn_kegg$V2, kegg_rxnID_list)
  
  return(rxn_frame)
}



getRheaRxnInfFromKEGGID <- function(rxn_frame=union_rxn_kegg, kegg_rxnID_list=union_rxn_kegg$rxnID) {
  # This function is used to get the detailed reaction annotation from Rhea database based on the KEGG ID of the reaction
  # Input
  # rxn_frame: A dataframe contains a column of KEGG ID
  # kegg_rxnID_list: A vector contains the above KEGG ID in the rxn_frame
  
  # note: in rhea database, a master rhea id could have three reaction rheaid, thus we not use rhea information to estimate the reaction
  # direction
  if (str_detect(kegg_rxnID_list[1], "kegg")) {
    kegg_rxnID_list <- str_replace_all(kegg_rxnID_list, "kegg:", "")
  }
  rhea2kegg <- read_tsv("data/rhea2kegg_reaction_2019_8.tsv")
  rhea2kegg$RHEA_ID <- paste("RHEA:", rhea2kegg$RHEA_ID, sep = "")
  rhea2kegg$MASTER_ID <- paste("RHEA:", rhea2kegg$MASTER_ID, sep = "")
  rxn_frame$RHEA_ID <- getSingleReactionFormula(rhea2kegg$MASTER_ID, rhea2kegg$ID, kegg_rxnID_list)
  
  # input the latest reaction annotation from rhea
  #rhea_reaction_summary <- read_csv("data/rhea reaction summary.csv")
  #rhea_reaction_summary$masterID <- paste("RHEA:", rhea_reaction_summary$masterID, sep = "")
  rhea_reaction_summary <- read_excel("data/rhea_reaction_summary.xlsx")
  rxn_frame$formula_rhea <- getMultipleReactionFormula(rhea_reaction_summary$formula, rhea_reaction_summary$rheaID, rxn_frame$RHEA_ID)
  rxn_frame$equation_rhea <- getMultipleReactionFormula(rhea_reaction_summary$equation, rhea_reaction_summary$rheaID, rxn_frame$RHEA_ID)
  return(rxn_frame)
}



getBiocycRxnInfFromBiocycID <- function(rxn_frame, biocyc_rxnID_list) {
  # This function is used to get the detailed reaction annotation from biocyc database based on the biocyc ID of the reaction
  # Input
  # rxn_frame: A dataframe contains a column of biocyc ID
  # biocyc_rxnID_list: A vector contains the above biocyc ID in the rxn_frame
  
  if (str_detect(biocyc_rxnID_list[1], "metacyc")) {
    biocyc_rxnID_list <- str_replace_all(biocyc_rxnID_list, "metacyc:", "")
  }
  biocyc_reaction <- read_excel("data/metcyc_rxn_RAVEN.xlsx")
  rxn_frame$rxn_name <- getSingleReactionFormula(biocyc_reaction$name, biocyc_reaction$Reactions, biocyc_rxnID_list)
  rxn_frame$reverse_biocyc <- getSingleReactionFormula(biocyc_reaction$Reaction.Direction, biocyc_reaction$Reactions, biocyc_rxnID_list)
  rxn_frame$formula_biocyc <- getSingleReactionFormula(biocyc_reaction$formula, biocyc_reaction$Reactions, biocyc_rxnID_list)
  rxn_frame$EC_biocyc <- getSingleReactionFormula(biocyc_reaction$EC.Number, biocyc_reaction$Reactions, biocyc_rxnID_list)
  
  return(rxn_frame)
}



getRheaRxnInfFromBiocycID <- function(rxn_frame, biocyc_rxnID_list) {
  # This function is used to get the detailed reaction annotation from biocyc database based on the biocyc ID of the reaction
  # Input
  # rxn_frame: A dataframe contains a column of biocyc_rxnID
  # biocyc_rxnID_list: A vector contains the above biocyc ID in the rxn_frame
  
  if (str_detect(biocyc_rxnID_list[1], "metacyc")) {
    biocyc_rxnID_list <- str_replace_all(biocyc_rxnID_list, "metacyc:", "")
  }
  
  rhea2metacyc <- read_tsv("data/rhea2metacyc.tsv")
  rxn_frame$rheaID <- getSingleReactionFormula(rhea2metacyc$RHEA_ID, rhea2metacyc$ID, biocyc_rxnID_list)
  rxn_frame$DIRECTION_rhea <- getSingleReactionFormula(rhea2metacyc$DIRECTION, rhea2metacyc$ID, biocyc_rxnID_list)
  
  return(rxn_frame)
}



metacycRxnIDmapping <- function() {
  # This fucntion is used to summarize  the reaction ID mapping between metacyc and kegg or rhea
  # The original ID mapping is downloaded from Metacyc database on August 25, 2019
  rxn_MetaCyc <- read_excel("data/All_reactions_of_MetaCyc.xlsx")
  # establish the ID mapping between metacyc, keggID and rheaID
  rxn_MetaCyc$KEGGID <- str_extract_all(rxn_MetaCyc$KEGG, "rn:[:graph:]+")
  
  rxn_MetaCyc$RHEAID <- str_extract_all(rxn_MetaCyc$RHEA, "id=[:graph:]+")
  rxn_MetaCyc0 <- rxn_MetaCyc[, c("Reaction", "KEGGID", "RHEAID")]
  
  rxn_MetaCyc0$KEGGID <- str_replace_all(rxn_MetaCyc0$KEGGID, "</a>", "") %>%
    str_replace_all(., "'", "") %>%
    str_replace_all(., "rn:", "") %>%
    str_replace_all(., ">", ";")
  rxn_MetaCyc0$RHEAID <- str_replace_all(rxn_MetaCyc0$RHEAID, "</a>", "") %>%
    str_replace_all(., "'", "") %>%
    str_replace_all(., "id=", "") %>%
    str_replace_all(., ">", ";")
  
  
  for (i in 1:nrow(rxn_MetaCyc0)) {
    print(i)
    k1 <- unique(unlist(str_extract_all(rxn_MetaCyc0$KEGGID[i], "R[:digit:]+")))
    k1 <- paste0(k1, collapse = ";")
    rxn_MetaCyc0$KEGGID0[i] <- k1
    
    r1 <- unique(unlist(str_extract_all(rxn_MetaCyc0$RHEAID[i], "[:digit:]+")))
    r1 <- paste0(r1, collapse = ";")
    rxn_MetaCyc0$RHEAID0[i] <- r1
  }
  
  # metacyc-keggid
  meta_kegg <- splitAndCombine(rxn_MetaCyc0$KEGGID0, rxn_MetaCyc0$Reaction, sep0 = ";")
  colnames(meta_kegg) <- c("KEGGID", "metacycID")
  meta_kegg$KEGGID <- str_trim(meta_kegg$KEGGID, side = "both")
  meta_kegg$metacycID <- str_trim(meta_kegg$metacycID, side = "both")
  meta_kegg <- meta_kegg[!(meta_kegg$KEGGID == "NA"), ]
  meta_kegg$KEGGID <- paste("kegg:", meta_kegg$KEGGID, sep = "")
  meta_kegg$metacycID <- paste("metacyc:", meta_kegg$metacycID, sep = "")
  
  meta_rhea <- splitAndCombine(rxn_MetaCyc0$RHEAID0, rxn_MetaCyc0$Reaction, sep0 = ";")
  colnames(meta_rhea) <- c("rheaID", "metacycID")
  meta_rhea$rheaID <- str_trim(meta_rhea$rheaID, side = "both")
  meta_rhea$metacycID <- str_trim(meta_rhea$metacycID, side = "both")
  meta_rhea <- meta_rhea[!(meta_rhea$rheaID == "NA"), ]
  meta_rhea$rheaID <- paste("rhea:", meta_rhea$rheaID, sep = "")
  meta_rhea$metacycID <- paste("metacyc:", meta_rhea$metacycID, sep = "")
  
  
  # rename the column and combine them together
  colnames(meta_kegg) <- c("kegg_or_rhea", "metacycID")
  colnames(meta_rhea) <- c("kegg_or_rhea", "metacycID")
  
  xref_biocyc <- rbind.data.frame(meta_kegg, meta_rhea)
  
  return(xref_biocyc)
}















metStandard.mnx <- function(met_frame, met_MNXID_list) {
  # This function is used to standard the metabolite based MNXID from MNX database
  # Input
  # met_frame: A dataframe contain a column -'MetID'
  # met_MNXID_list: A vector of met_MNXID_list
  
  # rxn_final0 <- select(union_rxn_mnx, rxnID, Equation_MNX)
  # colnames(rxn_final0) <- c("ID0", "Equation")
  
  #rxn_final0 <- data.frame(ID0=MNXID_list, Equation=MNX_equation_list)
  # remove the compartment information from metnet database
  #rxn_met <- splitRxnToMetabolite(rxn_final0, sep0 = " = ")
  #met_frame=rxn_biocyc1
  #met_MNXID_list=rxn_biocyc1$MetID
  
  rxn_met <- met_frame
  met_sum <- data.frame(MetID = unique(met_MNXID_list), stringsAsFactors = FALSE)
  chem_prop <- read_tsv("data/chem_prop_metanetx.tsv")
  met_sum$Formula <- getMultipleReactionFormula(chem_prop$Formula, chem_prop$MNX_ID, met_sum$MetID)
  met_sum$Charge <- getMultipleReactionFormula(chem_prop$Charge, chem_prop$MNX_ID, met_sum$MetID)
  met_sum$Description <- getMultipleReactionFormula(chem_prop$Description, chem_prop$MNX_ID, met_sum$MetID)
  
  # find the keggid and chebiID based on metnetID
  chem_xref <- read_tsv("data/chem_xref.tsv")
  chem_xref_kegg <- chem_xref[str_detect(chem_xref$XREF, "kegg:"), ]
  chem_xref_chebi <- chem_xref[str_detect(chem_xref$XREF, "chebi:"), ]
  met_sum$chebi_multiple <- getMultipleReactionFormula(chem_xref_chebi$XREF, chem_xref_chebi$MNX_ID, met_sum$MetID)
  met_sum$kegg_multiple <- getMultipleReactionFormula(chem_xref_kegg$XREF, chem_xref_kegg$MNX_ID, met_sum$MetID)
  
  chem_prop_chebi <- chem_prop[str_detect(chem_prop$Source, "chebi"), ]
  chem_prop_kegg <- chem_prop[str_detect(chem_prop$Source, "kegg"), ]
  
  met_sum$chebi <- getMultipleReactionFormula(chem_prop_chebi$Source, chem_prop_chebi$MNX_ID, met_sum$MetID)
  met_sum$kegg <- getMultipleReactionFormula(chem_prop_kegg$Source, chem_prop_kegg$MNX_ID, met_sum$MetID)
  # find the keggid based on chebiID
  chebi <- read_excel("data/chebi_compound_comprehensive.xlsx", sheet = "Sheet1")
  chebi$chebiID <- str_replace_all(chebi$chebiID, "CHEBI", "chebi")
  met_sum$kegg2 <- getMultipleReactionFormula(chebi$keggID, chebi$chebiID, met_sum$chebi)
  met_sum$kegg2[is.na(met_sum$kegg2)] <- "NA"
  
  for (i in seq_along(met_sum$kegg2)) {
    if (met_sum$kegg2[i] == "NA") {
      met_sum$kegg2[i] <- met_sum$kegg_multiple[i]
    } else {
      met_sum$kegg2[i] <- met_sum$kegg2[i]
    }
  }
  met_sum$kegg2 <- str_replace_all(met_sum$kegg2, "kegg:", "")
  met_sum <- subset(met_sum, select = -c(kegg, kegg_multiple, chebi_multiple))
  rxn_met0 <- left_join(rxn_met, met_sum, by.x = "MetID", by.y = "MetID")
  return(rxn_met0)
}




metStandard.biocyc <- function(rxn_inf = union_rxn_biocyc) {
  # This function is used to get the metaboite annotation for reaction from biocyc database
  # Input
  # rxn_inf: A dataframe should contains the column -'rxnID', 'formula_biocyc'
  
  rxn_biocyc <- rxn_inf[, c("rxnID", "formula_biocyc")]
  colnames(rxn_biocyc) <- c("ID0", "Equation")
  rxn_biocyc$Equation <- str_replace_all(rxn_biocyc$Equation, "<=>", "=>")
  rxn_biocyc1 <- splitRxnToMetabolite(reationFrame = rxn_biocyc, sep0 = "=>", source0 = "biocyc")
  # find MNXID of metabolite based on biocycid
  chem_xref <- read_tsv("data/chem_xref.tsv")
  chem_xref_metacyc <- chem_xref[str_detect(chem_xref$XREF,'metacyc'),]
  chem_xref_metacyc$XREF <- str_replace_all(chem_xref_metacyc$XREF, "metacyc:", "")
  rxn_biocyc1$MetID <- getMultipleReactionFormula(chem_xref_metacyc$MNX_ID, chem_xref_metacyc$XREF, rxn_biocyc1$MetID)
  rxn_biocyc1 <- metStandard.mnx(met_frame=rxn_biocyc1, met_MNXID_list=rxn_biocyc1$MetID)
  return(rxn_biocyc1)
}






# An example to plot the vnn graph
library(VennDiagram)
#an example
#summary
#BIGG <- unique(newBigg$panID)
#RAVEN <- unique(newRaven0$panID)
#Uniprot <- unique(newUniprot$panID)
#TCDB <- unique(newTCDB$panID)
#plot the graph
#venn.diagram(x= list(BIGG = BIGG, RAVEN = RAVEN, Uniprot = Uniprot, TCDB = TCDB), 
#             filename = "new gene in main database.png", height = 1000, width = 1000,resolution =300, imagetype="png", col="transparent",
#             fill=c("cornflowerblue","green","darkgreen","darkorchid1"),alpha = 0.50, cex=0.45, cat.cex=0.45)


plotDensityGraph <- function (data_frame, para, xlab_name, ylab_name='Density', title=''){
  # This function is used to plot the density plot for a single parameter
  # Input
  # data_frame: A dataframe
  # para: the column name of the parameter
  
  j = para
  ggplot(data_frame, aes_string(j)) +
    geom_density(fill="lightblue", alpha=1) +
    xlab(xlab_name) + 
    ylab(ylab_name) +
    theme(axis.text=element_text(size=20, family="Arial"),
          axis.title=element_text(size=24, family="Arial"),
          legend.text = element_text(size=20, family="Arial")) +
    ggtitle(title) +
    theme(panel.background = element_rect(fill = "white", color="black", size = 1)) 
  #+ggsave(out <- paste('result/','Metabolic gene number distribution for strains specific model from biocyc','.eps', sep = ""), width=5, height=5, dpi=300)
}


plotXYdotGraph <- function (data_frame, paraX, paraY, xlab_name, ylab_name, title=''){
  # This function is used to plot the dot plot between two parameters
  # Input
  # data_frame: A dataframe
  # paraX: the column name of the parameter
  # paraY: the column name of another parameter
  
  ggplot(data_frame, aes_string(x=paraX, y=paraY)) +
    geom_point() +
    stat_smooth(method = "loess",
                col = "#C42126",
                se = FALSE,
                size = 1) +
    xlab(xlab_name) + 
    ylab(ylab_name) +
    theme(axis.text=element_text(size=20, family="Arial"),
          axis.title=element_text(size=24, family="Arial"),
          legend.text = element_text(size=20, family="Arial")) +
    ggtitle(title) +
    theme(panel.background = element_rect(fill = "white", color="black", size = 1),
          plot.margin = margin(1, 1, 1, 1, "cm")) 
   
  #+ggsave(out <- paste('result/','Metabolic gene number distribution for strains specific model from biocyc','.eps', sep = ""), width=5, height=5, dpi=300)
}


plotXYdotGraph2 <- function (data_frame, paraX, paraY, xlab_name, ylab_name, title=''){
  # This function is used to plot the dot plot between two parameters
  # A bar column on the x and y axis will be ueed to do the statistical analysis
  # Input
  # data_frame: A dataframe
  # paraX: the column name of the parameter
  # paraY: the column name of another parameter
  
  p <- ggplot(data_frame, aes_string(x=paraX, y=paraY)) +
    geom_point(color="#69b3a2", alpha=0.8) +
    xlab(xlab_name) + 
    ylab(ylab_name) +
    theme(axis.text=element_text(size=20, family="Arial"),
          axis.title=element_text(size=24, family="Arial"),
          legend.text = element_text(size=20, family="Arial")) +
    ggtitle(title) +
    theme(panel.background = element_rect(fill = "white", color="black", size = 1),
          plot.margin = margin(1, 1, 1, 1, "cm")) 
  ggExtra::ggMarginal(p, type = "histogram", color="grey")
  #+ggsave(out <- paste('result/','Metabolic gene number distribution for strains specific model from biocyc','.eps', sep = ""), width=5, height=5, dpi=300)
}
