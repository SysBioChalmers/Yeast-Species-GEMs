# Add pathway information for the new rxn
# One rxn could have multiple pathways
# Revised by Hongzhong 2020-03-12

# load library
source('function_general.R')
new_rxn <- read_excel("result/new rxn information from MNX database.xlsx")

kegg_rxn_pathway <- read.table("data/kegg/reactionID_pathway.txt", stringsAsFactors = FALSE)
kegg_rxn_pathway$V1 <- str_replace_all(kegg_rxn_pathway$V1, "rn:", "kegg:")
kegg_rxn_pathway <- kegg_rxn_pathway[str_detect(kegg_rxn_pathway$V2,"map"),]
pathway_list_kegg <- read_excel("data/kegg/pathway_list_kegg.xlsx", col_names = FALSE)
colnames(pathway_list_kegg) <- c("id", "name")
kegg_rxn_pathway$name <-  getSingleReactionFormula(pathway_list_kegg$name,pathway_list_kegg$id,kegg_rxn_pathway$V2)
new_rxn$kegg_pathway_name <- getMultipleReactionFormula(kegg_rxn_pathway$name,kegg_rxn_pathway$V1,new_rxn$keggID)



write.table(new_rxn, "result/new_rxn_with_pathway.txt", row.names = FALSE, sep = "\t")

