# Compare panYeast from different sources
# Revised by Hongzhong 2020-03-31

# load library
source('function_general.R')

# rxn comparison
# summarize the current panYeast reactions and choose the common reactions sets
rxn_aybraham <- parseReportedPanYeast()
# input the yeast species pan model of our group
yeast_species_pan_03_30 <- read_excel("data/yeast_species_pan_03_30.xls")
rxn_pan_03_30 <- yeast_species_pan_03_30$`MNX ID`
rxn_pan_03_30 <- rxn_pan_03_30[!is.na(rxn_pan_03_30)]

#plot the graph of rxn
venn.diagram(x= list(panYeast_33 = rxn_aybraham, panYeast343 = rxn_pan_03_30), 
             filename = "result/rxn comparison of different panYeast.png", height = 1000, width = 1000,resolution =300, imagetype="png", col="transparent",
             fill=c("blue","red"),alpha = 0.50, cex=0.45, cat.cex=0.45)



###################################################################################################
# metabolite comparasion
met_aybraham <- read_excel("data/aybraham.xlsx", sheet = "metabolites")
met_aybraham <- met_aybraham$`Metabolite KEGG ID`
met_aybraham <- unique(met_aybraham[!is.na(met_aybraham)])


met_yeast343 <- read_excel("data/yeast_species_pan_03_30.xls", sheet = "Metabolite List")
met_yeast343 <- met_yeast343$`KEGG ID`
met_yeast343 <- unique(met_yeast343[!is.na(met_yeast343)])
#plot the graph of rxn
venn.diagram(x= list(panYeast_33 = met_aybraham, panYeast343 = met_yeast343), 
             filename = "result/met comparison of different panYeast.png", height = 1000, width = 1000,resolution =300, imagetype="png", col="transparent",
             fill=c("blue","red"),alpha = 0.50, cex=0.45, cat.cex=0.45)

