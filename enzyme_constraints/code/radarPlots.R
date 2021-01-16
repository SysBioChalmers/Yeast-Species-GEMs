library(reshape)
library(ggbiplot)
library(ggplot2)
library(fmsb)
repoPath <- '/Users/ivand/Documents/GitHub/Yeast-Species-GEMs'
setwd(repoPath)
#Load plotting functions
#source('scripts/analyzeModels/plotData.R')

#Load datasets
dataPath <- paste(repoPath,'/enzyme_constraints/',sep='')
fileName <- paste(dataPath,'results/ecModels_metrics_improved.txt',sep='')
dataset  <- read.delim(fileName, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
dataset$enzCoverage <- dataset$ec_nEnz/dataset$nGenes
dataset$enzymes <- dataset$ec_enzymes
dataset$bioYield <- dataset$gRate/(dataset$GUR*0.18)
dataset$GUR_Error <- abs(dataset$GUR-dataset$GUR_exp)/dataset$GUR_exp

dataset$bioYield_exp <- dataset$gRate/(dataset$GUR_exp*0.18)
dataset$bioError <- abs(dataset$bioYield-dataset$bioYield_exp)/dataset$bioYield_exp
dataset$EtExc[which(dataset$EtExc==0)] <- 1E-8
dataset$EtOH_exp[which(dataset$EtOH_exp==0)] <- 1E-8
dataset$EtOH_Error <- abs(dataset$EtExc-dataset$EtOH_exp)/dataset$EtOH_exp
dataset$gRate_error <- rep(0,length(dataset$EtOH_Error))

#features <- c('ec_nEnz','isoenzymes','promiscuous','complexes')
features <- c('bioYield','bioYield_exp')
features <- c('GUR','GUR_exp')
features <- c('EtExc','EtOH_exp')
#features <- c('gRate','gRate_exp')

maxLim <- 20
minLim <- 0
dataset$key  <- c('Dbr','Esi','Kla','Kmx','Ppa','Lfe','Lth','Ncs','Seu','Spo','Tbl','Tpf','Zro','Kdo')

newData <- c()
for (feat in features){
column <- which(colnames(dataset)==feat)
newData <- rbind(newData,t(dataset[,column]))
}
colnames(newData) <- dataset$key
rownames(newData) <- features
newData <- rbind(rep(maxLim,ncol(newData)) , rep(minLim,ncol(newData)) , newData)
newData <- as.data.frame(newData)
newData <- (as.data.frame(newData))

# Color vector
if (length(features)==4){
colors_border=c(rgb(0.,0.,0.,0.9), rgb(0.8,0.2,0.2,0.9) , rgb(0.2,0.9,0.6,0.9),rgb(0.3,0.4,0.8,0.9) )
colors_in=c( rgb(0.,0.,0.,0.3), rgb(0.8,0.2,0.2,0.3) , rgb(0.2,0.9,0.6,0.3),rgb(0.3,0.4,0.8,0.3) )

}
if (length(features)==3){ 
  colors_border=c(rgb(0.,0.,0.,0.9),rgb(0.8,0.2,0.2,0.9),rgb(0.3,0.4,0.8,0.9))
  colors_in=c( rgb(0.,0.,0.,0.3),rgb(0.8,0.2,0.2,0.3),rgb(0.3,0.4,0.8,0.3) )
}
if (length(features)==2){ 
  colors_border=c(rgb(0.,0.,0.,0.9),rgb(0.8,0.2,0.2,0.9))
  colors_in=c( rgb(0.,0.,0.,0.3),rgb(0.8,0.2,0.2,0.3))
}
if (length(features)==1){ 
colors_border=c(rgb(0.,0.,0.,0.9))
colors_in=c( rgb(0.,0.,0.,0.3))
}
#colors_border=c(rgb(0.,0.,0.,0.9),rgb(0.8,0.2,0.2,0.9))
#colors_in=c( rgb(0.,0.,0.,0.3),rgb(0.8,0.2,0.2,0.3))

# plot with default options:
radarchart( newData  , axistype=1 , 
            #custom polygon
            pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="black", caxislabels=seq(minLim,maxLim*100,(maxLim-minLim)/4),cglwd=0.8,
            #custom labels
            vlcex=1 
)

# Add a legend
legendStr <- c('bioError','EtOH_Error')
legendStr <- c('Predicted','Experimental')
#legendStr <- c('PredictedEt','ExperimentalEt','PredictedGlc','ExperimentalGlc')
#legendStr <- c('enzymes','isoenzymes','promiscuous','complexes')
legend(x=1, y=1.3, legend =legendStr, bty = "n", pch=20 , col=colors_in , text.col = "black", cex=1.2, pt.cex=3)
