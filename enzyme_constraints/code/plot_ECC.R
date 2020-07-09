library(pheatmap)
#Get heatmap with growth control coefficients for different enzymes across strains
if (exists("RStudio.Version")){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(getSrcDirectory()[1])
}
files <- list.files(path = "../results/gCC")
i <- 1
protNumber <- 3
topRxns    <- c()
organisms  <- c()
#Get unique top enzymes
for (file in files){
  if (grepl('.txt',file)){
    name <- gsub('_limGrowth.txt','',file)
    fileName  <- paste("../results/gCC/",file,sep='')
    dataset   <- read.delim(fileName, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
    organisms <- c(organisms,name)
    dataset   <- dataset[1:protNumber,]
    rxnNames  <- dataset[,6]
    ECC       <- dataset[,5]
    topRxns   <- c(topRxns,rxnNames)
    i <- i + 1
  }
}
topRxns    <- unique(topRxns)
shortNames <- c('HMG2','ARO2','TSA1','ARG7','ACC1','ILV5','GPH1','ERG13','ADE6','IMD2') #All unconstrained
shortNames <- c('HMG2','ARO2','SAH1','FBA1','TPS1','TPS2','ACO1','ADE12','TSA1','ILV5','TDH3','MDH1','PSD1') #All unconstrained

#shortNames <- c('HMG2','ARO2','TSA1','ARG7','ACC1','MDH1','GPH1','ERG13','IMD2','TDH3','ADE6','IDH2')

keyCodes   <- c('Dbr','Esi','Kla','Kmx','Kpa','Lfe','Lth','Nca','Seu','Spo','Tbl','Tph','Zro','Kdo')
#Get the ECC of such enzymes
i <- 1
ECCs <- matrix(0, nrow = length(topRxns), ncol = length(organisms))
for (file in files){
  if (grepl('.txt',file)){
    #Open each file
    name <- gsub('_limGrowth.txt','',file)
    fileName <- paste("../results/gCC/",file,sep='')
    dataset  <- read.delim(fileName, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
    rxnNames <- dataset[,6]
    ECC      <- log10(dataset[,4])
    #numericla values
    #ECC      <- (dataset[,5])
    #
    j <- 1
    #look for each numerical value according to the order of top reactions
    for (rxn in topRxns){
      index <- which(rxnNames==rxn)
      #print(rxnNames[index])
      if (length(index)>0){
        #In case that there are several reactions with the same name, then take the maximum numerical value
        ECCs[j,i] <- max(ECC[index])
      } 
      else{
        ECCs[j,i] <- 0
      }
      j <- j+1
    }
  }
  i <- i+1
}
ECCs <- as.data.frame(ECCs)
colnames(ECCs) <- keyCodes 
ECCs <- t(ECCs)
ECC_sums <- colSums(ECCs)
orderIndx <- order(ECC_sums,decreasing = TRUE)
orderIndx <- orderIndx[1:10]
ECCs <- ECCs[,orderIndx]
topRxns <- topRxns[orderIndx]
colnames(ECCs) <- c('TRP5','HMG1','ERG13','ILV5','IMD2','ADE6','ARG7','GPH1','PSD1','TPH3')#,'ADE4','PSD1','FBA1','ILV6')
colnames(ECCs) <- c('TRP5','HMG1','ERG13','ILV5','IMD2','ARG7','GPH1','TDH3','HMG2','ACC1')#,'ADE4','PSD1','FBA1','ILV6')

#plot heatmap
fileName  <- '../results/figures/heatMap_topECC_kcats.png'
png(fileName,width=600, height=400)
pheatmap(ECCs,cluster_cols = F,cluster_rows = T, show_rownames = TRUE,scale='none',fontsize = 12)
dev.off()