library(pheatmap)
#Get heatmap with growth control coefficients for different enzymes across strains
if (exists("RStudio.Version")){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(getSrcDirectory()[1])
}
files <- list.files(path = "../results/gCC")
i <- 1
protNumber <- 2
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
shortNames <- c('ACO2','ARO2','TSA1','MDH1','PSD1','ARG7','ILV5','HXK1','GPH1','IMD2','VAS1','HOM2')
keyCodes   <- c('Dbr','Esi','Kla','Kmx','Kpa','Lfe','Lth','Nca','Seu','Spo','Tbl','Tph','Zro','Kdo')
#Get the ECC of such enzymes
i <- 1
ECCs <- matrix(0, nrow = length(topRxns), ncol = length(organisms))
for (file in files){
  if (grepl('.txt',file)){
    name <- gsub('_limGrowth.txt','',file)
    fileName <- paste("../results/gCC/",file,sep='')
    dataset  <- read.delim(fileName, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
    rxnNames <- dataset[,6]
    ECC      <- dataset[,5]
    #
    j <- 1
    for (rxn in topRxns){
      index <- which(rxnNames==rxn)
      if (length(index)>0){
        ECCs[j,i] <- ECC[index]
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
rownames(ECCs) <- shortNames 
ECCs <- t(ECCs)
#plot heatmap
#png(fileName,width=Width, height=Height)
pheatmap(ECCs,cluster_cols = F,cluster_rows = F, show_rownames = TRUE)
#dev.off()