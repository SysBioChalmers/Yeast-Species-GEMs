  library(pheatmap)
  library(viridis)
  library(matrixStats)
  
  #Get heatmap with growth control coefficients for different enzymes across strains
  if (exists("RStudio.Version")){
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  } else {
    setwd(getSrcDirectory()[1])
  }
  files <- list.files(path = "../results/gCC")
  i <- 1
  
  topRxns    <- c()
  organisms  <- c()
  proteins   <- c()
  allPanGenes <- c()
  allPanNames <- c()
  #Get unique top enzymes
  for (file in files){
    if (grepl('.txt',file)){
      name <- gsub('_limGrowth.txt','',file)
      fileName  <- paste("../results/gCC/",file,sep='')
      dataset   <- read.delim(fileName, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
      protNumber <- nrow(dataset)
      organisms <- c(organisms,name)
      dataset   <- dataset[1:protNumber,]
      rxnNames  <- dataset$rxnNames
      ECC       <- dataset$ECC
      prot      <- dataset$enzyme
      proteins  <- c(proteins,prot)
      topRxns   <- c(topRxns,rxnNames)
      allPanGenes <- c(allPanGenes,dataset$enz_panIDs)
      i <- i + 1
    }
  }
  commonPanIDs    <- unique(allPanGenes)
  #shortNames <- c('HMG2','ARO2','TSA1','ARG7','ACC1','ILV5','GPH1','ERG13','ADE6','IMD2') #All unconstrained
  #shortNames <- c('HMG2','ACO1','ARO2','GAPDH','TPS1','IMD2','ERG13','PSD1','MET17','ILV5','ADE6','HOM2','PGM1') #All unconstrained
  
  keyCodes   <- c('Dbr','Esi','Kla','Kmx','Ppa','Lfe','Lth','Ncs','Seu','Spo','Tbl','Tph','Zro','Kdo')
  #Get the ECC of such enzymes
  i <- 1
  ECCs <- matrix(0, nrow = length(commonPanIDs), ncol = length(organisms))
  for (file in files){
    if (grepl('.txt',file)){
      #Open each file
      name <- gsub('_limGrowth.txt','',file)
      fileName <- paste("../results/gCC/",file,sep='')
      dataset  <- read.delim(fileName, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
      panIDs <- dataset$enz_panIDs
      kcat <- log10(dataset$Kcat)
      #numericla values
      ECC <- dataset$ECC#dataset$ECC#log10(dataset$Kcat)#
      #
      j <- 1
      #look for each numerical value according to the order of top reactions
      for (panID in commonPanIDs){
        index <- which(panIDs==panID)
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
  org_sums <- rowSums(ECCs)
  org_sums <- data.frame(rownames(ECCs),org_sums)
  #org_sums$row_std = rowSds(as.matrix(df1[,c(2,3,4)]))
  orderIndx <- order(ECC_sums,decreasing = TRUE)
  orderIndx <- orderIndx[1:20]
  ECCs <- ECCs[,orderIndx]
  commonPanIDs <- commonPanIDs[orderIndx]
  colnames(ECCs) <- c('ERG13','TRP5','HMG1','TDH2','IMD4',
                      'ADE4','ILV5','FBA1','ARO1','ARO2',
                      'PSD1','ALA1','PFK1','SAH1','PYK2',
                      'OLE1','LYS4','HXK2','ACO2','YLR446W')
  #colnames(ECCs) <- c('TRP5','HMG1','ERG13','ILV5','IMD2','ARG7','GPH1','TDH3','HMG2','ACC1')#,'ADE4','PSD1','FBA1','ILV6')
  #colnames(ECCs) <- shortNames[orderIndx]

  #plot heatmap
  fileName  <- '../results/figures/heatMap_topECC.pdf'
  pdf(fileName,width=16, height=8)
  pheatmap(ECCs,cluster_cols = F,color=cividis(100),cluster_rows = T, show_rownames = TRUE,scale='none',fontsize = 12)
  dev.off()