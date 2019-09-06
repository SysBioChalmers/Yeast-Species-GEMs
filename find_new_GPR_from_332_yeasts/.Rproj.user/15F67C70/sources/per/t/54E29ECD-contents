library(readxl)
library(stringr)
library(tidyverse)
library(hongR)


# RAVEM-biocyc
# This is result using the defalut paramters, that is pidenty = 45%, bitscore = 100
pan_biocyc <- read.delim2('data/panYeast2_biocyc/excelRxns.txt', header =TRUE, sep = "\t", stringsAsFactors = FALSE)
sce_biocyc <- read.delim2('data/s288c_biocyc/excelRxns.txt', header =TRUE, sep = "\t", stringsAsFactors = FALSE)
newRxn_biocyc <- setdiff(pan_biocyc$ID, sce_biocyc$ID)
newRxn_biocyc0 <- pan_biocyc[pan_biocyc$ID %in% newRxn_biocyc,]
which(duplicated(newRxn_biocyc0$ID))
# when the pidenty = 55%, bitscore = 110
pan_biocyc <- read.delim2('data/panYeast2_biocyc_55_110/excelRxns.txt', header =TRUE, sep = "\t", stringsAsFactors = FALSE)
newRxn_biocyc <- setdiff(pan_biocyc$ID, sce_biocyc$ID)
newRxn_biocyc0_55_110 <- pan_biocyc[pan_biocyc$ID %in% newRxn_biocyc,]




# RAVEN-KEGG
#pan_kegg <- scan('data/panYeast2_kegg/excelRxns.txt', sep = "\n", what = "complex")
#pan_kegg[1] <- str_replace_all(pan_kegg[1], "#", '')
#writeLines(pan_kegg, file("data/panYeast2_kegg/excelRxns"))
infile1 = 'data/panYeast2_kegg/excelRxns.txt'
infile2 = 'data/s288c_kegg/excelRxns.txt'

InputKeggModelRaven <- function(infile) {
  # this function is used to preprocess the kegg model obtained by RAVEN
  # input
  # infile: the directory

  input_kegg <- read.delim2(infile, row.names = NULL, sep = "\t", stringsAsFactors = FALSE)
  # remove the first column and rename all columns
  input_kegg <- input_kegg[-c(1)]
  old_columns <- colnames(input_kegg)
  new_columns <- old_columns[-c(1, 2)]
  new_columns <- c(c("ID"), new_columns, c("other"))
  colnames(input_kegg) <- new_columns
  return(input_kegg)
}


pan_kegg <- InputKeggModelRaven(infile1)
sce_kegg <- InputKeggModelRaven(infile2)
newRxn_kegg <- setdiff(pan_kegg$ID, sce_kegg$ID)
newRxn_kegg0 <- pan_kegg[pan_kegg$ID %in% newRxn_kegg,]
# give the mnxid for all the new reactions

write.table(newRxn_biocyc0, "result/newRxn_biocyc_RAVEN.txt", row.names = FALSE, sep = "\t")
write.table(newRxn_biocyc0_55_110, "result/newRxn_biocyc_RAVEN_55_110.txt", row.names = FALSE, sep = "\t")
write.table(newRxn_kegg0, "result/newRxn_kegg_RAVEN.txt", row.names = FALSE, sep = "\t")


