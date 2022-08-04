jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}
#a <- c(1,5,8,10,22,14,15,16,2,7)
#b <- c(10,12,13,2,7,9,2,7,23,15)
#jaccard(a,b)
# for python version, please find it in https://stackoverflow.com/questions/37003272/how-to-compute-jaccard-similarity-from-a-pandas-dataframe



# test input two models
#---------------------------------------------------
# input the model from kegg
strain <- list.files('../ComplementaryData/draft_yeast_GEMs/strain_specific_model_from_RAVEN_kegg')


# get the combination
strain_combine <-  as.data.frame(combn(strain, 2), stringsAsFactors = FALSE)
result_df <- data.frame()
for (i in 1:ncol(strain_combine)){
  print(i)
  strain_combine[,c(i)]
  ss <-  strain_combine[,c(i)]
  s1 <- ss[1]
  s2 <- ss[2]
  file1 <- paste('../ComplementaryData/draft_yeast_GEMs/strain_specific_model_from_RAVEN_kegg/',s1, '/excelRxns.txt', sep = "")
  file2 <- paste('../ComplementaryData/draft_yeast_GEMs/strain_specific_model_from_RAVEN_kegg/',s2, '/excelRxns.txt', sep = "")
  rxn1 <- read_table2(file1)
  rxn2 <- read_table2(file2)
  similirity <- jaccard(rxn1[["#"]],rxn2[["#"]])
  
  df1 = data.frame(s1,s2,similirity)
  
  # adding names to the row values
  names(df1)=c("s1","s2","similirity") 
  
  # passing the original data frame and 
  # new data frame into the rbind() function 
  result_df=rbind(result_df,df1) 
}

write.table(result_df, "result/model_similirity_from_kegg.txt", row.names = FALSE, sep = "\t")






# input the model from biocyc
strain <- list.files('../ComplementaryData/draft_yeast_GEMs/strain_specific_model_from_RAVEN_biocyc_55_110')


# get the combination
strain_combine <-  as.data.frame(combn(strain, 2), stringsAsFactors = FALSE)
result_df <- data.frame()
for (i in 1:ncol(strain_combine)){
  print(i)
  strain_combine[,c(i)]
  ss <-  strain_combine[,c(i)]
  s1 <- ss[1]
  s2 <- ss[2]
  file1 <- paste('../ComplementaryData/draft_yeast_GEMs/strain_specific_model_from_RAVEN_biocyc_55_110/',s1, '/excelRxns.txt', sep = "")
  file2 <- paste('../ComplementaryData/draft_yeast_GEMs/strain_specific_model_from_RAVEN_biocyc_55_110/',s2, '/excelRxns.txt', sep = "")
  rxn1 <- read_table2(file1)
  rxn2 <- read_table2(file2)
  similirity <- jaccard(rxn1[["#"]],rxn2[["#"]])
  
  df1 = data.frame(s1,s2,similirity)
  
  # adding names to the row values
  names(df1)=c("s1","s2","similirity") 
  
  # passing the original data frame and 
  # new data frame into the rbind() function 
  result_df=rbind(result_df,df1) 
}

write.table(result_df, "result/model_similirity_from_biocyc.txt", row.names = FALSE, sep = "\t")

