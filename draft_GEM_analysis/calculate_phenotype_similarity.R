library(readxl)

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}
#a <- c(1,5,8,10,22,14,15,16,2,7)
#b <- c(10,12,13,2,7,9,2,7,23,15)
#jaccard(a,b)
# for python version, please find it in https://stackoverflow.com/questions/37003272/how-to-compute-jaccard-similarity-from-a-pandas-dataframe

yeast_traits <- read_excel("../ComplementaryData/phenotype/yeast_traits.xlsx", 
                           sheet = "Substrate utilization")



yeast_traits <- yeast_traits[2:60, ]
colnames(yeast_traits) <- yeast_traits[1,]
yeast_traits <- yeast_traits[2:59, ]


# input phenotype information
strain <- colnames(yeast_traits)[5:333]
# get the combination
strain_combine <-  as.data.frame(combn(strain, 2), stringsAsFactors = FALSE)
result_df <- data.frame()
for (i in 1:ncol(strain_combine)){
  print(i)
  strain_combine[,c(i)]
  ss <-  strain_combine[,c(i)]
  s1 <- ss[1]
  s2 <- ss[2]
  r1 <- yeast_traits[[s1]]
  r2 <- yeast_traits[[s2]]
  
  index1 <- which(r1=="n")
  index2 <- which(r2=="n")
  index_remove <- union(index1, index2)
  
  r10 <- r1[-index_remove]
  r20 <- r2[-index_remove]
  
  similirity <- length(which(r10 == r20))/length(r10) # similar to jaccard function
  
  df1 = data.frame(s1,s2,similirity)
  
  # adding names to the row values
  names(df1)=c("s1","s2","similirity") 
  
  # passing the original data frame and 
  # new data frame into the rbind() function 
  result_df=rbind(result_df,df1) 
}

write.table(result_df, "result/trait_similirity.txt", row.names = FALSE, sep = "\t")

