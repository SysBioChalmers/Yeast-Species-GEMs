# this code is used to find the new GPR from panID based on EC number predicted by deepec
# https://bitbucket.org/account/user/kaistsystemsbiology/projects/DEEP
# summary: the ec number predicted from deepec seems not so good as we can found a lot of ec number from s288c can't found from deepec result

# Then we compared the predicted EC with annotation of EC for s288c from SGD and uniprot
# Revised by hongzhong 2019-8-12

#library(tidyverse)
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(readr)
library(hongR)
source('main_function.R')


annotation_s288c <- read_excel("data/all_gene_yeast with annotation from different database.xlsx")
ec_combine <- paste(annotation_s288c$ec_SGD, annotation_s288c$ec_uniprot, sep = ";")
ec_combine0 <- unlist(str_split(ec_combine, ";"))
ec_combine0  <- str_trim(ec_combine0, side = "both")
ec_combine0 <- unique(ec_combine0[ec_combine0 !="NA"]) # total 957 unique enzymes


ec_p1 <- read.delim2('data/EC_number_annotation_deepec/part1/DeepEC_Result.txt', stringsAsFactors = FALSE)
ec_p2 <- read.delim2('data/EC_number_annotation_deepec/part2/DeepEC_Result.txt', stringsAsFactors = FALSE)
ec_p3 <- read.delim2('data/EC_number_annotation_deepec/part3/DeepEC_Result.txt', stringsAsFactors = FALSE)
ec_p4 <- read.delim2('data/EC_number_annotation_deepec/part4/DeepEC_Result.txt', stringsAsFactors = FALSE)
ec_p5 <- read.delim2('data/EC_number_annotation_deepec/part5/DeepEC_Result.txt', stringsAsFactors = FALSE)
ec_p6 <- read.delim2('data/EC_number_annotation_deepec/part6/DeepEC_Result.txt', stringsAsFactors = FALSE)

ec_deepec <- rbind.data.frame(ec_p1,ec_p2,ec_p3,ec_p4,ec_p5,ec_p6)
ec_deepec$Predicted.EC.number <- str_replace_all(ec_deepec$Predicted.EC.number, "EC:", "")

newEC <- setdiff(ec_deepec$Predicted.EC.number, ec_combine0)
geneNewEC <- ec_deepec[ec_deepec$Predicted.EC.number %in% newEC, ]
# further remove the new EC from s288c
geneNewEC <- geneNewEC[!str_detect(geneNewEC$Query.ID, "Saccharomyces_cerevisiae"), ]

write.table(geneNewEC, "result/newEC_predicted_by_deep_ec_for_pan_genome.txt", row.names = FALSE, sep = "\t")



