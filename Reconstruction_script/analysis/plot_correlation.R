library("PerformanceAnalytics")

group_file <- read.table("/Users/feiranl/Documents/GitHub/Yeast-Species-GEMs/Reconstruction_script/analysis/table.txt",header = T,sep = ",") # file name should be 'Evolution_event_table.txt'

M<-cor(group_file)
M2 = M
corrplot.mixed(M2, lower = "ellipse", upper = "circle")
corrplot.mixed(M2, lower = "square", upper = "circle")
corrplot.mixed(M2, lower = "shade", upper = "circle")
corrplot.mixed(M2, tl.pos = "lt")
corrplot.mixed(M2, tl.pos = "lt", diag = "u")
corrplot.mixed(M2, tl.pos = "lt", diag = "l")
corrplot.mixed(M2, tl.pos = "n")


chart.Correlation(group_file, histogram=TRUE, pch=10)
