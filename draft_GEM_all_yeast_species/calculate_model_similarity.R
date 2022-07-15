
# a functional example
a <- c(1,5,8,10,22,14,15,16,2,7)
b <- c(10,12,13,2,7,9,2,7,23,15)
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

jaccard(a,b)

# for python version, please find it in https://stackoverflow.com/questions/37003272/how-to-compute-jaccard-similarity-from-a-pandas-dataframe



# test input two models














