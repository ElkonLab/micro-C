
# rowCV -------------------------------------------------------------

# Input:
#   A numeric matrix.
# 
# Output:
#   An atomic vector containing the coefficient of variation for each row in the input matrix.


rowCV <- function(x){
  ## Takes a matrix and return a vector of coefficient of variations for each row
  require(genefilter)
  mat <- as.matrix(x, na.rm = T)
  stdv <- rowSds(mat, na.rm=T)
  mn <- rowMeans(mat, na.rm=T)
  cv <- stdv/mn
}

# write.tsv ---------------------------------------------------------------

# A wrapper for the function write.table() 

write.tsv <- function(x, file, header = T){
  write.table(x, file, append = F, quote = F, sep = '\t', row.names = F, col.names = header)
}

# write.tsv ---------------------------------------------------------------

specificity_matrix <- function(data){
  mat <- as.matrix(data)
  means <- rowMeans(mat, na.rm = T)
  fc_mat <- mat/means
  df <- data.frame(fc_mat)
  colnames(df) <- colnames(data)
  df
}