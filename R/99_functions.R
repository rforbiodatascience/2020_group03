# Define project functions
# ------------------------------------------------------------------------------ 
encode_peptide = function(x, matrix){
  matrix <- paste("/Users/laurasansc/github/2020_group03/data/_raw/",matrix, ".txt", sep="")
  m <- read_table2(file = matrix)
  X_enc = x %>%
    str_split('') %>%
    lapply(function(x_i){
      m[x_i,] %>%
        as.matrix %>%
        t %>%
        matrix(nrow = 1, byrow = TRUE) %>%
        return
    })
  X_enc = do.call(rbind, X_enc)
  rownames(X_enc) = x
  return(X_enc)
}