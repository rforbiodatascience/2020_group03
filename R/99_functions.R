# Define project functions
# ------------------------------------------------------------------------------ 
#' Encode protein or peptide sequences
#'
#' @description
#' `encode_peptide` returns a dataframe encoding the sequence in the method `m' specified.
#' All the peptides have to be the same length.
#'
#' @param x peptide, string
#' @param matrix method, string
#' 
#' @return if all the parameters are logical will return a dataframe with the row.names = peptide and each column the encoding. The number of columns is len(peptide) x 20 amino acids. 
#' 
#' @example
#' encode_peptide('ATRALPLTW', "blosum62")
#' 
encode_peptide = function(x, m){
  mat <- paste("./data/_raw/",m, ".txt", sep="")
  encoder <- read.table(file = mat, row.names = 1, header =TRUE)
  X_enc = x %>%
    str_split('') %>%
    lapply(function(x_i){
      encoder[x_i,] %>%
        as.matrix %>%
        t %>%
        matrix(nrow = 1, byrow = TRUE) %>%
        return
    })
  X_enc = do.call(rbind, X_enc)
  rownames(X_enc) = x
  return(X_enc)
}

