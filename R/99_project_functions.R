# Define project functions
# ------------------------------------------------------------------------------

# Simple function to load a given amino acid scale
load_scale = function(enc)
{
  if (enc == 'z5')
  {
    encoder = read.csv(file = './data/z_scales.csv', sep = ";", dec = ",", header=TRUE, row.names=1)
  }
  if (enc == 'bl62')
  {
    encoder = read.table(file = './data/BLOSUM62.txt')
  }
  return (encoder)
}

#' Simple encoding function for converting a given peptide 'x' to a numerical
#' representation using the 'm' matrix
encode_peptide = function(x, m)
{
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