# Define project functions
# ------------------------------------------------------------------------------

# Simple function to load a given amino acid scale
load_scale = function(enc)
{
  if (enc == 'z5')
  {
    encoder = read.csv(file = './_raw/z_scales.csv', sep = ";", dec = ",", header=TRUE, row.names=1)
  }
  if (enc == 'bl62')
  {
    encoder = read.table(file = './_raw/BLOSUM62.txt')
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

# Simple function to convert a mutation into a protein sequence based on a unprot id
make_sequence = function(Mutated_residue, Mutation_position, Mutation, uniprot_ID)
{
  WT_sequence <- as_tibble(GetSequences(uniprot_ID)) %>%
    pull(Sequence)
  sequence <- WT_sequence
  str_sub(sequence, Mutation_position, Mutation_position) <- Mutation
  return(sequence)
}
