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

# Simple encoding function for converting a given peptide 'x' to a numerical
# representation using the 'm' matrix
encode_peptide = function(x, m){
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
#
# Usage:
# make_sequence(Mutated_residue, Mutation_position, Mutation, uniprot_ID)
#
# Returns:
# modified sequence
make_sequence = function(Mutated_residue, Mutation_position, Mutation, uniprot_ID)
{
  WT_sequence <- as_tibble(GetSequences(uniprot_ID)) %>%
    pull(Sequence)
  sequence <- WT_sequence
  str_sub(sequence, Mutation_position, Mutation_position) <- Mutation
  return(sequence)
}

# Function to convert a mutation description into an mutated residue
#
# Usage:
# get_mutated_residue(variant)
#
# Returns:
# mutated residue
get_mutated_residue <- function(variant) {
  mutated_residue = str_extract(variant, "[A-Z]")
  return(mutated_residue)
}

# Function to convert a mutation description into an mutattion position
#
# Usage:
# get_mutation_position(variant)
#
# Returns:
# mutation position
get_mutation_position <- function(variant) {
  mutation_position = str_extract(variant, "[0-9]+")
  return(mutation_position)
}

# Function to convert a mutation description into an mutation
#
# Usage:
# get_mutation(variant)
#
# Returns:
# mutation
get_mutation <- function(variant) {
  mutation = str_extract(variant, "[A-Z,*]$")
  return(mutation)
}

# Function to convert variant to sequence
#
# Usage:
# convert_variant_to_sequence(tibble, uniprot_id)
# Returns:
# tibble augmented with sequence
convert_variant_to_sequence <- function(df, protein)  {
  clean_df <- df %>%
  mutate(
    mutated_residue = get_mutated_residue(variant), 
    mutation_position = get_mutation_position(variant),
    mutation = get_mutation(variant),
    sequence = make_sequence(mutated_residue, mutation_position, mutation, protein))
  return(clean_df)
}