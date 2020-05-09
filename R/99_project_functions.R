# Define project functions
# ------------------------------------------------------------------------------


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