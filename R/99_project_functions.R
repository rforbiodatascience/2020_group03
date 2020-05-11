# Define project functions
# ------------------------------------------------------------------------------

#' @description
#' `encode_peptide` returns a dataframe encoding the sequence in the method `m' specified.
#' All the peptides have to be the same length.
#'
#' @param x peptide, string
#' @param m method, string, options: blosum45, blosum63, blosum80, blosum90, pam30, pam70 and pam250
#' 
#' @return if all the parameters are logical will return a dataframe with the row.names = peptide and each column the encoding. The number of columns is len(peptide) x 20 amino acids. 
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
        t() %>%
        matrix(nrow = 1, byrow = TRUE) %>%
        return
    })
  X_enc = do.call(rbind, X_enc)
  rownames(X_enc) = x
  return(X_enc)
}


#' make_sequence
#' @description
#' `make_sequence` convert variant information and WT sequence to variant sequence
#' Supports only single mutations
#'
#' @param Mutated_residue, string residue to mutate
#' @param Mutation_position, numeric, residue number to mutate
#' @param Mutation, string mutation to which amino acid
#' @param uniprot_ID, string uniprot_ID
#' 
#' @return sequence, string protein sequence
#' 
#' @examples
#' make_sequence("D", 2, "K", "37636")
make_sequence = function(Mutated_residue, Mutation_position, Mutation, uniprot_ID)
{
  WT_sequence <- as_tibble(GetSequences(uniprot_ID)) %>%
    pull(Sequence)
  sequence <- WT_sequence
  str_sub(sequence, Mutation_position, Mutation_position) <- Mutation
  return(sequence)
}


#' @description
#' `get_mutated_residue` extracts inforationm about WT residue
#' Supports only single mutations
#'
#' @param variant, string format example "D23E"
#'
#' @return mutated_residue, string
#'
#' @examples
#' get_mutated_residue("D2K")
get_mutated_residue <- function(variant) {
  mutated_residue = str_extract(variant, "[A-Z]")
  return(mutated_residue)
}


#' @description
#' `get_mutation_position` extracts information about mutation position
#'
#' @param variant, string
#'
#' @return mutation_position, numeric
#'
#' @examples
#' get_mutation_position("D2K")
get_mutation_position <- function(variant) {
  mutation_position = str_extract(variant, "[0-9]+")
  return(mutation_position)
}


#' @description
#' `get_mutation` extracts information about mutation to which amino_acid
#'
#' @param variant, string
#'
#' @return mutation, string
#' @export
#'
#' @examples
#' get_mutation("D2K")
get_mutation <- function(variant) {
  mutation = str_extract(variant, "[A-Z,*]$")
  return(mutation)
}


#' @description
#' `convert_variant_to_sequence` converts a tibble of variant_IDs to a tibble with corrensponding
#' sequences based on a uniprot ID
#' 
#' @param df, tibble of variant_IDs column name = variant
#' @param protein, string uniprot ID
#'
#' @return tibble with sequences
#'
#' @examples convert_variant_to_sequence(df, "74742")
convert_variant_to_sequence <- function(df, protein)  {
  clean_df <- df %>%
  mutate(
    mutated_residue = get_mutated_residue(variant), 
    mutation_position = get_mutation_position(variant),
    mutation = get_mutation(variant),
    sequence = make_sequence(mutated_residue, mutation_position, mutation, protein))
  return(clean_df)
}

#' @description 
#' `heatmap_score` makes a ggplot2 tile plot based on a tibble
#' 
#'
#' @param df, tibble with columns score, mutation_position, mutation info
#' @param score_min, numeric lower cutoff
#' @param score_max, numeric higher cutoff
#'
#' @return heatmap object
#'
#' @examples heatmap_score(df, 0.2, 4.3)
heatmap_score <- function(df, score_min, score_max) {
  heatmap <- ggplot(df, aes(mutation_position, mutation, fill = score)) +
    scale_x_continuous(expand=c(0,0)) +
    scale_fill_continuous(limits=c(score_min, score_max),type ="viridis", oob= scales::squish) +
    theme(panel.background = element_blank()) +
    geom_tile() +
    theme_classic() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      axis.line = element_line(colour = "black", size = 0),
      axis.text.x = element_text(face = "bold", color = "#000000"),
      axis.text.y = element_text(face = "bold", color = "#000000")
    )
  return(heatmap)
}


#' @description 
#' `density_plot` makes a ggplot2 density plot based on a tibble
#' 
#' 
#' @param df, tibble with score column 
#'
#' @return density plot object
#'
#' @examples density_plot(df)
density_plot <- function(df, title) {
  density <- ggplot(df, aes(score)) +
    geom_density(alpha=0.4, fill="lightblue") +
    ggtitle(title) +
    theme_classic() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      axis.line = element_line(colour = "black", size = 0),
      axis.text.x = element_text(face = "bold", color = "#000000"),
      axis.text.y = element_text(face = "bold", color = "#000000")
    )
  return(density)
}

#' @description 
#' `density_plot_residue` makes a ggplot2 density plot based on a tibble
#' 
#' 
#' @param df, tibble with score column and mutated_residue columns
#'
#' @return density plot object
#'
#' @examples density_plot_residue(df)
density_plot_residue <- function(df){
  density_residue <- ggplot(df, aes(score, color=mutated_residue, fill=mutated_residue)) +
    geom_density(alpha=0.1) +
    theme_classic() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      axis.line = element_line(colour = "black", size = 0),
      axis.text.x = element_text(face = "bold", color = "#000000"),
      axis.text.y = element_text(face = "bold", color = "#000000")
    )
  return(density_residue)
}

#' @description 
#' `quick_sar` makes a ggplot2 line plot based on a tibble
#' 
#' @param df, tibble with score column and mutation_position columns
#'
#' @return density plot object
#'
#' @examples quick_sar(df)
quick_sar <- function(df, cutoff) {
  df_eff <- df %>%
    mutate(effective = case_when(score <= cutoff ~ 0,
                                 score > cutoff ~ 1)) %>%
    group_by(mutation_position) %>%
    summarise(N_eff = sum(effective))
  
  neff_plot <- df_eff %>%
    ggplot(aes(x=mutation_position, y = N_eff)) +
    geom_line() +
    theme_classic() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      axis.line = element_line(colour = "black", size = 0),
      axis.text.x = element_text(face = "bold", color = "#000000"),
      axis.text.y = element_text(face = "bold", color = "#000000")
    )
  return(neff_plot)
}