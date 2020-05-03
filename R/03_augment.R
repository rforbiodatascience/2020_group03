# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library("UniprotR")

# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_functions.R")

# Load data
# ------------------------------------------------------------------------------
my_data_clean_aug <- read_table2(file = "data/clean_test.txt", col_names = FALSE)
amino_acids <- read_csv("data/_raw/amino_acids.csv", col_names = TRUE)

# Wrangle data
# ------------------------------------------------------------------------------

# Remove non-important columns
# my_data_clean_aug <- my_data_clean_aug %>% 
#   rename(
#     peptide = X1,
#     activity = X2,
#   ) %>%  select(c(1,2))

# Add new columns (encoding)

#pepts <- my_data_clean_aug %>% 
#  pull(peptide) %>% 
#  encode_peptide("blosum62")

# Write data
# ------------------------------------------------------------------------------
#write_tsv(x = my_data_clean_aug,
#         path = "data/03_my_data_clean_aug.tsv")


# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data
# ------------------------------------------------------------------------------
clean_data_set_1 <- read_tsv(file = "./data/02_clean_data_set_1.tsv")
clean_data_set_2 <- read_tsv(file = "./data/02_clean_data_set_2.tsv")
clean_data_set_3 <- read_tsv(file = "./data/02_clean_data_set_3.tsv")
clean_data_set_4 <- read_tsv(file = "./data/02_clean_data_set_4.tsv")


# Wrangle data
# ------------------------------------------------------------------------------

# need a function

aug_data_set_1 <- clean_data_set_1 %>%
  mutate(
    mutated_residue = get_mutated_residue(variant), 
    mutation_position = get_mutation_position(variant),
    mutation = get_mutation(variant),
    sequence = make_sequence(mutated_residue, mutation_position, mutation, "P38398"))

aug_data_set_2 <- clean_data_set_2 %>%
  mutate(
    mutated_residue = get_mutated_residue(variant), 
    mutation_position = get_mutation_position(variant),
    mutation = get_mutation(variant),
    sequence = make_sequence(mutated_residue, mutation_position, mutation, "P28482"))

aug_data_set_3 <- clean_data_set_3 %>%
  mutate(
    mutated_residue_3_letter = str_extract(variant, "[A-Z][a-z][a-z]"), 
    mutation_position = str_extract(variant, "[0-9]+"),
    mutation_3_letter = str_extract(variant, "[A-Z,a-z,*]+$"),
    mutation_3_letter = str_replace(mutation_3_letter, "[*]", mutated_residue_3_letter)) %>%
  inner_join(
    amino_acids, by = c("mutation_3_letter" = "threeletter")) %>%
  rename(mutation = oneletter) %>%
  inner_join(
    amino_acids, by = c("mutated_residue_3_letter" = "threeletter")) %>%
  rename(mutated_residue = oneletter) %>%
  mutate(
    Sequence = make_sequence(mutated_residue, mutation_position, mutation, "Q5SW96"))

aug_data_set_4 <- clean_data_set_4 %>%
  mutate(
    mutated_residue = get_mutated_residue(variant), 
    mutation_position = get_mutation_position(variant),
    mutation = get_mutation(variant),
    sequence = make_sequence(mutated_residue, mutation_position, mutation, "P04147"))


# Write data
# ------------------------------------------------------------------------------
write_tsv(x = aug_data_set_1,
          path = "./data/03_aug_data_set_1.tsv")

write_tsv(x = aug_data_set_2,
          path = "./data/03_aug_data_set_2.tsv")

write_tsv(x = aug_data_set_3,
          path = "./data/03_aug_data_set_3.tsv")

write_tsv(x = aug_data_set_4,
          path = "./data/03_aug_data_set_4.tsv")