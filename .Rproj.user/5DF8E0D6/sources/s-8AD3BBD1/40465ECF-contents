# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library("UniprotR")

# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data
# ------------------------------------------------------------------------------
clean_data_set_1 <- read_tsv(file = "./data/02_clean_data_set_1.tsv")
clean_data_set_2 <- read_tsv(file = "./data/02_clean_data_set_2.tsv")


# Wrangle data
# ------------------------------------------------------------------------------
aug_data_set_1 <- clean_data_set_1 %>% 
  mutate(
    Mutated_residue = str_extract(Variant_ID, "[A-Z]"),
    Mutation_position = str_extract(Variant_ID, "[0-9]+"),
    Mutation = str_extract(Variant_ID, "[A-Z,*]$"),
    Mutation = str_replace(Mutation, "[*]", Mutated_residue),
    Sequence = make_sequence(Mutated_residue, Mutation_position, Mutation, "P38398"))

# aug_data_set_2 <- clean_data_set_2 %>% 
#   mutate(
#     Mutated_residue = str_extract(ERK2_Mutant, "[A-Z]"),
#     Mutation_position = str_extract(ERK2_Mutant, "[0-9]+"),
#     Mutation = str_extract(ERK2_Mutant, "[A-Z,*]$"),
#     Mutation = str_replace(Mutation, "[*]", Mutated_residue),
#     Sequence = make_sequence(Mutated_residue, Mutation_position, Mutation, "P38398"))


# Write data
# ------------------------------------------------------------------------------
write_tsv(x = aug_data_set_1,
          path = "./data/03_aug_data_set_1.tsv")