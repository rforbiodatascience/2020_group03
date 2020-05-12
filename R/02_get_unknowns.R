# 02_get_unknowns.R
# ------------------------------------------------------------------------------
print('02_clean.R -> cleaning data')

# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library("dplyr")

# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
data_set_1 <- read_tsv(file = "./data/01_load_data_set_1.tsv")
data_set_2 <- read_tsv(file = "./data/01_load_data_set_2.tsv")
data_set_3 <- read_tsv(file = "./data/01_load_data_set_3.tsv")
data_set_4 <- read_tsv(file = "./data/01_load_data_set_4.tsv")
amino_acids <- read_csv("data/_raw/amino_acids.csv", col_names = TRUE)

# Wrangle data
# ------------------------------------------------------------------------------

# Data_set 1
# select variant_ID and activity to predict, and rename columns
data_set_1_unknowns <- data_set_1 %>%
  select(Variant_ID, E3_score) %>%
  rename(
    variant = Variant_ID,
    score = E3_score
  ) %>%
  filter(!str_detect(variant, "\\*")) %>%
  filter(is.na(score))



# Data_set 2
# select ERK2_Mutant and activity to predict, and rename columns
data_set_2_unknowns <- data_set_2 %>%
  select(ERK2_Mutant, SCH_Average) %>%
  rename(
    variant = ERK2_Mutant,
    score = SCH_Average
  ) %>%
  filter(is.na(score))


# Data_set 3
# convert three letter amino acid representation to oneletter, select hgvs_pro and score to predict, and rename columns
data_set_3_unknowns <- data_set_3 %>%
  select(hgvs_pro, score) %>%
  rename(variant = hgvs_pro) %>%
  mutate(
    mutated_residue_3_letter = str_extract(variant, "[A-Z][a-z][a-z]"),
    mutation_position = str_extract(variant, "[0-9]+"),
    mutation_3_letter = str_extract(variant, "[A-Z,a-z,=]+$"),
    mutation_3_letter = str_replace(mutation_3_letter, "[=]", mutated_residue_3_letter)
  ) %>%
  full_join(
    amino_acids,
    by = c("mutation_3_letter" = "threeletter")
  ) %>%
  rename(mutation = oneletter) %>%
  full_join(
    amino_acids,
    by = c("mutated_residue_3_letter" = "threeletter")
  ) %>%
  rename(mutated_residue = oneletter) %>%
  unite(variant, mutated_residue, mutation_position, mutation, sep = "") %>%
  select(variant, score) %>%
  filter(!grepl('NA', variant)) %>%
  filter(is.na(score))

# Data_set 4
# pivot to longer table, rename columns
data_set_4_unknowns <- data_set_4 %>%
  pivot_longer(-X1, names_to = "substitution", values_to = "score") %>%
  unite("variant", X1:substitution, sep = "") %>%
  filter(!str_detect(variant, "\\*"))  %>%
  filter(is.na(score))


# Write data
# ------------------------------------------------------------------------------
write_tsv(x = data_set_1_unknowns, path = "./data/02_data_set_1_unknowns.tsv")
write_tsv(x = data_set_2_unknowns, path = "./data/02_data_set_2_unknowns.tsv")
write_tsv(x = data_set_3_unknowns, path = "./data/02_data_set_3_unknowns.tsv")
write_tsv(x = data_set_4_unknowns, path = "./data/02_data_set_4_unknowns.tsv")
