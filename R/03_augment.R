# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library("UniprotR")
library("protr")

# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_functions.R")

# Load data
# ------------------------------------------------------------------------------
my_data_clean_aug <- read_table2(file = "data/clean_test.txt", col_names = FALSE)

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
aug_data_set_1 <- clean_data_set_1 %>%
  mutate(variant = case_when(variant=="NA-NA" ~ "D2D",
                             !variant=="NA-NA" ~ variant))
# augmenting with sequence
aug_data_set_1 <- convert_variant_to_sequence(aug_data_set_1, "P38398")
aug_data_set_2 <- convert_variant_to_sequence(clean_data_set_2, "P28482")
aug_data_set_3 <- convert_variant_to_sequence(clean_data_set_3, "Q5SW96")
aug_data_set_4 <- convert_variant_to_sequence(clean_data_set_3, "P04147")

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