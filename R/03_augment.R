# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")

# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_functions.R")

# Load data
# ------------------------------------------------------------------------------
my_data_clean_aug <- read_table2(file = "data/clean_test.txt", col_names = FALSE)

# Wrangle data
# ------------------------------------------------------------------------------

# Remove non-important columns
my_data_clean_aug <- my_data_clean_aug %>% 
  rename(
    peptide = X1,
    activity = X2,
  ) %>%  select(c(1,2))

# Add new columns (encoding)

pepts <- my_data_clean_aug %>% 
  pull(peptide) %>% 
  encode_peptide("blosum62")

# Write data
# ------------------------------------------------------------------------------
#write_tsv(x = my_data_clean_aug,
#         path = "data/03_my_data_clean_aug.tsv")