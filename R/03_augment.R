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
data <- read_table2(file = "data/clean_test.txt", col_names = FALSE)

# Wrangle data
# ------------------------------------------------------------------------------

# Remove non-important columns
data <- data %>% 
  rename(
    peptide = X1,
    activity = X2,
  ) %>%  select(c(1,2))

data$len = str_length(data$peptide)
data <- data %>% group_by(len) %>% filter(len == 15) 

# Add new columns (encoding)
encoded_seq = data %>%
  pull(peptide) %>% 
  encode_peptide(m = "blosum62")

# Write data
# ------------------------------------------------------------------------------
#write_tsv(x = my_data_clean_aug,
#         path = "data/03_my_data_clean_aug.tsv")