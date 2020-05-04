# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library("ggseqlogo")
library("protr")
library("tensorflow")
library("keras")

# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
aug_clean_data_set_2 <- read_tsv(file = "./data/03_aug_data_set_2.tsv")

# Wrangle data
# ------------------------------------------------------------------------------
df <- aug_clean_data_set_2 %>%
  filter(score>10)


bl62 = read.table(file = 'data/_raw/BLOSUM62_.txt')

X_train <- df %>%
  pull(sequence) %>%
  encode_peptide(m = bl62)

y_train = df %>%
  pull(score)



# Visualise data
# ---------------------------
# ht <- df %>%
#     ggplot(aes(x = score)) +
#     geom_histogram()

# Write data
# ------------------------------------------------------------------------------
#write_tsv(...)
#ggsave(...)
