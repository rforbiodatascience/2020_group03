# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library("ggseqlogo")

# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
aug_clean_data_set_1 <- read_tsv(file = "./data/03_aug_data_set_1.tsv")

# Wrangle data
# ------------------------------------------------------------------------------
df <- aug_clean_data_set_1 %>%
  filter(score > 1)

# Model data
# ------------------------------------------------------------------------------
#my_data_clean_aug %>% ...

# Visualise data
# ---------------------------
ht <- df %>%
    ggplot(aes(x = score)) +
    geom_histogram()

seq <- df %>%
  select(sequence)

# Write data
# ------------------------------------------------------------------------------
#write_tsv(...)
#ggsave(...)
