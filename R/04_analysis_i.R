# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library('Peptides')

# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_functions.R")

# Load data
# ------------------------------------------------------------------------------
aug_clean_data_set_1 <- read_tsv(file = "./data/03_aug_data_set_1.tsv")

# Wrangle data
# ------------------------------------------------------------------------------
# Create new dataframe with ID, Seq, E3 score only
seq_properties <- aug_clean_data_set_1[,c(1,2,6)]

# Calculate molecular weight
seq_properties$MW <- aug_clean_data_set_1 %>%
  pull(Sequence) %>% 
  mw(monoisotopic = FALSE)

# Calculate pI of the sequences
seq_properties$pI <- aug_clean_data_set_1 %>%
  pull(Sequence) %>% 
  pI(pKscale = "EMBOSS")

# Calculate z-scale of the sequences
seq_properties[c("z1", "z2", "z3", "z4", "z5")] <- aug_clean_data_set_1 %>%
  pull(Sequence) %>% 
  zScales()


# Calculate proportion of peptides






# Visualise data
# ------------------------------------------------------------------------------