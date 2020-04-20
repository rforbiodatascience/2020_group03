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
# load amino acid scales
aa_scale_z5 <- load_scale('z5')
aa_scale_bl62 <- load_scale('bl62')

# load peptide sequence and results
peptides <- read_csv('./data/peptides.txt')


# Wrangle data
# ------------------------------------------------------------------------------
# peptide sequences
X <- peptides %>%
  pull(Sequence) %>%
  encode_peptide(m = aa_scale_bl62)

y <- peptides %>%
  pull(Result)

Sequences <- peptides %>%
  pull(Sequence)

p1 = ggseqlogo(Sequences)
print(p1)

# Write data
# ------------------------------------------------------------------------------
#write_tsv(x = my_data,
  #        path = "data/01_my_data.tsv")
