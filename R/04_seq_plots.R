# 04_seq_plots.R
# ------------------------------------------------------------------------------

# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library("ggplot2")
library("UniprotR")
library("ggseqlogo")

# Define functions
# ------------------------------------------------------------------------------
source(file = "./R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
data_set_1 <- read_tsv(file = "./data/03_aug_data_set_1.tsv")
data_set_2 <- read_tsv(file = "./data/03_aug_data_set_2.tsv")
data_set_3 <- read_tsv(file = "./data/03_aug_data_set_3.tsv")
data_set_4 <- read_tsv(file = "./data/03_aug_data_set_4.tsv")

# Make sequence logo plots
# ------------------------------------------------------------------------------
seq_logo_data_set_1 <- sequence_logo(data_set_1, start_position = 40, end_position = 80)
seq_logo_data_set_2 <- sequence_logo(data_set_1, start_position = 40, end_position = 80)
seq_logo_data_set_3 <- sequence_logo(data_set_1, start_position = 110, end_position = 150)
seq_logo_data_set_4 <- sequence_logo(data_set_1, start_position = 135, end_position = 185)

# Save sequence logos
# ------------------------------------------------------------------------------
ggsave(plot = seq_logo_data_set_1, "./results/04_sequence_logos/seq_logo_data_set_1.png",width = 10, height = 3, dpi=300)
ggsave(plot = seq_logo_data_set_2, "./results/04_sequence_logos/seq_logo_data_set_2.png",width = 10, height = 3, dpi=300)
ggsave(plot = seq_logo_data_set_3, "./results/04_sequence_logos/seq_logo_data_set_3.png",width = 10, height = 3, dpi=300)
ggsave(plot = seq_logo_data_set_4, "./results/04_sequence_logos/seq_logo_data_set_4.png",width = 10, height = 3, dpi=300)