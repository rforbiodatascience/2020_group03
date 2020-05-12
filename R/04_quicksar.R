# 04_quicksar.R
# ------------------------------------------------------------------------------
print('04_quicksar.R -> quickSARing data')

# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library("ggplot2")
library("UniprotR")

# Define functions
# ------------------------------------------------------------------------------
source(file = "./R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
data_set_1 <- read_tsv(file = "./data/03_aug_data_set_1.tsv")
data_set_2 <- read_tsv(file = "./data/03_aug_data_set_2.tsv")
data_set_3 <- read_tsv(file = "./data/03_aug_data_set_3.tsv")
data_set_4 <- read_tsv(file = "./data/03_aug_data_set_4.tsv")

data_set_4 <- data_set_4 %>%
  drop_na(score)

data_set_3 <- data_set_3 %>%
  drop_na(score)

# Plot effective amino acids at each position with a cutoff
# ------------------------------------------------------------------------------
quick_sar_data_set_1 <- quick_sar(data_set_1, cutoff = 0.5)
quick_sar_data_set_2 <- quick_sar(data_set_2, cutoff = 3)
quick_sar_data_set_3 <- quick_sar(data_set_3, cutoff = 0.6)
quick_sar_data_set_4 <- quick_sar(data_set_4, cutoff = 0)

# Save quickSAR plots
# ------------------------------------------------------------------------------
ggsave(plot = quick_sar_data_set_1, "./results/04_quick_SAR/quick_sar_data_set_1.png",width = 10, height = 3, dpi=300)
ggsave(plot = quick_sar_data_set_2, "./results/04_quick_SAR/quick_sar_data_set_2.png",width = 10, height = 3, dpi=300)
ggsave(plot = quick_sar_data_set_3, "./results/04_quick_SAR/quick_sar_data_set_3.png",width = 10, height = 3, dpi=300)
ggsave(plot = quick_sar_data_set_4, "./results/04_quick_SAR/quick_sar_data_set_4.png",width = 10, height = 3, dpi=300)