# 04_heatmaps.R
# ------------------------------------------------------------------------------
print('04_heatmap.R -> heatmapping data')


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

# Generate heatmaps
# ------------------------------------------------------------------------------
heatmap_data_set_1 <- heatmap_score(data_set_1, score_min = 0, score_max = 1.5) 
heatmap_data_set_2 <- heatmap_score(data_set_2, score_min = 1, score_max = 10)
heatmap_data_set_3 <- heatmap_score(data_set_3, score_min = 0, score_max = 1)
heatmap_data_set_4 <- heatmap_score(data_set_4, score_min = -2, score_max = 2)

# Save heatmaps
# ------------------------------------------------------------------------------
ggsave(plot = heatmap_data_set_1, "./doc/heatmaps/heatmap_data_set_1.png")
ggsave(plot = heatmap_data_set_2, "./doc/heatmaps/heatmap_data_set_2.png")
ggsave(plot = heatmap_data_set_3, "./doc/heatmaps/heatmap_data_set_3.png")
ggsave(plot = heatmap_data_set_4, "./doc/heatmaps/heatmap_data_set_4.png")