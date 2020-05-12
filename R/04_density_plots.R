# 04_density_plots.R
# ------------------------------------------------------------------------------
print('04_density_plots.R -> density plotting data')

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

# Generate density plots
# ------------------------------------------------------------------------------
density_data_set_1 <- density_plot(data_set_1, title="Score density plot for Data Set 1")
density_data_set_2 <- density_plot(data_set_2, title="Score density plot for Data Set 2")
density_data_set_3 <- density_plot(data_set_3, title="Score density plot for Data Set 3")
density_data_set_4 <- density_plot(data_set_4, title="Score density plot for Data Set 4")

# Generate density plot per mutation
# ------------------------------------------------------------------------------
density_residue_data_set_1 <- density_plot_residue(data_set_1)
density_residue_data_set_2 <- density_plot_residue(data_set_2)
density_residue_data_set_3 <- density_plot_residue(data_set_3)
density_residue_data_set_4 <- density_plot_residue(data_set_4)

# Save density plots
# ------------------------------------------------------------------------------
ggsave(plot = density_data_set_1, "./replace/density_plots/density_data_set_1.png", width = 7, height = 5, dpi=300)
ggsave(plot = density_data_set_2, "./replace/density_plots/density_data_set_2.png", width = 7, height = 5, dpi=300)
ggsave(plot = density_data_set_3, "./replace/density_plots/density_data_set_3.png", width = 7, height = 5, dpi=300)
ggsave(plot = density_data_set_4, "./replace/density_plots/density_data_set_4.png", width = 7, height = 5, dpi=300)

# Save density plots per mutation
# ------------------------------------------------------------------------------
ggsave(plot = density_residue_data_set_1, "./replace/density_plots/density_residue_data_set_1.png", width = 7, height = 6.25, dpi=300)
ggsave(plot = density_residue_data_set_2, "./replace/density_plots/density_residue_data_set_2.png", width = 7, height = 6.25, dpi=300)
ggsave(plot = density_residue_data_set_3, "./replace/density_plots/density_residue_data_set_3.png", width = 7, height = 6.25, dpi=300)
ggsave(plot = density_residue_data_set_4, "./replace/density_plots/density_residue_data_set_4.png", width = 7, height = 6.25, dpi=300)