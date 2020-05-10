# 04_density_plots.R
# ------------------------------------------------------------------------------

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
density_data_set_1 <- density_plot(data_set_1)
density_data_set_2 <- density_plot(data_set_2)
density_data_set_3 <- density_plot(data_set_3)
density_data_set_4 <- density_plot(data_set_4)

# Generate density plot per mutation
# ------------------------------------------------------------------------------
density_residue_data_set_1 <- density_plot_residue(data_set_1)
density_residue_data_set_2 <- density_plot_residue(data_set_2)
density_residue_data_set_3 <- density_plot_residue(data_set_3)
density_residue_data_set_4 <- density_plot_residue(data_set_4)