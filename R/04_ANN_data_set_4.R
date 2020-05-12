# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library("caret")
library("UniprotR")
library("ANN2")
library("yardstick")


# Define functions
# ------------------------------------------------------------------------------
source(file = "./R/99_project_functions.R")


# Load data
# ------------------------------------------------------------------------------
data_set_1 <- read_tsv(file = "./data/03_aug_data_set_1.tsv")
data_set_2 <- read_tsv(file = "./data/03_aug_data_set_2.tsv")
data_set_3 <- read_tsv(file = "./data/03_aug_data_set_3.tsv")
data_set_4 <- read_tsv(file = "./data/03_aug_data_set_4.tsv")




# Select dataset and wrangle
# ------------------------------------------------------------------------------
test_set = data_set_4
filename = "./results/06_ann/ann_corr_data_set_4.png"
RMSE_path = "./data/ANN_zs_RMSE_data_set_4.tsv"


data <- test_set %>%
  drop_na(score)


# ANN2 hyperparameters
# ------------------------------------------------------------------------------
epochs = 10 # number of epochs
hidden_layers = c(100, 100) # number and size of hidden layers
scale = "z_scales" # amino acid descriptor scale
train_size = 0.75 # train sample size
seed_value = 42 # random seed

rm <- ANN2(data, filename, epochs, hidden_layers, scale = "z_scales", train_size, seed_value)