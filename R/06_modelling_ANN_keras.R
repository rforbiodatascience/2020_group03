# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library("devtools")
library("keras")
library("ggplot2")

# Define functions
# ------------------------------------------------------------------------------
source(file = "./R/99_project_functions.R")

# Data set 1 
# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------
data_set_1 <- read_tsv(file = "./data/03_aug_data_set_1.tsv")

ANN_keras_func(
  name = "data_set_1",
  df = data_set_1,
  folderplot = "./results/06_ANN_plots/",
  folderdata = "./results/06_ANN_data/",
  epochs = 10,  # number of epochs
  scale = "blosum62",
  partitions = 0.30,
  batch_size = 200,
  learning_rate = 0.1
)

# Data set 2
# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------
data_set_2 <- read_tsv(file = "./data/03_aug_data_set_2.tsv")

ANN_keras_func(
  name = "data_set_2",
  df = data_set_2,
  folderplot = "./results/06_ANN_plots/",
  folderdata = "./results/06_ANN_data/",
  epochs = 10,  # number of epochs
  scale = "blosum62",
  partitions = 0.30,
  batch_size = 200,
  learning_rate = 0.1
)

# Data set 3 
# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------
data_set_3 <- read_tsv(file = "./data/03_aug_data_set_3.tsv")

ANN_keras_func(
  name = "data_set_3",
  df = data_set_3,
  folderplot = "./results/06_ANN_plots/",
  folderdata = "./results/06_ANN_data/",
  epochs = 10,  # number of epochs
  scale = "blosum62",
  partitions = 0.30,
  batch_size = 200,
  learning_rate = 0.1
)

# Data set 4 
# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------
data_set_4 <- read_tsv(file = "./data/03_aug_data_set_4.tsv")

ANN_keras_func(
  name = "data_set_4",
  df = data_set_4,
  folderplot = "./results/06_ANN_plots/",
  folderdata = "./results/06_ANN_data/",
  epochs = 10,  # number of epochs
  scale = "blosum62",
  partitions = 0.30,
  batch_size = 200,
  learning_rate = 0.1
)
