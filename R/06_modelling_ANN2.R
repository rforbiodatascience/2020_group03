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


# Hello from ANN2 and set folder
# ------------------------------------------------------------------------------
folder = "./results/06_ANN"
print('06_modelling_ANN2.R...')


# Data set 1 - ANN2 artificial neural network 
# ------------------------------------------------------------------------------
ANN2(
  name = "data_set_1",
  df = data_set_1,
  folder = folder,
  epochs = 10,  # number of epochs
  hidden_layers = c(100, 100),
  scale = "z_scales",
  train_size = 0.75,
  seed_value = 42)

ANN2(
  name = "data_set_1",
  df = data_set_1,
  folder = folder,
  epochs = 10,  # number of epochs
  hidden_layers = c(100, 100),
  scale = "blosum62",
  train_size = 0.75,
  seed_value = 42)


# Data set 2 - ANN2 artificial neural network 
# ------------------------------------------------------------------------------
ANN2(
  name = "data_set_2",
  df = data_set_2,
  folder = folder,
  epochs = 10,  # number of epochs
  hidden_layers = c(100, 100),
  scale = "z_scales",
  train_size = 0.75,
  seed_value = 42)

ANN2(
  name = "data_set_2",
  df = data_set_2,
  folder = folder,
  epochs = 10,  # number of epochs
  hidden_layers = c(100, 100),
  scale = "blosum62",
  train_size = 0.75,
  seed_value = 42)


# Data set 3 - ANN2 artificial neural network 
# ------------------------------------------------------------------------------
ANN2(
  name = "data_set_3",
  df = data_set_3,
  folder = folder,
  epochs = 10,  # number of epochs
  hidden_layers = c(100, 100),
  scale = "z_scales",
  train_size = 0.75,
  seed_value = 42)

ANN2(
  name = "data_set_3",
  df = data_set_3,
  folder = folder,
  epochs = 10,  # number of epochs
  hidden_layers = c(100, 100),
  scale = "blosum62",
  train_size = 0.75,
  seed_value = 42)


# Data set 4 - ANN2 artificial neural network 
# ------------------------------------------------------------------------------
ANN2(
  name = "data_set_4",
  df = data_set_4,
  folder = folder,
  epochs = 10,  # number of epochs
  hidden_layers = c(100, 100),
  scale = "z_scales",
  train_size = 0.75,
  seed_value = 42)

ANN2(
  name = "data_set_4",
  df = data_set_4,
  folder = folder,
  epochs = 10,  # number of epochs
  hidden_layers = c(100, 100),
  scale = "blosum62",
  train_size = 0.75,
  seed_value = 42)
