# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library("glmnet")
library("caret")
library("UniprotR")
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


# Hello from glmnet  and set folder
# ------------------------------------------------------------------------------
folder = "./results/06_glmnet"
print('06_modelling_glmnet.R...')


# Data set 1 - glmnet 
# ------------------------------------------------------------------------------
glmnet_reg(
  name = "data_set_1",
  df = data_set_1,
  folder = folder,
  alpha_ = 0.2, # mix between ridge and lasso
  s_ = 0.01, # regularization strenght
  scale = "z_scales", # amino acid descriptor scale
  train_size = 0.75, # train sample size
  seed_value = 42) # random seed

glmnet_reg(
  name = "data_set_1",
  df = data_set_1,
  folder = folder,
  alpha_ = 0.2, # mix between ridge and lasso
  s_ = 0.01, # regularization strenght
  scale = "blosum62", # amino acid descriptor scale
  train_size = 0.75, # train sample size
  seed_value = 42) # random seed


# Data set 2 - glmnet 
# ------------------------------------------------------------------------------
glmnet_reg(
  name = "data_set_2",
  df = data_set_2,
  folder = folder,
  alpha_ = 0.2, # mix between ridge and lasso
  s_ = 0.01, # regularization strenght
  scale = "z_scales", # amino acid descriptor scale
  train_size = 0.75, # train sample size
  seed_value = 42) # random seed

glmnet_reg(
  name = "data_set_2",
  df = data_set_2,
  folder = folder,
  alpha_ = 0.2, # mix between ridge and lasso
  s_ = 0.01, # regularization strenght
  scale = "blosum62", # amino acid descriptor scale
  train_size = 0.75, # train sample size
  seed_value = 42) # random seed


# Data set 3 - glmnet 
# ------------------------------------------------------------------------------
glmnet_reg(
  name = "data_set_3",
  df = data_set_3,
  folder = folder,
  alpha_ = 0.2, # mix between ridge and lasso
  s_ = 0.01, # regularization strenght
  scale = "z_scales", # amino acid descriptor scale
  train_size = 0.75, # train sample size
  seed_value = 42) # random seed

glmnet_reg(
  name = "data_set_3",
  df = data_set_3,
  folder = folder,
  alpha_ = 0.2, # mix between ridge and lasso
  s_ = 0.01, # regularization strenght
  scale = "blosum62", # amino acid descriptor scale
  train_size = 0.75, # train sample size
  seed_value = 42) # random seed


# Data set 4 - glmnet 
# ------------------------------------------------------------------------------
glmnet_reg(
  name = "data_set_4",
  df = data_set_4,
  folder = folder,
  alpha_ = 0.2, # mix between ridge and lasso
  s_ = 0.01, # regularization strenght
  scale = "z_scales", # amino acid descriptor scale
  train_size = 0.75, # train sample size
  seed_value = 42) # random seed

glmnet_reg(
  name = "data_set_4",
  df = data_set_4,
  folder = folder,
  alpha_ = 0.2, # mix between ridge and lasso
  s_ = 0.01, # regularization strenght
  scale = "blosum62", # amino acid descriptor scale
  train_size = 0.75, # train sample size
  seed_value = 42) # random seed
