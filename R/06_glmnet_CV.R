# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library("glmnet")
library("caret")
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



# Select dataset and wrangle
# ------------------------------------------------------------------------------
folder= "./results/06_glmnet"
print('06_modelling_CV_glmnet.R...')



# Get sequence window and truncate
# ------------------------------------------------------------------------------
glmnet_CV(
  df = data_set_1,
  name = "data_set_1",
  folder = folder,
  scale = "z_scales",
  train_size = 0.75,
  seed_value = 42,
  alpha_ = 0.2)

glmnet_CV(
  df = data_set_1,
  name = "data_set_1",
  folder = folder,
  scale = "blosum62",
  train_size = 0.75,
  seed_value = 42,
  alpha_ = 0.2)



glmnet_CV(
  df = data_set_2,
  name = "data_set_2",
  folder = folder,
  scale = "z_scales",
  train_size = 0.75,
  seed_value = 42,
  alpha_ = 0.2)

glmnet_CV(
  df = data_set_2,
  name = "data_set_2",
  folder = folder,
  scale = "blosum62",
  train_size = 0.75,
  seed_value = 42,
  alpha_ = 0.2)



glmnet_CV(
  df = data_set_3,
  name = "data_set_3",
  folder = folder,
  scale = "z_scales",
  train_size = 0.75,
  seed_value = 42,
  alpha_ = 0.2)

glmnet_CV(
  df = data_set_3,
  name = "data_set_3",
  folder = folder,
  scale = "blosum62",
  train_size = 0.75,
  seed_value = 42,
  alpha_ = 0.2)



glmnet_CV(
  df = data_set_4,
  name = "data_set_4",
  folder = folder,
  scale = "z_scales",
  train_size = 0.75,
  seed_value = 42,
  alpha_ = 0.2)

glmnet_CV(
  df = data_set_4,
  name = "data_set_4",
  folder = folder,
  scale = "z_scales",
  train_size = 0.75,
  seed_value = 42,
  alpha_ = 0.2)
