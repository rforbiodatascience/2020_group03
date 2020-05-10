# 04_heatmaps.R
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

# calculate effective amino acids

data_set_1_eff <- data_set_1 %>%
  mutate(effective = case_when(score<=0.3 ~ 0,
                               score>0.3 ~1)) %>%
  group_by(mutation_position) %>%
  summarise(N_eff = sum(effective))


neff_1 <- data_set_1_eff %>%
  ggplot(aes(x=mutation_position, y = N_eff)) +
  geom_line()
        

data_set_2_eff <- data_set_2 %>%
  mutate(effective = case_when(score<=4 ~ 0,
                               score>4 ~1)) %>%
  group_by(mutation_position) %>%
  summarise(N_eff = sum(effective))


neff_2 <- data_set_2_eff %>%
  ggplot(aes(x=mutation_position, y = N_eff)) +
  geom_line()
                        

