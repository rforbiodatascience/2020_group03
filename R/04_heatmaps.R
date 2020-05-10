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

# Heatmap 
heatmap_data_set_1 <- ggplot(data_set_1, aes(mutation_position, mutation, fill = score)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_fill_continuous(limits=c(0, 1.5),type ="viridis", oob= scales::squish) +
  theme(panel.background = element_blank()) +
  geom_tile()

heatmap_data_set_2 <- ggplot(data_set_2, aes(mutation_position, mutation, fill = score)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_fill_continuous(type = "viridis") +
  scale_fill_continuous(limits=c(1, 10),type ="viridis", oob= scales::squish) +
  geom_tile()

heatmap_data_set_3 <- ggplot(data_set_3, aes(mutation_position, mutation, fill = score)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_fill_continuous(limits=c(0, 1),type ="viridis", oob= scales::squish) +
  theme(panel.background = element_blank()) +
  geom_tile()

heatmap_data_set_4 <- ggplot(data_set_4, aes(mutation_position, mutation, fill = score)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_fill_continuous(limits=c(-2, 2),type ="viridis", oob= scales::squish) +
  theme(panel.background = element_blank()) +
  geom_tile()
  
distribution_data_set_1 <- ggplot(data_set_1, aes(score)) +
  geom_density(alpha=0.4, fill="lightblue")

distribution_data_set_2 <- ggplot(data_set_2, aes(score)) +
  geom_density(alpha=0.4, fill="lightblue")

distribution_data_set_3 <- ggplot(data_set_3, aes(score)) +
  geom_density(alpha=0.4, fill="lightblue")

distribution_data_set_4 <- ggplot(data_set_4, aes(score)) +
  geom_density(alpha=0.4, fill="lightblue")
  
distribution_data_set_1_c <- ggplot(data_set_1, aes(score, color=mutated_residue)) +
  geom_density(alpha=0.4)

distribution_data_set_2_c <- ggplot(data_set_2, aes(score, color=mutated_residue)) +
  geom_density(alpha=0.4)

distribution_data_set_3_c <- ggplot(data_set_3, aes(score, color=mutated_residue)) +
  geom_density(alpha=0.4)

distribution_data_set_4_c <- ggplot(data_set_4, aes(score, color=mutated_residue)) +
  geom_density(alpha=0.4)

amino_acids <- read_csv("data/_raw/amino_acid_types.csv", col_names = TRUE)

heatmap_type_1 <- data_set_1 %>%
  full_join(amino_acids, by =c("mutation" = "Amino_acid"))

heatmap_type_data_set_1 <- ggplot(heatmap_type_1, aes(mutation_position, Type, fill = score)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_fill_continuous(type ="viridis", oob= scales::squish) +
  theme(panel.background = element_blank()) +
  geom_tile()


