# 04_seq_plots.R
# ------------------------------------------------------------------------------



# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library("ggplot2")
library("UniprotR")
library("ggseqlogo")

# Define functions
# ------------------------------------------------------------------------------
source(file = "./R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
data_set_1 <- read_tsv(file = "./data/03_aug_data_set_1.tsv")
data_set_2 <- read_tsv(file = "./data/03_aug_data_set_2.tsv")
data_set_3 <- read_tsv(file = "./data/03_aug_data_set_3.tsv")
data_set_4 <- read_tsv(file = "./data/03_aug_data_set_4.tsv")

data_set <- data_set_3 %>%
  drop_na() %>%
  mutate(score = (score - min(score)) / (max((score)-min(score))))

score_sum <- data_set %>%
  group_by(mutation_position) %>%
  summarise(score_sum=sum(score))

gg_seq_1 <- data_set %>%
  full_join(score_sum, by = "mutation_position") %>%
  filter(mutation_position >110 & mutation_position < 150) %>%
  mutate(score = score / score_sum) %>%
  select(mutation_position, mutation, score) %>%
  replace(is.na(.), 0) %>%
  pivot_wider(names_from = mutation_position, values_from = score)

custom_mat = as.matrix(gg_seq_1)
aa=custom_mat[,1]
colnames(custom_mat) <- NULL
row.names(custom_mat) <- aa
custom_mat<-custom_mat[,-1]

# Generate sequence logo
logo <- ggseqlogo(custom_mat, method='custom', seq_type='aa') + ylab('scaled scoring')
  
ggseq <- logo + theme(axis.text.x = element_blank())

ggsave(plot = ggseq, "./doc/sequence_logos/seq_logo_data_set_3.png")