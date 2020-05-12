# 04_density_plots.R
# ------------------------------------------------------------------------------
print('04_geatmaps.R -> plotting heatmaps')

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

# Score heatmaps
# ------------------------------------------------------------------------------
heatmap_data_set_score_1 <- heatmap_score(data_set_1, score_min = 0, score_max = 1.5)
heatmap_data_set_score_2 <- heatmap_score(data_set_2, score_min = 1, score_max = 10)
heatmap_data_set_score_3 <- heatmap_score(data_set_3, score_min = 0, score_max = 1)
heatmap_data_set_score_4 <- heatmap_score(data_set_4, score_min = -2, score_max = 2)

# Save heatmaps
# ------------------------------------------------------------------------------
ggsave(plot = heatmap_data_set_score_1, "./results/heatmaps/heatmap_data_set_score_1.png",width = 10, height = 3, dpi=300)
ggsave(plot = heatmap_data_set_score_2, "./results/heatmaps/heatmap_data_set_score_2.png",width = 10, height = 3, dpi=300)
ggsave(plot = heatmap_data_set_score_3, "./results/heatmaps/heatmap_data_set_score_3.png",width = 10, height = 3, dpi=300)
ggsave(plot = heatmap_data_set_score_4, "./results/heatmaps/heatmap_data_set_score_4.png",width = 10, height = 3, dpi=300)


# Add aminoacid type / call function aminoacid_type
data_set_1 <- data_set_1 %>%
  mutate(aminoacid_class = aminoacid_type(mutation))
data_set_2 <- data_set_2 %>%
  mutate(aminoacid_class = aminoacid_type(mutation))
data_set_3 <- data_set_3 %>%
  mutate(aminoacid_class = aminoacid_type(mutation))
data_set_4 <- data_set_4 %>%
  mutate(aminoacid_class = aminoacid_type(mutation))


# Score grouped heatmaps
# ------------------------------------------------------------------------------
heatmap_data_set_score_g_1 <- data_set_1 %>% 
  group_by(aminoacid_class) %>%
  heatmap_g_score(score_min = 1, score_max = 10 )
heatmap_data_set_score_g_2 <- data_set_2 %>% 
  group_by(aminoacid_class) %>% 
  heatmap_g_score(score_min = 1, score_max = 10)
heatmap_data_set_score_g_3 <- data_set_3 %>% 
  group_by(aminoacid_class) %>% 
  heatmap_g_score(score_min = 0, score_max = 1)
heatmap_data_set_score_g_4 <- data_set_4 %>% 
  group_by(aminoacid_class) %>% 
  heatmap_g_score(score_min = -2, score_max = 2)

# Save heatmaps
# ----------------------------------------------------------------------------
ggsave(plot = heatmap_data_set_score_g_1, "./results/heatmaps/heatmap_data_set_score_g_1.png",width = 10, height = 3, dpi=300)
ggsave(plot = heatmap_data_set_score_g_2, "./results/heatmaps/heatmap_data_set_score_g_2.png",width = 10, height = 3, dpi=300)
ggsave(plot = heatmap_data_set_score_g_3, "./results/heatmaps/heatmap_data_set_score_g_3.png",width = 10, height = 3, dpi=300)
ggsave(plot = heatmap_data_set_score_g_4, "./results/heatmaps/heatmap_data_set_score_g_4.png",width = 10, height = 3, dpi=300)
