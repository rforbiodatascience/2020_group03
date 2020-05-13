# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library("stats")
library("broom")
library("dplyr")
library("ggpubr")

# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
# data_set_1 <- read_tsv(file = "./data/03_aug_data_set_1.tsv")
# data_set_2 <- read_tsv(file = "./data/03_aug_data_set_2.tsv")
# data_set_3 <- read_tsv(file = "./data/03_aug_data_set_3.tsv")
data_set_4 <- read_tsv(file = "./data/03_aug_data_set_4.tsv")

# Wrangle data
# ------------------------------------------------------------------------------
# Encode the sequences
# encoded_data_set_1 <- data_set_1 %>%
#  pull(sequence) %>%
#  encode_peptide(m = "blosum62")
# encoded_data_set_2 <- data_set_2 %>%
#  pull(sequence) %>%
#  encode_peptide(m = "blosum62")
# encoded_data_set_3 <- data_set_3 %>%
#  pull(sequence) %>%
#  encode_peptide(m = "blosum62")
encoded_data_set_4 <- data_set_4 %>%
  pull(sequence) %>%
  encode_peptide(m = "blosum62")

# Creation of two group labels based on activity measurements to use as a colour option later in plots
# data_set_1 <- data_set_1 %>%
#  mutate(activity_group = case_when(score >= 1.5 ~ 'High activity',
#                                    score < 1.5 ~ 'Low activity'))
# data_set_2 <- data_set_2 %>%
#  mutate(activity_group = case_when(score >= 11 ~ 'High activity',
#                                    score < 11 ~ 'Low activity'))
# data_set_3 <- data_set_3 %>%
#  mutate(activity_group = case_when(score >= 0.975 ~ 'High activity',
#                                    score < 0.975 ~ 'Low activity'))
data_set_4 <- data_set_4 %>%
  mutate(activity_group = case_when(
    score >= 1.5 ~ "High activity",
    score < 1.5 ~ "Low activity"
  ))

# Visualise data before PCA
# ------------------------------------------------------------------------------
# Activity distribution before and after the mean computation for common strings
distribution_4 <- data_set_4 %>%
  ggplot(aes(x = score)) +
  geom_density() +
  geom_vline(aes(xintercept = mean(score)),
    linetype = "dashed", size = 1
  ) +
  theme_classic() +
  labs(x = "Activity", y = "Density") +
  ggtitle("Activity distribution") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black", size = 0),
    axis.text.x = element_text(face = "bold", color = "#000000"),
    axis.text.y = element_text(face = "bold", color = "#000000")
  )

# Save distribution plots
# ------------------------------------------------------------------------------
# ggsave(plot = distribution_1, "./results/05_PCA_plots/PCA_distribution_1.png", width = 5.5, height = 5, dpi = 300)
# ggsave(plot = distribution_2, "./results/05_PCA_plots/PCA_distribution_2.png", width = 5.5, height = 5, dpi = 300)
# ggsave(plot = distribution_3, "./results/05_PCA_plots/PCA_distribution_3.png", width = 5.5, height = 5, dpi = 300)
ggsave(plot = distribution_4, "./results/05_PCA_plots/PCA_distribution_4.png", width = 5.5, height = 5, dpi = 300)

# PCA
# ------------------------------------------------------------------------------
# Compute PCA on encoded sequences -- zero centered values and unit variance
# pca_1 <- encoded_data_set_1 %>%
#  prcomp(center = TRUE)
# pca_2 <- encoded_data_set_2 %>%
#  prcomp(center = TRUE)
# pca_3 <- encoded_data_set_3 %>%
#  prcomp(center = TRUE)
pca_4 <- encoded_data_set_4 %>%
  prcomp(center = TRUE)

# Check results
# pca_4 %>% tidy("pcs")

# PCA -Scree plot
scree_plot_pca_4 <- pca_4 %>%
  tidy("pcs") %>%
  ggplot() +
  geom_line(aes(x = PC, y = percent)) +
  ggtitle("Variance explained for the principal components") +
  theme_classic() +
  labs(x = "Principal components", y = "Variance percentage") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black", size = 0),
    axis.text.x = element_text(face = "bold", color = "#000000"),
    axis.text.y = element_text(face = "bold", color = "#000000")
  )

# Save scree plots
# ------------------------------------------------------------------------------
# ggsave(plot = scree_plot_pca_1, "./results/05_PCA_plots/scree_plot_pca_1.png", width = 6, height = 5, dpi = 300)
# ggsave(plot = scree_plot_pca_2, "./results/05_PCA_plots/scree_plot_pca_2.png", width = 6, height = 5, dpi = 300)
# ggsave(plot = scree_plot_pca_3, "./results/05_PCA_plots/scree_plot_pca_3.png", width = 6, height = 5, dpi = 300)
ggsave(plot = scree_plot_pca_4, "./results/05_PCA_plots/scree_plot_pca_4.png", width = 10, height = 4, dpi = 300)


# PC projections using the two groups to color
# ------------------------------------------------------------------------------
# Check results
# pca_4 %>% tidy("samples")

# Augment PCA
# aug_pca_1 <- pca_1 %>% augment(encoded_data_set_1)
# aug_pca_2 <- pca_2 %>% augment(encoded_data_set_2)
# aug_pca_3 <- pca_3 %>% augment(encoded_data_set_3)
aug_pca_4 <- pca_4 %>% augment(data_set_4)

biplot_pca_4 <- aug_pca_4 %>%
  ggplot() +
  geom_point(aes(x = .fittedPC1, y = .fittedPC2, colour = activity_group)) +
  theme_linedraw() +
  labs(x = "PC1", y = "PC2") +
  ggtitle("Projection into the first two principal components") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank(),
    axis.line = element_line(colour = "black", size = 0),
    axis.text.x = element_text(face = "bold", color = "#000000"),
    axis.text.y = element_text(face = "bold", color = "#000000")
  )

# Save projection plots from the PCA
# ------------------------------------------------------------------------------
# ggsave(plot = biplot_pca_1, "./results/05_PCA_plots/biplot_pca_1.png", width = 7, height = 6, dpi = 300)
# ggsave(plot = biplot_pca_2, "./results/05_PCA_plots/biplot_pca_2.png", width = 7, height = 6, dpi = 300)
# ggsave(plot = biplot_pca_3, "./results/05_PCA_plots/biplot_pca_3.png", width = 7, height = 6, dpi = 300)
ggsave(plot = biplot_pca_4, "./results/05_PCA_plots/biplot_pca_4.png", width = 7, height = 6, dpi = 300)
