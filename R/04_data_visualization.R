# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library("Peptides")
library("ggseqlogo")
library("dplyr")


# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
data_set_1 <- read_tsv(file = "./data/03_aug_data_set_1.tsv")
data_set_2 <- read_tsv(file = "./data/03_aug_data_set_2.tsv")
data_set_3 <- read_tsv(file = "./data/03_aug_data_set_3.tsv")
data_set_4 <- read_tsv(file = "./data/03_aug_data_set_4.tsv")

# Wrangle data
# ------------------------------------------------------------------------------
# Calculate molecular weight
data_set_1$MW <- data_set_1 %>%
  pull(sequence) %>%
  mw(monoisotopic = FALSE)
data_set_2$MW <- data_set_2 %>%
  pull(sequence) %>%
  mw(monoisotopic = FALSE)
data_set_3$MW <- data_set_3 %>%
  pull(sequence) %>%
  mw(monoisotopic = FALSE)
data_set_4$MW <- data_set_4 %>%
  pull(sequence) %>%
  mw(monoisotopic = FALSE)

# Calculate pI of the sequences
data_set_1$pI <- data_set_1 %>%
  pull(sequence) %>%
  pI(pKscale = "EMBOSS")
data_set_2$pI <- data_set_2 %>%
  pull(sequence) %>%
  pI(pKscale = "EMBOSS")
data_set_3$pI <- data_set_3 %>%
  pull(sequence) %>%
  pI(pKscale = "EMBOSS")
data_set_4$pI <- data_set_4 %>%
  pull(sequence) %>%
  pI(pKscale = "EMBOSS")

# Calculate z-scale of the sequences
data_set_1$charge <- data_set_1 %>%
  pull(sequence) %>%
  charge(pH = 7, pKscale = "EMBOSS")
data_set_2$charge <- data_set_2 %>%
  pull(sequence) %>%
  charge(pH = 7, pKscale = "EMBOSS")
data_set_3$charge <- data_set_3 %>%
  pull(sequence) %>%
  charge(pH = 7, pKscale = "EMBOSS")
data_set_4$charge <- data_set_4 %>%
  pull(sequence) %>%
  charge(pH = 7, pKscale = "EMBOSS")

# Calculate z-scale of the sequences
zscales_1 <- data_set_1 %>%
  pull(sequence) %>%
  zScales()
zscales_2 <- data_set_2 %>%
  pull(sequence) %>%
  zScales()
zscales_3 <- data_set_3 %>%
  pull(sequence) %>%
  zScales()
zscales_4 <- data_set_4 %>%
  pull(sequence) %>%
  zScales()

# separate Z-scales in columns
data_set_1 <- data_set_1 %>%
  zscale2col(zscales_1)
data_set_2 <- data_set_2 %>%
  zscale2col(zscales_2)
data_set_3 <- data_set_3 %>%
  zscale2col(zscales_3)
data_set_4 <- data_set_4 %>%
  zscale2col(zscales_4)

# Add aminoacid type / call function aminoacid_type
data_set_1 <- data_set_1 %>%
  mutate(aminoacid_class = aminoacid_type(mutation))
data_set_2 <- data_set_2 %>%
  mutate(aminoacid_class = aminoacid_type(mutation))
data_set_3 <- data_set_3 %>%
  mutate(aminoacid_class = aminoacid_type(mutation))
data_set_4 <- data_set_4 %>%
  mutate(aminoacid_class = aminoacid_type(mutation))

# Visualise data
# ------------------------------------------------------------------------------
# Generate pI heatmaps
# ------------------------------------------------------------------------------
heatmap_data_set_pI_1 <- heatmap(data_set_1, x = "mutation_position", y = "mutation", fill = "pI")
heatmap_data_set_pI_2 <- heatmap(data_set_2, x = "mutation_position", y = "mutation", fill = "pI")
heatmap_data_set_pI_3 <- heatmap(data_set_3, x = "mutation_position", y = "mutation", fill = "pI")
heatmap_data_set_pI_4 <- heatmap(data_set_4, x = "mutation_position", y = "mutation", fill = "pI")

# Save heatmaps
# ------------------------------------------------------------------------------
ggsave(plot = heatmap_data_set_pI_1, "./results/04_heatmaps/heatmap_data_set_pI_1.png", width = 10, height = 3, dpi = 300)
ggsave(plot = heatmap_data_set_pI_2, "./results/04_heatmaps/heatmap_data_set_pI_2.png", width = 10, height = 3, dpi = 300)
ggsave(plot = heatmap_data_set_pI_3, "./results/04_heatmaps/heatmap_data_set_pI_3.png", width = 10, height = 3, dpi = 300)
ggsave(plot = heatmap_data_set_pI_4, "./results/04_heatmaps/heatmap_data_set_pI_4.png", width = 10, height = 3, dpi = 300)

# Generate z-scale 5 heatmaps
# ------------------------------------------------------------------------------
heatmap_data_set_Z5_1 <- heatmap(data_set_1, x = "mutation_position", y = "mutation", fill = "Z5")
heatmap_data_set_Z5_2 <- heatmap(data_set_2, x = "mutation_position", y = "mutation", fill = "Z5")
heatmap_data_set_Z5_3 <- heatmap(data_set_3, x = "mutation_position", y = "mutation", fill = "Z5")
heatmap_data_set_Z5_4 <- heatmap(data_set_4, x = "mutation_position", y = "mutation", fill = "Z5")

# Save heatmaps
# ------------------------------------------------------------------------------
ggsave(plot = heatmap_data_set_Z5_1, "./results/04_heatmaps/heatmap_data_set_Z5_1.png", width = 10, height = 3, dpi = 300)
ggsave(plot = heatmap_data_set_Z5_2, "./results/04_heatmaps/heatmap_data_set_Z5_2.png", width = 10, height = 3, dpi = 300)
ggsave(plot = heatmap_data_set_Z5_3, "./results/04_heatmaps/heatmap_data_set_Z5_3.png", width = 10, height = 3, dpi = 300)
ggsave(plot = heatmap_data_set_Z5_4, "./results/04_heatmaps/heatmap_data_set_Z5_4.png", width = 10, height = 3, dpi = 300)


# Generate distribution scores plots
# ------------------------------------------------------------------------------
distribution_data_set_1 <- score_distribution(data_set_1, threshold = 1.5)
distribution_data_set_2 <- score_distribution(data_set_2, threshold = 11)
distribution_data_set_3 <- score_distribution(data_set_3, threshold = 0.975)
distribution_data_set_4 <- score_distribution(data_set_4, threshold = 1.5)

# Save distribution scores plots
# ------------------------------------------------------------------------------
ggsave(plot = distribution_data_set_1, "./results/04_distribution_plots/distribution_data_set_1.png", width = 7, height = 10, dpi = 300)
ggsave(plot = distribution_data_set_2, "./results/04_distribution_plots/distribution_data_set_2.png", width = 7, height = 10, dpi = 300)
ggsave(plot = distribution_data_set_3, "./results/04_distribution_plots/distribution_data_set_3.png", width = 7, height = 10, dpi = 300)
ggsave(plot = distribution_data_set_4, "./results/04_distribution_plots/distribution_data_set_4.png", width = 7, height = 10, dpi = 300)

# Score heatmaps. With a cuff-off scores
# ------------------------------------------------------------------------------
heatmap_data_set_score_1 <- heatmap_score(data_set_1, score_min = 0, score_max = 1.5)
heatmap_data_set_score_2 <- heatmap_score(data_set_2, score_min = 1, score_max = 10)
heatmap_data_set_score_3 <- heatmap_score(data_set_3, score_min = 0, score_max = 1)
heatmap_data_set_score_4 <- heatmap_score(data_set_4, score_min = -2, score_max = 2)

# Save heatmaps
# ------------------------------------------------------------------------------
ggsave(plot = heatmap_data_set_score_1, "./results/04_heatmaps/heatmap_data_set_score_1.png", width = 10, height = 3, dpi = 300)
ggsave(plot = heatmap_data_set_score_2, "./results/04_heatmaps/heatmap_data_set_score_2.png", width = 10, height = 3, dpi = 300)
ggsave(plot = heatmap_data_set_score_3, "./results/04_heatmaps/heatmap_data_set_score_3.png", width = 10, height = 3, dpi = 300)
ggsave(plot = heatmap_data_set_score_4, "./results/04_heatmaps/heatmap_data_set_score_4.png", width = 10, height = 3, dpi = 300)


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
  heatmap_g_score(score_min = 1, score_max = 10)
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
ggsave(plot = heatmap_data_set_score_g_1, "./results/04_heatmaps/heatmap_data_set_score_g_1.png", width = 10, height = 3, dpi = 300)
ggsave(plot = heatmap_data_set_score_g_2, "./results/04_heatmaps/heatmap_data_set_score_g_2.png", width = 10, height = 3, dpi = 300)
ggsave(plot = heatmap_data_set_score_g_3, "./results/04_heatmaps/heatmap_data_set_score_g_3.png", width = 10, height = 3, dpi = 300)
ggsave(plot = heatmap_data_set_score_g_4, "./results/04_heatmaps/heatmap_data_set_score_g_4.png", width = 10, height = 3, dpi = 300)


# Generate density plots
# ------------------------------------------------------------------------------
density_data_set_1 <- density_plot(data_set_1, title = "Score density plot for Data Set 1")
density_data_set_2 <- density_plot(data_set_2, title = "Score density plot for Data Set 2")
density_data_set_3 <- density_plot(data_set_3, title = "Score density plot for Data Set 3")
density_data_set_4 <- density_plot(data_set_4, title = "Score density plot for Data Set 4")

# Generate density plot per mutation
# ------------------------------------------------------------------------------
density_residue_data_set_1 <- density_plot_residue(data_set_1)
density_residue_data_set_2 <- density_plot_residue(data_set_2)
density_residue_data_set_3 <- density_plot_residue(data_set_3)
density_residue_data_set_4 <- density_plot_residue(data_set_4)

# Save density plots
# ------------------------------------------------------------------------------
ggsave(plot = density_data_set_1, "./results/04_density_plots/density_data_set_1.png", width = 7, height = 5, dpi = 300)
ggsave(plot = density_data_set_2, "./results/04_density_plots/density_data_set_2.png", width = 7, height = 5, dpi = 300)
ggsave(plot = density_data_set_3, "./results/04_density_plots/density_data_set_3.png", width = 7, height = 5, dpi = 300)
ggsave(plot = density_data_set_4, "./results/04_density_plots/density_data_set_4.png", width = 7, height = 5, dpi = 300)

# Save density plots per mutation
# ------------------------------------------------------------------------------
ggsave(plot = density_residue_data_set_1, "./results/04_density_plots/density_residue_data_set_1.png", width = 7, height = 6.25, dpi = 300)
ggsave(plot = density_residue_data_set_2, "./results/04_density_plots/density_residue_data_set_2.png", width = 7, height = 6.25, dpi = 300)
ggsave(plot = density_residue_data_set_3, "./results/04_density_plots/density_residue_data_set_3.png", width = 7, height = 6.25, dpi = 300)
ggsave(plot = density_residue_data_set_4, "./results/04_density_plots/density_residue_data_set_4.png", width = 7, height = 6.25, dpi = 300)



# Make sequence logo plots
# ------------------------------------------------------------------------------
seq_logo_data_set_1 <- sequence_logo(data_set_1, start_position = 40, end_position = 80)
seq_logo_data_set_2 <- sequence_logo(data_set_2, start_position = 40, end_position = 80)
seq_logo_data_set_3 <- sequence_logo(data_set_3, start_position = 110, end_position = 150)
seq_logo_data_set_4 <- sequence_logo(data_set_4, start_position = 135, end_position = 185)

# Save sequence logos
# ------------------------------------------------------------------------------
ggsave(plot = seq_logo_data_set_1, "./results/04_sequence_logos/seq_logo_data_set_1.png", width = 10, height = 3, dpi = 300)
ggsave(plot = seq_logo_data_set_2, "./results/04_sequence_logos/seq_logo_data_set_2.png", width = 10, height = 3, dpi = 300)
ggsave(plot = seq_logo_data_set_3, "./results/04_sequence_logos/seq_logo_data_set_3.png", width = 10, height = 3, dpi = 300)
ggsave(plot = seq_logo_data_set_4, "./results/04_sequence_logos/seq_logo_data_set_4.png", width = 10, height = 3, dpi = 300)



# Drop NA for the next plot
# ------------------------------------------------------------------------------
data_set_4 <- data_set_4 %>%
  drop_na(score)

data_set_3 <- data_set_3 %>%
  drop_na(score)

# Plot effective amino acids at each position with a cutoff
# ------------------------------------------------------------------------------
quick_sar_data_set_1 <- quick_sar(data_set_1, cutoff = 0.5)
quick_sar_data_set_2 <- quick_sar(data_set_2, cutoff = 3)
quick_sar_data_set_3 <- quick_sar(data_set_3, cutoff = 0.6)
quick_sar_data_set_4 <- quick_sar(data_set_4, cutoff = 0)

# Save quickSAR plots
# ------------------------------------------------------------------------------
ggsave(plot = quick_sar_data_set_1, "./results/04_quick_SAR/quick_sar_data_set_1.png", width = 10, height = 3, dpi = 300)
ggsave(plot = quick_sar_data_set_2, "./results/04_quick_SAR/quick_sar_data_set_2.png", width = 10, height = 3, dpi = 300)
ggsave(plot = quick_sar_data_set_3, "./results/04_quick_SAR/quick_sar_data_set_3.png", width = 10, height = 3, dpi = 300)
ggsave(plot = quick_sar_data_set_4, "./results/04_quick_SAR/quick_sar_data_set_4.png", width = 10, height = 3, dpi = 300)