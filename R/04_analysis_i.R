# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library("Peptides")
library("ggseqlogo")
library("dplyr")
library("hrbrthemes")
library("RColorBrewer")
library("corrr")

# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
aug_clean_data_set_1 <- read_tsv(file = "./data/03_aug_data_set_1.tsv")
aug_clean_data_set_2 <- read_tsv(file = "./data/03_aug_data_set_2.tsv")
aug_clean_data_set_3 <- read_tsv(file = "./data/03_aug_data_set_3.tsv")
aug_clean_data_set_4 <- read_tsv(file = "./data/03_aug_data_set_4.tsv")

# Wrangle data
# ------------------------------------------------------------------------------
# Create new dataframe with ID, Seq, E3 score only
seq_properties <- aug_clean_data_set_1[, c(1, 2, 4, 5, 6)]

# Calculate molecular weight
seq_properties$MW <- aug_clean_data_set_1 %>%
  pull(sequence) %>%
  mw(monoisotopic = FALSE)

# Calculate pI of the sequences
seq_properties$pI <- aug_clean_data_set_1 %>%
  pull(sequence) %>%
  pI(pKscale = "EMBOSS")

# Calculate z-scale of the sequences
seq_properties$charge <- aug_clean_data_set_1 %>%
  pull(sequence) %>%
  charge(pH = 7, pKscale = "EMBOSS")

# Calculate z-scale of the sequences
zscales <- aug_clean_data_set_1 %>%
  pull(sequence) %>%
  zScales()

zscales_df <- data.frame(zscales)
zscales <- as_tibble(cbind(nms = names(zscales_df), t(zscales_df))) %>%
  select(.data$Z1, .data$Z2, .data$Z3, .data$Z4, .data$Z5)
seq_properties <- bind_cols(seq_properties, zscales)
seq_properties <- seq_properties %>%
  mutate(
    Z1 = as.numeric(Z1),
    Z2 = as.numeric(Z2),
    Z3 = as.numeric(Z3),
    Z4 = as.numeric(Z4),
    Z5 = as.numeric(Z5)
  )

# Add aminoacid type
seq_properties <- seq_properties %>%
  mutate(aminoacid_class = case_when(
    mutation == "R" | mutation == "H" | mutation == "K" ~ "basic",
    mutation == "D" | mutation == "E" ~ "acid",
    mutation == "S" | mutation == "T" ~ "hydroxilic",
    mutation == "Q" | mutation == "N" ~ "amidic",
    mutation == "A" | mutation == "G" | mutation == "I" | mutation == "L" | mutation == "P" | mutation == "V" ~ "aliphatic",
    mutation == "F" | mutation == "W" | mutation == "Y" ~ "aromatic",
    mutation == "C" | mutation == "M" ~ "sulfur-containing"
  ))

# Visualise data
# ------------------------------------------------------------------------------
## Logo we need shorter sequences

# Visualize mutation pI per mutation and per position
seq_properties %>%
  group_by(mutation) %>%
  ggplot() +
  geom_point(aes(x = mutation, y = pI, color = aminoacid_class)) +
  theme_classic() +
  labs(x = "Residue mutated", y = "pI", color = "Residue character") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black", size = 0),
    axis.text.x = element_text(face = "bold", color = "#000000"),
    axis.text.y = element_text(face = "bold", color = "#000000")
  )

seq_properties %>%
  group_by(mutation_position) %>%
  ggplot() +
  geom_point(aes(x = mutation_position, y = pI, color = aminoacid_class)) +
  theme_classic() +
  labs(x = "Residue mutated", y = "pI", color = "Residue character") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black", size = 0),
    axis.text.x = element_text(face = "bold", color = "#000000"),
    axis.text.y = element_text(face = "bold", color = "#000000")
  )


# Visualize mutation mW per mutation and per position
seq_properties %>%
  group_by(mutation) %>%
  ggplot() +
  geom_boxplot(aes(x = mutation, y = MW, fill = aminoacid_class), color = "black", width = 0.5, position = position_dodge(0.9)) +
  theme_classic() +
  labs(x = "Residue mutated", y = "mW", color = "Residue character") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black", size = 0),
    axis.text.x = element_text(face = "bold", color = "#000000"),
    axis.text.y = element_text(face = "bold", color = "#000000")
  )

seq_properties %>%
  group_by(mutation_position) %>%
  ggplot() +
  geom_point(aes(x = mutation_position, y = MW, color = aminoacid_class)) +
  theme_classic() +
  labs(x = "Residue mutated", y = "mW", color = "Residue character") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black", size = 0),
    axis.text.x = element_text(face = "bold", color = "#000000"),
    axis.text.y = element_text(face = "bold", color = "#000000")
  )

# Visualize charge per mutation and per position
seq_properties %>%
  group_by(mutation) %>%
  ggplot() +
  geom_point(aes(x = mutation, y = charge, color = aminoacid_class)) +
  theme_classic() +
  labs(x = "Residue mutated", y = "Net charge", color = "Residue character") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black", size = 0),
    axis.text.x = element_text(face = "bold", color = "#000000"),
    axis.text.y = element_text(face = "bold", color = "#000000")
  )

seq_properties %>%
  group_by(mutation_position) %>%
  ggplot() +
  geom_point(aes(x = mutation_position, y = charge, color = aminoacid_class)) +
  theme_classic() +
  labs(x = "Residue mutated", y = "Net Charge", color = "Residue character") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black", size = 0),
    axis.text.x = element_text(face = "bold", color = "#000000"),
    axis.text.y = element_text(face = "bold", color = "#000000")
  )

# Visualize Z5 per mutation and per position
seq_properties %>%
  group_by(mutation) %>%
  ggplot() +
  geom_point(aes(x = mutation, y = Z5, color = aminoacid_class)) +
  theme_classic() +
  labs(x = "Residue mutated", y = "Z5", color = "Residue character") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black", size = 0),
    axis.text.x = element_text(face = "bold", color = "#000000"),
    axis.text.y = element_text(face = "bold", color = "#000000")
  )

seq_properties %>%
  group_by(mutation_position) %>%
  ggplot() +
  geom_point(aes(x = mutation_position, y = Z5, color = aminoacid_class)) +
  theme_classic() +
  labs(x = "Residue mutated", y = "Z5", color = "Residue character") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black", size = 0),
    axis.text.x = element_text(face = "bold", color = "#000000"),
    axis.text.y = element_text(face = "bold", color = "#000000")
  )

# Heat Map
# pi
seq_properties %>% ggplot() +
  geom_tile(aes(y = mutation, x = mutation_position, fill = pI), color = "black") +
  scale_fill_distiller(palette = "RdPu") +
  theme_classic() +
  scale_x_discrete(breaks = c(0, 305)) +
  labs(x = "Mutation position", y = "Residue mutated", color = "Score") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black", size = 0),
    axis.text.x = element_text(face = "bold", color = "#000000"),
    axis.text.y = element_text(face = "bold", color = "#000000")
  )


# mw
seq_properties %>% ggplot() +
  geom_tile(aes(y = mutation, x = mutation_position, fill = MW), color = "black") +
  scale_fill_distiller(palette = "RdPu") +
  theme_classic() +
  scale_x_discrete(breaks = c(0, 305)) +
  labs(x = "Mutation position", y = "Residue mutated", color = "Score") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black", size = 0),
    axis.text.x = element_text(face = "bold", color = "#000000"),
    axis.text.y = element_text(face = "bold", color = "#000000")
  )


# z5
seq_properties %>% ggplot() +
  geom_tile(aes(y = mutation, x = mutation_position, fill = Z5), color = "black") +
  scale_fill_distiller(palette = "RdPu") +
  theme_classic() +
  scale_x_discrete(breaks = c(0, 305)) +
  labs(x = "Mutation position", y = "Residue mutated", color = "Score") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black", size = 0),
    axis.text.x = element_text(face = "bold", color = "#000000"),
    axis.text.y = element_text(face = "bold", color = "#000000")
  )


# z5 grouped by type of aa
seq_properties %>%
  group_by(aminoacid_class) %>%
  ggplot() +
  geom_tile(aes(y = aminoacid_class, x = mutation_position, fill = Z5), color = "black") +
  scale_fill_distiller(palette = "RdPu") +
  theme_classic() +
  scale_x_discrete(breaks = c(0, 305)) +
  labs(x = "Mutation position", y = "Residue mutated", color = "Score") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black", size = 0),
    axis.text.x = element_text(face = "bold", color = "#000000"),
    axis.text.y = element_text(face = "bold", color = "#000000")
  )


# Score
seq_properties %>% ggplot() +
  geom_tile(aes(y = mutation, x = mutation_position, fill = score), color = "black") +
  scale_fill_distiller(palette = "RdPu") +
  theme_classic() +
  scale_x_discrete(breaks = c(0, 305)) +
  labs(x = "Mutation position", y = "Residue mutated", color = "Score") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black", size = 0),
    axis.text.x = element_text(face = "bold", color = "#000000"),
    axis.text.y = element_text(face = "bold", color = "#000000")
  )


# Score grouped by type of aa
seq_properties %>%
  group_by(aminoacid_class) %>%
  ggplot() +
  geom_tile(aes(y = aminoacid_class, x = mutation_position, fill = score), color = "black") +
  scale_fill_distiller(palette = "RdPu") +
  theme_classic() +
  scale_x_discrete(breaks = c(0, 305)) +
  labs(x = "Mutation position", y = "Residue mutated", color = "Score") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black", size = 0),
    axis.text.x = element_text(face = "bold", color = "#000000"),
    axis.text.y = element_text(face = "bold", color = "#000000")
  )
# Distribution
# Get sd  and means

summary <- seq_properties %>%
  group_by(mutation_position) %>%
  summarise(mean = mean(score, na.rm = T), sd = sd(score, na.rm = T))

summary %>%
  ggplot(aes(x = mutation_position, y = mean)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
    width = .2,
    position = position_dodge(0.05)
  ) +
  theme_classic() +
  labs(x = "Residue mutated", y = "score") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black", size = 0),
    axis.text.x = element_text(face = "bold", color = "#000000"),
    axis.text.y = element_text(face = "bold", color = "#000000")
  )


colourCount <- 20
getPalette <- colorRampPalette(brewer.pal(12, "Accent"))

# Distribution
# applying filters
seq_properties %>%
  group_by(mutation) %>%
  filter(score > 1.5) %>%
  ggplot(aes(x = mutation_position, y = score, color = mutation)) +
  scale_color_manual(values = getPalette(colourCount)) +
  geom_point() +
  theme_classic() +
  labs(x = "Residue mutated", y = "score") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black", size = 0),
    axis.text.x = element_text(face = "bold", color = "#000000"),
    axis.text.y = element_text(face = "bold", color = "#000000")
  ) +
  facet_grid(aminoacid_class ~ .)

# Same but with heatmap
seq_properties %>%
  group_by(mutation) %>%
  filter(score > 1.5) %>%
  ggplot() +
  geom_tile(aes(y = aminoacid_class, x = mutation_position, fill = score), color = "black") +
  scale_fill_distiller(palette = "RdPu") +
  theme_classic() +
  scale_x_discrete(breaks = c(0, 305)) +
  labs(x = "Mutation position", y = "Residue mutated", color = "Score") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black", size = 0),
    axis.text.x = element_text(face = "bold", color = "#000000"),
    axis.text.y = element_text(face = "bold", color = "#000000")
  )


# Box plot
seq_properties %>%
  group_by(mutation) %>%
  ggplot(aes(y = score, color = mutation, fill = mutation)) +
  geom_boxplot(width = 0.5, position = position_dodge(0.9), alpha = 0.5) +
  scale_fill_manual(values = getPalette(colourCount)) +
  scale_color_manual(values = getPalette(colourCount)) +
  theme_classic() +
  labs(x = "Residue mutated", y = "score") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black", size = 0),
    axis.text.x = element_text(face = "bold", color = "#000000"),
    axis.text.y = element_text(face = "bold", color = "#000000"),
    legend.position = c(0.7, 0.223)
  ) +
  facet_wrap(. ~ aminoacid_class) +
  scale_x_discrete(breaks = c(0, 0)) +
  guides(fill = guide_legend(ncol = 5))


