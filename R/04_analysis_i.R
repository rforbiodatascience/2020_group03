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
# Create new dataframe with ID, Seq, E3 score only
data_set_1 <- data_set_1[, c(1, 2, 4, 5, 6)]
data_set_2 <- data_set_2[, c(1, 2, 4, 5, 6)]
data_set_3 <- data_set_3[, c(1, 2, 4, 5, 6)]
data_set_4 <- data_set_4[, c(1, 2, 4, 5, 6)]

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
# Heat Map
# pi
data_set_4 %>% ggplot() +
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


# Distribution
# applying filters
data_set_4 %>%
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


