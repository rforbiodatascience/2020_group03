# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())


# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library("caret")
library("UniprotR")
library("ANN2")
library("yardstick")
library("glmnet")

# Load unknown data
# ------------------------------------------------------------------------------
unknown_data_set_4 <- read_tsv(file = "./data/03_aug_unknown_data_set_4.tsv")
known_data_set_4 <- read_tsv(file = "./data/03_aug_data_set_4.tsv")


# Define functions
# ------------------------------------------------------------------------------
source(file = "./R/99_project_functions.R")


# Load model
# ------------------------------------------------------------------------------
load(file = "./results/06_glmnet_data/06_glmnet_z_scales_data_set_4.RData")

data = unknown_data_set_4

# Get sequence window and truncate
# ------------------------------------------------------------------------------
sequence_window <- data %>%
  mutate(mutation_position = as.integer(mutation_position)) %>%
  pull(mutation_position)

data <- data %>%
  mutate(sequence_to_model = str_sub(sequence, min(sequence_window), max(sequence_window)))

# Encode sequences and make test and train sets
# ------------------------------------------------------------------------------
to_predict <- data %>%
  pull(sequence_to_model) %>%
  encode_peptide(m = "z_scales")

to_predict <- unname(as.matrix(to_predict))

prediction <- predict(fit.elasticnet, newx = to_predict, s = 0.01, alpha = 0.2)

knowns <- known_data_set_4 %>%
  select(score) %>%
  add_column(status = "known")

prediction <- tibble(prediction) %>%
  add_column(status = "unknown") %>%
  rename(score = prediction)


comp <- prediction %>%
  full_join(knowns, by="score") %>%
  mutate(status = case_when(
    is.na(status.x) ~ status.y,
    is.na(status.y) ~ status.x))

density <- ggplot(comp, aes(score, color=status, fill=status)) +
  geom_density(alpha=0.4) +
  theme_classic() +
  labs(x="Score", y="Density",title="Score density plot for known and unknown in data set 4") +
  theme(
    legend.position = "bottom",
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black", size = 0),
    axis.text.x = element_text(face = "bold", color = "#000000"),
    axis.text.y = element_text(face = "bold", color = "#000000")
  )

ggsave(plot = density, "./results/07_predict/07_density_data_set_4.png", width = 12, height = 5, dpi=300)