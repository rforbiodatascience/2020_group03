# 03_augment.R
# ------------------------------------------------------------------------------
print('03_augment.R -> augmenting data')


# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library("UniprotR")

# Define functions
# ------------------------------------------------------------------------------
source(file = "./R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
clean_data_set_1 <- read_tsv(file = "./data/02_clean_data_set_1.tsv")
clean_data_set_2 <- read_tsv(file = "./data/02_clean_data_set_2.tsv")
clean_data_set_3 <- read_tsv(file = "./data/02_clean_data_set_3.tsv")
clean_data_set_4 <- read_tsv(file = "./data/02_clean_data_set_4.tsv")

# Wrangle data
# ------------------------------------------------------------------------------
# little extra wrangling step
aug_data_set_1 <- clean_data_set_1 %>%
  mutate(variant = case_when(
    variant == "NA-NA" ~ "D2D",
    !variant == "NA-NA" ~ variant
  ))

# Augmenting with sequence
# ------------------------------------------------------------------------------
aug_data_set_1 <- convert_variant_to_sequence(aug_data_set_1, "P38398")
aug_data_set_2 <- convert_variant_to_sequence(clean_data_set_2, "P28482")
aug_data_set_3 <- convert_variant_to_sequence(clean_data_set_3, "Q5SW96")
aug_data_set_4 <- convert_variant_to_sequence(clean_data_set_4, "P04147")

# Write data
# ------------------------------------------------------------------------------
<<<<<<< HEAD
write_tsv(
  x = aug_data_set_1,
  path = "./data/03_aug_data_set_1.tsv"
)

write_tsv(
  x = aug_data_set_2,
  path = "./data/03_aug_data_set_2.tsv"
)

write_tsv(
  x = aug_data_set_3,
  path = "./data/03_aug_data_set_3.tsv"
)

write_tsv(
  x = aug_data_set_4,
  path = "./data/03_aug_data_set_4.tsv"
)
=======
write_tsv(x = aug_data_set_1, path = "./data/03_aug_data_set_1.tsv")
write_tsv(x = aug_data_set_2, path = "./data/03_aug_data_set_2.tsv")
write_tsv(x = aug_data_set_3, path = "./data/03_aug_data_set_3.tsv")
write_tsv(x = aug_data_set_4, path = "./data/03_aug_data_set_4.tsv")
>>>>>>> 1d23e18731429af2f697b28b94cbdfecf178f800
