# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")

# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
#my_data <- read_tsv(file = "data/01_my_data.tsv")

# Wrangle data
# ------------------------------------------------------------------------------
#my_data_clean <- my_data # %>% ...

# Write data
# ------------------------------------------------------------------------------
#write_tsv(x = my_data_clean,
#          path = "data/02_my_data_clean.tsv")
data_set_1 <- read_tsv(file = "./data/01_load_data_set_1.tsv")
data_set_2 <- read_tsv(file = "./data/01_load_data_set_2.tsv")

# Wrangle data
# ------------------------------------------------------------------------------
data_set_1_clean <- data_set_1  %>% 
  select(Variant_ID, E3_score)

data_set_2_clean <- data_set_2  %>% 
  select(ERK2_Mutant, SCH_Average)

# Write data
# ------------------------------------------------------------------------------
write_tsv(x = data_set_1_clean,
          path = "./data/02_clean_data_set_1.tsv")

write_tsv(x = data_set_2_clean,
          path = "./data/02_clean_data_set_2.tsv")