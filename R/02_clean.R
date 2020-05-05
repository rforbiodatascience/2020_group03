# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library("dplyr")

# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
data_set_1 <- read_tsv(file = "./data/01_load_data_set_1.tsv")
data_set_2 <- read_tsv(file = "./data/01_load_data_set_2.tsv")
data_set_3 <- read_tsv(file = "./data/01_load_data_set_3.tsv")
data_set_4 <- read_tsv(file = "./data/01_load_data_set_4.tsv")

# Wrangle data
# ------------------------------------------------------------------------------

# select variant_ID and activity to predict, and rename columns
data_set_1_clean <- data_set_1  %>% 
  select(Variant_ID, E3_score) %>%
  rename(variant = Variant_ID,
         score = E3_score)
  
# select ERK2_Mutant and activity to predict, and rename columns
data_set_2_clean <- data_set_2  %>% 
  select(ERK2_Mutant, SCH_Average) %>%
  rename(variant = ERK2_Mutant,
         score = SCH_Average)

data_set_3_clean <- data_set_3  %>% 
  select(hgvs_pro, score) %>%
  rename(variant = hgvs_pro)

data_set_4_clean <- data_set_4 %>%
  pivot_longer(-X1, names_to = "substitution", values_to = "score") %>%
  unite("variant", X1:substitution, sep="")
  

# Write data
# ------------------------------------------------------------------------------
write_tsv(x = data_set_1_clean,
          path = "./data/02_clean_data_set_1.tsv")

write_tsv(x = data_set_2_clean,
          path = "./data/02_clean_data_set_2.tsv")

write_tsv(x = data_set_3_clean,
          path = "./data/02_clean_data_set_3.tsv")

write_tsv(x = data_set_4_clean,
          path = "./data/02_clean_data_set_4.tsv")