<<<<<<< HEAD
# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")

# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_functions.R")

# Load data
# ------------------------------------------------------------------------------
my_data_clean_aug <- read_table2(file = "data/clean_test.txt", col_names = FALSE)

# Wrangle data
# ------------------------------------------------------------------------------

# Remove non-important columns
my_data_clean_aug <- my_data_clean_aug %>% 
  rename(
    peptide = X1,
    activity = X2,
  ) %>%  select(c(1,2))

# Add new columns (encoding)

pepts <- my_data_clean_aug %>% 
  pull(peptide) %>% 
  encode_peptide("blosum62")

# Write data
# ------------------------------------------------------------------------------
#write_tsv(x = my_data_clean_aug,
#         path = "data/03_my_data_clean_aug.tsv")
library("UniprotR")

# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data
# ------------------------------------------------------------------------------
clean_data_set_1 <- read_tsv(file = "./_raw/data/02_clean_data_set_1.tsv")
clean_data_set_2 <- read_tsv(file = "./_raw/data/02_clean_data_set_2.tsv")


# Wrangle data
# ------------------------------------------------------------------------------
aug_data_set_1 <- clean_data_set_1 %>% 
  mutate(
    Mutated_residue = str_extract(Variant_ID, "[A-Z]"),
    Mutation_position = str_extract(Variant_ID, "[0-9]+"),
    Mutation = str_extract(Variant_ID, "[A-Z,*]$"),
    Mutation = str_replace(Mutation, "[*]", Mutated_residue),
    Sequence = make_sequence(Mutated_residue, Mutation_position, Mutation, "P38398"))

# aug_data_set_2 <- clean_data_set_2 %>% 
#   mutate(
#     Mutated_residue = str_extract(ERK2_Mutant, "[A-Z]"),
#     Mutation_position = str_extract(ERK2_Mutant, "[0-9]+"),
#     Mutation = str_extract(ERK2_Mutant, "[A-Z,*]$"),
#     Mutation = str_replace(Mutation, "[*]", Mutated_residue),
#     Sequence = make_sequence(Mutated_residue, Mutation_position, Mutation, "P38398"))


# Write data
# ------------------------------------------------------------------------------
write_tsv(x = aug_data_set_1,
          path = "./_raw/data/03_aug_data_set_1.tsv")

=======
# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library("UniprotR")
library("protr")

# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_functions.R")

# Load data
# ------------------------------------------------------------------------------
data <- read_table2(file = "data/clean_test.txt", col_names = FALSE)

# Wrangle data
# ------------------------------------------------------------------------------

# Remove non-important columns
<<<<<<< HEAD
# my_data_clean_aug <- my_data_clean_aug %>% 
#   rename(
#     peptide = X1,
#     activity = X2,
#   ) %>%  select(c(1,2))
=======
data <- data %>% 
  rename(
    peptide = X1,
    activity = X2,
  ) %>%  select(c(1,2))
>>>>>>> master

data$len = str_length(data$peptide)
data <- data %>% group_by(len) %>% filter(len == 15) 

<<<<<<< HEAD
# Add new columns (encoding)
encoded_seq = data %>%
  pull(peptide) %>% 
  encode_peptide(m = "blosum62")
=======
#pepts <- my_data_clean_aug %>% 
#  pull(peptide) %>% 
#  encode_peptide("blosum62")
>>>>>>> master

# Write data
# ------------------------------------------------------------------------------
#write_tsv(x = my_data_clean_aug,
<<<<<<< HEAD
#         path = "data/03_my_data_clean_aug.tsv")
=======
#         path = "data/03_my_data_clean_aug.tsv")


# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data
# ------------------------------------------------------------------------------
clean_data_set_1 <- read_tsv(file = "./data/02_clean_data_set_1.tsv")
clean_data_set_2 <- read_tsv(file = "./data/02_clean_data_set_2.tsv")
clean_data_set_3 <- read_tsv(file = "./data/02_clean_data_set_3.tsv")
clean_data_set_4 <- read_tsv(file = "./data/02_clean_data_set_4.tsv")


# Wrangle data
# ------------------------------------------------------------------------------
aug_data_set_1 <- clean_data_set_1 %>%
  mutate(variant = case_when(variant=="NA-NA" ~ "D2D",
                             !variant=="NA-NA" ~ variant))
# augmenting with sequence
aug_data_set_1 <- convert_variant_to_sequence(aug_data_set_1, "P38398")
aug_data_set_2 <- convert_variant_to_sequence(clean_data_set_2, "P28482")
aug_data_set_3 <- convert_variant_to_sequence(clean_data_set_3, "Q5SW96")
aug_data_set_4 <- convert_variant_to_sequence(clean_data_set_3, "P04147")

# Write data
# ------------------------------------------------------------------------------
write_tsv(x = aug_data_set_1,
          path = "./data/03_aug_data_set_1.tsv")

<<<<<<< HEAD
write_tsv(x = aug_data_set_2,
          path = "./data/03_aug_data_set_2.tsv")

write_tsv(x = aug_data_set_3,
          path = "./data/03_aug_data_set_3.tsv")

write_tsv(x = aug_data_set_4,
          path = "./data/03_aug_data_set_4.tsv")
=======
>>>>>>> master
>>>>>>> master
>>>>>>> master
