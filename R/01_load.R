<<<<<<< HEAD
# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")

# Define functions
# ------------------------------------------------------------------------------
#source(file = "R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
#my_data_raw <- read_tsv(file = "data/_raw/my_raw_data.tsv")

# Wrangle data
# ------------------------------------------------------------------------------
#my_data <- my_data_raw # %>% ...

# Write data
# ------------------------------------------------------------------------------
#write_tsv(x = my_data,
#          path = "data/01_my_data.tsv")

library("ggseqlogo")
library("readxl")
library("ggplot2")
library('UniprotR') 

# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
# load amino acid scales
aa_scale_z5 <- load_scale('z5')

aa_scale_bl62 <- load_scale('bl62')

# load peptide sequence and results
data_set_1 <- read_excel('./data/_raw/genetics.115.175802-6.xls')

data_set_2 <- read_excel('./data/_raw/1-s2.0-S2211124716313171-mmc2.xlsx',
                         sheet = 'Supplemental_Table_1')


# Write data
# ------------------------------------------------------------------------------
write_tsv(x = data_set_1,
          path = "./_raw/data/01_load_data_set_1.tsv")

write_tsv(x = data_set_2,
          path = "./_raw/data/01_load_data_set_2.tsv")


=======
# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
#install.packages("tidyverse")
#install.packages("ggplot2")
#install.packages("UniprotR")
library("tidyverse")

# Define functions
# ------------------------------------------------------------------------------
#source(file = "R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
#my_data_raw <- read_tsv(file = "data/_raw/my_raw_data.tsv")

# Wrangle data
# ------------------------------------------------------------------------------
#my_data <- my_data_raw # %>% ...

# Write data
# ------------------------------------------------------------------------------
#write_tsv(x = my_data,
#          path = "data/01_my_data.tsv")

#library("ggseqlogo")
library("readxl")
#library("ggplot2")
#library('UniprotR') 

# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
# load amino acid scales
#aa_scale_z5 <- load_scale('z5')

#aa_scale_bl62 <- load_scale('bl62')

# load peptide sequence and results

# paper reference here: xxxx
data_set_1 <- read_excel('./data/_raw/genetics.115.175802-6.xls')

# paper reference here: xxxx
data_set_2 <- read_excel('./data/_raw/1-s2.0-S2211124716313171-mmc2.xlsx',
                         sheet = 'Supplemental_Table_1')

# paper reference here: xxxx
data_set_3 <- read_csv('./data/_raw/urn_mavedb_00000036-a-1_scores.csv', skip = 4)

# data_set_4: https://doi.org/10.1371/journal.pgen.1004918
data_set_4 <- read_tsv('./data/_raw/S1_Table.txt')


# Write data
# ------------------------------------------------------------------------------
write_tsv(x = data_set_1,
          path = "./data/01_load_data_set_1.tsv")

write_tsv(x = data_set_2,
          path = "./data/01_load_data_set_2.tsv")

write_tsv(x = data_set_3,
          path = "./data/01_load_data_set_3.tsv")

write_tsv(x = data_set_4,
          path = "./data/01_load_data_set_4.tsv")
>>>>>>> master
