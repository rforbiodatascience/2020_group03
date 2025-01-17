# 01_load.R
# ------------------------------------------------------------------------------
print('01_load.R -> loading data')


# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library("readxl")
library("ggplot2")
library('UniprotR') 

# Define functions
# ------------------------------------------------------------------------------
source(file = "./R/99_project_functions.R")


# Load data
# ------------------------------------------------------------------------------
# paper reference here: https://doi.org/10.1534/genetics.115.175802
data_set_1 <- read_excel('./data/_raw/genetics.115.175802-6.xls')

# paper reference here: https://doi.org/10.1016/j.celrep.2016.09.061
data_set_2 <- read_excel('./data/_raw/1-s2.0-S2211124716313171-mmc2.xlsx',
                         sheet = 'Supplemental_Table_1')

# paper reference here: https://www.mavedb.org/scoreset/urn:mavedb:00000036-a-1/
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
