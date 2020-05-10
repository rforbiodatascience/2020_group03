# 00_do_it.R
# ------------------------------------------------------------------------------


# Install packages ONLY if needed
# ------------------------------------------------------------------------------
# install.packages("tidyverse")
# install.packages("ggplot2")
# install.packages("UniprotR")
# install.packages("devtools")
# install.packages('neuralnet')
# library("devtools")
# install_github("rstudio/keras")
# 
# Would you like to install miniconda? Y
# library(keras)
# install_keras(tensorflow = "1.13.1")


# Run all scripts
# ------------------------------------------------------------------------------
source(file = "./R/01_load.R")
source(file = "./R/02_clean.R")
source(file = "./R/03_augment.R")
source(file = "./R/04_heatmaps.R")
source(file = "./R/04_quicksar.R")
source(file = "./R/04_density_plots.R")