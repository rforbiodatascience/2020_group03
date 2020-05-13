# 00_do_it.R
# ------------------------------------------------------------------------------
# Install packages ONLY if needed
# ------------------------------------------------------------------------------
# install.packages("tidyverse")
# install.packages("ggplot2")
# install.packages("UniprotR")
# install.packages("devtools")
# install.packages('neuralnet')
# install.packages('glmnet')
# install.packages('caret')
# install.packages('neuralnet')
# install.packages('yardstick')
#  install.packages('ANN2')
# library("devtools")
# install_github("rstudio/keras")
# 
# Would you like to install miniconda? Y
library(keras)
install_keras(tensorflow = "1.13.1")


# Run all scripts
# ------------------------------------------------------------------------------
source(file = "./R/01_load.R")
source(file = "./R/02_clean.R")
source(file = "./R/02_get_unknowns.R")
source(file = "./R/03_augment.R")
source(file = "./R/03_get_unknowns.R")
source(file = "./R/04_data_visualization.R")
source(file = "./R/05_PCA.R") 
source(file = "./R/06_modelling_ANN_keras.R") # be aware lengthy process
source(file = "./R/06_modelling_ANN2.R")      # be aware lengthy process
source(file = "./R/06_modelling_glmnet.R")    # be aware lengthy process
source(file = "./R/06_glmnet_CV.R")           # be aware lengthy process
source(file = "./R/06_get_unkowns.R") 
source(file = "./R/06_model_summary_both.R")  # be aware lengthy process
source(file = "./R/07_predict.R") 