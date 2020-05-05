<<<<<<< HEAD
# Install packages
# ------------------------------------------------------------------------------
#install.packages("tidyverse")
#install.packages("ggplot2")
#install.packages("UniprotR")
install.packages("devtools")
#install.packages("keras")
install.packages('neuralnet')
library(keras)
install_keras(tensorflow = "1.13.1")
# Run all scripts
# ------------------------------------------------------------------------------
source(file = "R/01_load.R")
source(file = "R/02_clean.R")
source(file = "R/03_augment.R")
source(file = "R/04_analysis_i.R")

=======
# Install packages
# ------------------------------------------------------------------------------
#install.packages("tidyverse")
#install.packages("ggplot2")
#install.packages("UniprotR")

# Run all scripts
# ------------------------------------------------------------------------------
source(file = "R/01_load.R")
source(file = "R/02_clean.R")
source(file = "R/03_augment.R")
#source(file = "R/04_analysis_i.")

>>>>>>> 288df4c889b2ff9ca2c8d8c0ed71ba630dee6161
