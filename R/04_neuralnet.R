# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library("neuralnet")
library("caret")
library("GauPro")
library("UniprotR")
library("ANN2")


# Define functions
# ------------------------------------------------------------------------------
source(file = "./R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
#clean_data_set_1 <- read_tsv(file = "./data/02_clean_data_set_1.tsv")
#clean_data_set_2 <- read_tsv(file = "./data/02_clean_data_set_2.tsv")
clean_data_set_3 <- read_tsv(file = "./data/02_clean_data_set_3.tsv")
#clean_data_set_4 <- read_tsv(file = "./data/02_clean_data_set_4.tsv")

# Wrangle data
# ------------------------------------------------------------------------------
#aug_data_set_1 <- clean_data_set_1 %>%
#  mutate(variant = case_when(variant=="NA-NA" ~ "D2D",
#                             !variant=="NA-NA" ~ variant))

# augmenting with sequence
#aug_data_set_1 <- convert_variant_to_sequence(aug_data_set_1, "P38398")
#aug_data_set_2 <- convert_variant_to_sequence(clean_data_set_2, "P28482")
aug_data_set_3 <- convert_variant_to_sequence(clean_data_set_3, "Q5SW96")
#aug_data_set_4 <- convert_variant_to_sequence(clean_data_set_4, "P04147")

data <- aug_data_set_3 %>%
  drop_na(score)

#
# wrap in function
#

sequence_window <- data %>%
  mutate(mutation_position = as.integer(mutation_position)) %>%
  pull(mutation_position)

data <- data %>%
  mutate(sequence_to_model = str_sub(sequence, min(sequence_window), max(sequence_window)))

Y <- data  %>%
  select(score)
 
Y <- unname(as.matrix(Y))

set.seed(42)
partition <- createDataPartition(y = Y, p = 0.75, list = F)
training = data[partition, ]
test <- data[-partition, ]

#Add new columns (encoding)
X_train <- training %>%
  pull(sequence_to_model) %>%
  encode_peptide(m = "blosum62")

X_test <- test %>%
  pull(sequence_to_model) %>%
  encode_peptide(m = "blosum62")




Y_test <- test %>%
  select(score)

Y_train <- training %>%
  select(score)

Y_test <- unname(as.matrix(Y_test))
Y_train <- unname(as.matrix(Y_train))

X_train <- unname(as.matrix(X_train))
X_test <- unname(as.matrix(X_test))

print('NN running...')

NN <- neuralnetwork(X = X_train, y = Y_train, hidden.layers = c(500, 300, 100),
                    optim.type = 'adam', learn.rates = 0.01, val.prop = 0)

pred <- predict(NN, newdata = X_test)
Y_Pred <- as.numeric(pred$predictions)
res <- tibble(Y_pred, Y_test)
gg <- ggplot(res, aes(x=y_Pred, y=Y_test)) + geom_point()
