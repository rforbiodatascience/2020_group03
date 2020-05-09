# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library("glmnet")
library("caret")
library("GauPro")
library("UniprotR")

# Define functions
# ------------------------------------------------------------------------------
source(file = "./R/99_functions.R")


# Define functions
# ------------------------------------------------------------------------------
source(file = "./R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
#clean_data_set_1 <- read_tsv(file = "./data/02_clean_data_set_1.tsv")
clean_data_set_2 <- read_tsv(file = "./data/02_clean_data_set_2.tsv")
#clean_data_set_3 <- read_tsv(file = "./data/02_clean_data_set_3.tsv")
#clean_data_set_4 <- read_tsv(file = "./data/02_clean_data_set_4.tsv")

# Wrangle data
# ------------------------------------------------------------------------------
#aug_data_set_1 <- clean_data_set_1 %>%
#  mutate(variant = case_when(variant=="NA-NA" ~ "D2D",
#                             !variant=="NA-NA" ~ variant))

# augmenting with sequence
#aug_data_set_1 <- convert_variant_to_sequence(aug_data_set_1, "P38398")
aug_data_set_2 <- convert_variant_to_sequence(clean_data_set_2, "P28482")
#aug_data_set_3 <- convert_variant_to_sequence(clean_data_set_3, "Q5SW96")
#aug_data_set_4 <- convert_variant_to_sequence(clean_data_set_4, "P04147")

data <- aug_data_set_2 %>%
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
  encode_peptide(m = "z_scales")

X_test <- test %>%
  pull(sequence_to_model) %>%
  encode_peptide(m = "z_scales")

Y_test <- test %>%
  select(score)

Y_train <- training %>%
  select(score)

Y_test <- unname(as.matrix(Y_test))
Y_train <- unname(as.matrix(Y_train))

X_train <- unname(as.matrix(X_train))

print("Lasso")
fit.lasso.cv <- cv.glmnet(X_train, Y_train, type.measure="mae", alpha=1, 
                          family="gaussian",
                          standardize=TRUE)

#print("Ridge")
#fit.ridge.cv <- cv.glmnet(X_train, Y_train, type.measure="mse", alpha=0,
#                          family="gaussian",
#                          standardize=TRUE)

# print("ElasticNet")
# fit.elnet.cv <- cv.glmnet(X_train, Y_train, type.measure="mae", alpha=.5,
#                           family="gaussian",
#                           standardize=TRUE)

X_test <- unname(as.matrix(X_test))

Y_pred_Lasso <- predict(fit.lasso.cv, newx = X_test, s = "lambda.min")
#Y_pred_Ridge <- predict(fit.ridge.cv, newx = X_test, s = "lambda.min")
# Y_pred_ElasticNet <- predict(fit.elnet.cv, newx = X_test, s = "lambda.min")

results = tibble(Y_test, Y_pred_Ridge)

results <- results %>%
  mutate(score = Y_test) %>%
  full_join(test, by = "score") %>%
  select(variant, score, Y_test, Y_pred_Lasso) %>%
  pivot_longer(c(-variant, -score), names_to = "model", values_to = "prediction")
  

gg <- ggplot(results, aes(x=score, y=prediction, color=model)) + geom_point()