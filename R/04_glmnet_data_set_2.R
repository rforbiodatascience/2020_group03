# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library("glmnet")
library("caret")
library("UniprotR")
library("yardstick")


# Define functions
# ------------------------------------------------------------------------------
source(file = "./R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
data_set_1 <- read_tsv(file = "./data/03_aug_data_set_1.tsv")
data_set_2 <- read_tsv(file = "./data/03_aug_data_set_2.tsv")
data_set_3 <- read_tsv(file = "./data/03_aug_data_set_3.tsv")
data_set_4 <- read_tsv(file = "./data/03_aug_data_set_4.tsv")

# Wrangle data
# ------------------------------------------------------------------------------
data_set_1 <- data_set_1 %>%
  mutate(variant = case_when(variant=="NA-NA" ~ "D2D",
                             !variant=="NA-NA" ~ variant))



# Select dataset and wrangle
# ------------------------------------------------------------------------------
test_set = data_set_2
filename = "./doc/glmnet/corr_data_set_2.png"
RMSE_path = "./data/glmnet_RMSE_data_set_2.tsv"


data <- test_set %>%
  drop_na(score)


# Elastic net hyperparameters
# ------------------------------------------------------------------------------
alpha_ = 0.2 # mix between ridge and lasso
s_ = 0.01 # regularization strenght
scale = "z_scales" # amino acid descriptor scale
train_size = 0.75 # train sample size
seed_value = 42 # random seed


# Get sequence window and truncate
# ------------------------------------------------------------------------------
sequence_window <- data %>%
  mutate(mutation_position = as.integer(mutation_position)) %>%
  pull(mutation_position)

data <- data %>%
  mutate(sequence_to_model = str_sub(sequence, min(sequence_window), max(sequence_window)))

# Make score response variable
# ------------------------------------------------------------------------------
Y <- data  %>%
  select(score)
 
Y <- unname(as.matrix(Y))

# Partion data
# ------------------------------------------------------------------------------
set.seed(seed_value)
partition <- createDataPartition(y = Y, p = train_size , list = F)
training = data[partition, ]
test <- data[-partition, ]

# Encode sequences and make test and train sets
# ------------------------------------------------------------------------------
X_train <- training %>%
  pull(sequence_to_model) %>%
  encode_peptide(m = scale)

X_test <- test %>%
  pull(sequence_to_model) %>%
  encode_peptide(m = scale)

Y_test <- test %>%
  select(score)

Y_train <- training %>%
  select(score)

Y_test <- unname(as.matrix(Y_test))
Y_train <- unname(as.matrix(Y_train))
X_train <- unname(as.matrix(X_train))


# Fit elasticnet model
# ------------------------------------------------------------------------------
fit.elasticnet <- glmnet(X_train, Y_train, type.measure="mse", alpha = alpha_, s = s_,
                         family="gaussian",
                         standardize=TRUE)

# Preduct using elasticnet model
# ------------------------------------------------------------------------------
Y_pred_ElasticNet <- predict(fit.elasticnet, newx = X_test, s = s_, alpha=alpha_)
results = tibble(Y_test, Y_pred_ElasticNet)
RMSE_table = results


# Pivot results and test response variables
# ------------------------------------------------------------------------------
results <- results %>%
  mutate(score = Y_test) %>%
  full_join(test, by = "score") %>%
  select(variant, score, Y_test, Y_pred_ElasticNet) %>%
  pivot_longer(c(-variant, -score), names_to = "model", values_to = "prediction")
  
# Plot correlation between score and prediction
# ------------------------------------------------------------------------------
corr_plot <- ggplot(results, aes(x=score, y=prediction, color=model)) +
  geom_point() +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black", size = 0),
    axis.text.x = element_text(face = "bold", color = "#000000"),
    axis.text.y = element_text(face = "bold", color = "#000000")
  )

ggsave(plot=corr_plot, filename)


RMSE_table <- RMSE_table %>%
  select(Y_test, Y_pred_ElasticNet)

RMSE_ = rmse(RMSE_table, Y_test[,1], Y_pred_ElasticNet[,"1"])
write_tsv(x = RMSE_, path = RMSE_path)