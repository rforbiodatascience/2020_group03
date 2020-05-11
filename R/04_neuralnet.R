# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library("caret")
library("UniprotR")
library("ANN2")
library("yardstick")


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
  drop_na(score)# %>%
  #mutate(score_norm = (score - mean(score))/sd(score)) %>%
  #mutate(score_stad = (score - min(score)) / (max(score)-min(score)))



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
data_set_1 <- data_set_4 %>%
  mutate(variant = case_when(variant=="NA-NA" ~ "D2D",
                             !variant=="NA-NA" ~ variant))



# Select dataset and wrangle
# ------------------------------------------------------------------------------
test_set = data_set_1
filename = "./doc/ann/ann_corr_data_set_4.png"


data <- test_set %>%
  drop_na(score)


# ANN2 hyperparameters
# ------------------------------------------------------------------------------
epochs = 10 # number of epochs
hidden_layers = c(500,500) # number and size of hidden layers
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

print('NN running...')

NN <- neuralnetwork(X = X_train, y = Y_train,
                    hidden.layers = hidden_layers,
                    loss.type="squared",
                    optim.type = 'adam',
                    regression = TRUE,
                    activ.functions = "tanh",
                    learn.rates = 0.01, val.prop = 0,
                    n.epochs = epochs)

pred <- predict(NN, newdata = X_test)
Y_pred <- as.numeric(pred$predictions)
res <- tibble(Y_pred, Y_test)



corr_plot <- ggplot(res, aes(x=Y_pred, y=Y_test)) +
  geom_point() +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black", size = 0),
    axis.text.x = element_text(face = "bold", color = "#000000"),
    axis.text.y = element_text(face = "bold", color = "#000000")
  )

plot(gg)

ggsave(plot=corr_plot, filename)
