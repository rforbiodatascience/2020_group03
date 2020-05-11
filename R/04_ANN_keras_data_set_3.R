# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library("devtools")
library("keras")
library("ggplot2")

# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
#data_set_1 <- read_tsv(file = "./data/03_aug_data_set_1.tsv")
#data_set_2 <- read_tsv(file = "./data/03_aug_data_set_2.tsv")
data_set_3 <- read_tsv(file = "./data/03_aug_data_set_3.tsv")
#data_set_4 <- read_tsv(file = "./data/03_aug_data_set_4.tsv")

# Wrangle data
# ------------------------------------------------------------------------------
data_set_3 <- data_set_3 %>%
  mutate(variant = case_when(variant=="NA-NA" ~ "D2D",
                             !variant=="NA-NA" ~ variant))

data_set_3 %>% glimpse()

data <- data_set_3 %>% 
  rename(
    peptide = sequence,
    activity = score,
  ) %>%  select(c(2,6))

max_value <- data_set_3 %>% 
  summarise(max = max(mutation_position)) %>% 
  as.integer()

data <- data %>% 
  mutate(peptide = strtrim(peptide, max_value))  %>%
  mutate(len = str_length(peptide))

# Add new columns (encoding)
# change
encoded_seq <- data %>%
  pull(peptide) %>% 
  encode_peptide(m = "blosum62")
nrow <- dim(encoded_seq)[1]
# The ouput of encode_peptide function is a matrix, we need to convert it to a dataframe
df_encoded_seq <- data.frame(matrix(unlist(encoded_seq), nrow=nrow, byrow=T),stringsAsFactors=FALSE)

# To add the activity
df_encoded_seq$activity <- data$activity

# for having the rownames too
sequences <- data.frame(row.names(encoded_seq))
df_encoded_rownames <- cbind(sequences, encoded_seq)
df_encoded_rownames <- df_encoded_seq[,-1]

# Model data
# ------------------------------------------------------------------------------

## Multiple linear model
model <- lm(activity ~., data = df_encoded_seq)
# check lm
par(mfrow= c(2,2))
plot(model)
summary(model)

# predict using encode function
newpeptide <- 'LNVTSEDLGKTLEAK'
newpep_enc <- encode_peptide(newpeptide, "blosum62")
df_newpep_enc <- data.frame(matrix(unlist(newpep_enc), nrow=1, byrow=T),stringsAsFactors=FALSE)
predict(model, df_newpep_enc)

## ANN
# Example of custom metric, here pearson's correlation coefficient, i.e.
# equivalent to cor(x, y, method = "pearson"), but note how we need to use the
# keras (tensorflow) methods 'k_*'
metric_pcc = custom_metric("pcc", function(y_true, y_pred) {
  mu_y_true = k_mean(y_true)
  mu_y_pred = k_mean(y_pred)
  r = k_sum( (y_true - mu_y_true) * (y_pred - mu_y_pred) ) /
    ( k_sqrt(k_sum( k_square(y_true - mu_y_true) )) *
        k_sqrt(k_sum( k_square(y_pred - mu_y_pred) )) )
  return(r)
})

# Set nn data and test/training partitions
test_f = 0.20
nn_dat = df_encoded_seq %>%
  mutate_if(is.ordered, as.numeric) %>% 
  mutate(partition = sample(x = c('test', 'train'),
                            size = nrow(.),
                            replace = TRUE,
                            prob = c(test_f, 1 - test_f)))

# as we are going to predict activity (y) we need to exclude it from X_train and X_test
# Set training data
X_train = nn_dat %>%
  filter(partition == "train") %>%
  select(-activity, -partition) %>%
  as.matrix
y_train = nn_dat %>%
  filter(partition == "train") %>%
  pull(activity)

# Set test data
X_test = nn_dat %>%
  filter(partition == "test") %>%
  select(-activity, -partition) %>%
  as.matrix
y_test = nn_dat %>%
  filter(partition == "test") %>%
  pull(activity)

# General lineal model given the partitions
model_2 <- glm.fit(X_train, y_train, family = gaussian())
summary.glm(model_2)
predict.glm(model_2)

# Linear model with partitions
train_df = nn_dat %>%
  filter(partition == "train") %>%
  select(-partition)
X_test_df = nn_dat %>%
  filter(partition == "test") %>%
  select(-activity, -partition)

model_3 <- lm(activity ~., data = train_df)
par(mfrow = c(2, 2))  # Split the plotting panel into a 2 x 2 grid
plot(model_3)
summary(model_3)
y_pred <- predict(model_3, X_test_df)

# plot
ggplot(train_df, aes(x = X1, y = activity, color = activity) ) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black", size = 0),
    axis.text.x = element_text(face = "bold", color = "#000000"),
    axis.text.y = element_text(face = "bold", color = "#000000")
  )

# Set hyperparameters
n_epochs      = 10
batch_size    = 200
loss          = 'mean_squared_error'
learning_rate = 0.1
optimzer      = optimizer_adam(lr = learning_rate)
h1_activation = 'relu'
h1_n_hidden   = 2
h2_activation = 'relu'
h2_n_hidden   = 2
h3_activation = 'relu'
h3_n_hidden   = 2
o_activation  = 'linear'

# Set architecture
model_ann = keras_model_sequential() %>% 
  layer_dense(units = h1_n_hidden,
              activation = h1_activation,
              input_shape = ncol(X_train)) %>%
  layer_dense(units = h2_n_hidden,
              activation = h2_activation) %>% 
  layer_dense(units = h3_n_hidden,
              activation = h3_activation) %>% 
  layer_dense(units = 1,
              activation = o_activation)

# Compile model
model_ann %>% compile(
  loss      = loss,
  optimizer = optimzer,
  metrics   = metric_pcc # Note the custom metric here
)


# View model
model_ann %>% summary %>% print

# Train model
# ------------------------------------------------------------------------------

# Fit model on training data
history = model_ann %>%
  fit(x = X_train,
      y = y_train,
      epochs = n_epochs,
      batch_size = batch_size,
      validation_split = 0
  )
# Evaluate model
# ------------------------------------------------------------------------------

# Calculate performance on test data
y_test_true = y_test
y_test_pred = model_ann %>% predict(X_test) %>% as.vector
pcc_test = round(cor(y_test_pred, y_test_true, method = "pearson"), 3)

# Calculate performance on training data
y_train_true = y_train
y_train_pred = model_ann %>% predict(X_train) %>% as.vector
pcc_train = round(cor(y_train_pred, y_train_true, method = "pearson"), 3)

# Compile data for plotting
d_perf = bind_rows(tibble(y_pred = y_test_pred,
                          y_true = y_test_true,
                          partition = 'test'),
                   tibble(y_pred = y_train_pred,
                          y_true = y_train_true,
                          partition = 'train'))

# Visualise performance
# ------------------------------------------------------------------------------
title = "Performance of ANN Regression model"
sub_title = paste0("Test PCC = ", pcc_test, ", training PCC = ", pcc_train, ".")
d_perf %>%
  ggplot(aes(x = y_pred, y = y_true)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  ggtitle(label = title, subtitle = sub_title) +
  xlab("Predicted Activity") +
  ylab("Actual Activity") +
  facet_wrap(~partition) +
  theme_bw()