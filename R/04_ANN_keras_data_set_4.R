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
#data_set_3 <- read_tsv(file = "./data/03_aug_data_set_3.tsv")
data_set_4 <- read_tsv(file = "./data/03_aug_data_set_4.tsv")

# Wrangle data
# ------------------------------------------------------------------------------
data_set_4 <- data_set_4 %>%
  mutate(variant = case_when(variant=="NA-NA" ~ "D2D",
                             !variant=="NA-NA" ~ variant))

data_set_4 %>% glimpse()

# we will make the model with the mutations 

max_value <- data_set_4 %>% 
  summarise(max = max(mutation_position)) %>% 
  as.integer()

if(max_value > 60){
  data_rwos <- data_set_4 %>%
    #filter(mutation_position < (max_value/3)*2 & mutation_position > (max_value/3)) %>%
    mutate(peptide = substr(sequence, (max_value/4), (max_value/4)*2))  %>%
    mutate(len = str_length(peptide))
  
  data <- data_rwos %>% 
    rename(
      activity = score,
    ) %>%  select(c(2,7))
} else {
  data <- data_set_4 %>% 
    rename(
      peptide = sequence,
      activity = score,
    ) %>%  select(c(2,6))
}

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


# Model data
# ------------------------------------------------------------------------------

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
test_f = 0.30
nn_dat = df_encoded_seq %>%
  filter_all(all_vars(. != 0)) %>% 
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

# Linear model with partitions
train_df = nn_dat %>%
  filter(partition == "train") %>%
  select(-partition)
X_test_df = nn_dat %>%
  filter(partition == "test") %>%
  select(-activity, -partition)

model <- lm(activity ~., data = train_df)
par(mfrow = c(2, 2))  # Split the plotting panel into a 2 x 2 grid
plot(model)
summary(model)
y_pred <- predict(model, X_test_df)

# plot
plot_lm <- ggplot(train_df, aes(x = X1, y = activity, color = activity) ) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black", size = 0),
    axis.text.x = element_text(face = "bold", color = "#000000"),
    axis.text.y = element_text(face = "bold", color = "#000000")
  )

# Save plot
ggsave(plot = plot_lm, "./results/06_ANN_plots/prelm_data_set_4.png")

# Save model
save(model, file="./results/06_ANN_data/prelm_data_set_4.Rdata")

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
# Accuracy 
model <- model_ann
acc_test <-  model %>% evaluate(X_test, y_test, batch_size=32)
acc_train <-  model %>% evaluate(X_train, y_train, batch_size=32)

# Calculate performance on test data with pearson 
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
ann_keras_data_set_4 <- d_perf %>%
  ggplot(aes(x = y_pred, y = y_true)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  ggtitle(label = title, subtitle = sub_title) +
  xlab("Predicted Activity") +
  ylab("Actual Activity") +
  facet_wrap(~partition) +
  theme_bw()


ggsave(plot = ann_keras_data_set_4, "./results/06_ANN_plot/ann_keras_data_set_4.png")

# Save model
# ------------------------------------------------------------------------------

save_model_hdf5(object = model,
                filepath = "./results/06_ANN_data/ann_keras_ds4.h5")

# Load the model (Note the use of the custom_objects argument)
loaded_model = load_model_hdf5(filepath = './results/06_ANN_data/ann_keras_ds4.h5',
                               custom_objects = list('pcc' = metric_pcc))
