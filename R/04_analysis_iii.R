## bego
# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library("devtools")
library("keras")
# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_functions.R")

# Load data
# ------------------------------------------------------------------------------
aug_clean_data_set_1 <- read_tsv(file = "./data/03_aug_data_set_1.tsv")

# Wrangle data
# ------------------------------------------------------------------------------
#clean_data <- read.table(file = './data/clean_test.txt')
data <- read_table2(file = "./data/test_clean.txt", col_names = FALSE)

data <- data %>% 
  rename(
    peptide = X1,
    activity = X2,
  ) %>%  select(c(1,2))

data$len = str_length(data$peptide)
data <- data %>% group_by(len) %>% filter(len == 15) 
#write.csv(data,"data/data_clean_ii.csv", row.names = FALSE)

# Add new columns (encoding)
encoded_seq = data %>%
  pull(peptide) %>% 
  encode_peptide(m = "blosum62")
# The ouput of encode_peptide function is a matrix, we need to convert it to a dataframe
df_encoded_seq <- data.frame(matrix(unlist(encoded_seq), nrow=8910, byrow=T),stringsAsFactors=FALSE)
#write.csv(df_encoded_seq, './data/df_encoded_seq.csv')
#write.csv(data, './data/data.csv')
pep_unique <- unique(data$peptide)
df_encoded_seq$activity <- data$activity


# for having the rownames too
sequences <- data.frame(row.names(encoded_seq))
df_encoded_rownames <- cbind(sequences, encoded_seq)
df_encoded_rownames <- df_encoded_seq[,-1]

#write.csv(df_encoded_seq, './data/df_encoded_seq_activity.csv')
#write.csv(df_encoded_seq, './data/df_encoded_seq_rownames.csv')

# the dataframe needs to be reduce in its rows, some peptides are repeated
#df <- read.csv('./data/df_encoded_seq_rownames.csv')
#df_laura <- read.csv('./data/encoded_seq.csv')
# Model data
# ------------------------------------------------------------------------------

# Multiple linear model
model <- lm(activity ~., data = df_encoded_seq)
# check lm
par(mfrow= c(2,2))
plot(model)
summary(model)
#par(mfrow = c(1,2))
#qqnorm(model$residuals)
#qqline(model$residuals)
#qqPlot(model)	

# predict using encode function
newpeptide <- 'LNVTSEDLGKTLEAK'
newpep_enc <- encode_peptide(newpeptide, "blosum62")
df_newpep_enc <- data.frame(matrix(unlist(newpep_enc), nrow=1, byrow=T),stringsAsFactors=FALSE)
predict(model, df_newpep_enc)

# Given we have a lot of variables I will try with the 'step' function to reduce dimensionality ??
# computationally demanding
# model_step <- step(model)

# Reduce the dimensionality (just a try)
df_coef <- summary(model)$coefficients
relevant_variables <- df_coef[df_coef[,4] < 0.05,]
model_rel <- lm(activity ~ X25 + X27 + X41+ X48+ X55+ X114+ X127+ X154+ X164+ X214+ X235+ X266+ X285+ X286+ X297+ X321+ X329+ X331+ X332, data = df_encoded_seq)
predict(model_rel, df_newpep_enc)

relevant_variables <- df_coef[df_coef[,4] < 0.01,]
model_rel_2 <- lm(activity ~ X114+ X214+ X266, data = df_encoded_seq)
predict(model_rel_2, df_newpep_enc)

AIC(model, model_rel, model_rel_2) # model_rel would be selected for less AIC
BIC(model, model_rel, model_rel_2)

# predict th first row was eliminated
peptide_to_predict = 'LNVTSEDLGKTFSVG'
newpep_enc = encode_peptide(peptide_to_predict, "blosum62")
df_newpep_enc <- data.frame(matrix(unlist(newpep_enc), nrow=1, byrow=T),stringsAsFactors=FALSE)
predict(model, df_newpep_enc)
predict(model_rel, df_newpep_enc)
predict(model_rel_2, df_newpep_enc)

# https://github.com/leonjessen/RPharma2019/blob/master/R/04_diamonds_regression.R
# ANN + other regression model
# https://www.rdocumentation.org/packages/neuralnet/versions/1.44.2/topics/neuralnet

#index <- sample(1:nrow(df_encoded_seq),round(0.75*nrow(df_encoded_seq)))
#train <- df_encoded_seq[index,]
#test <- df_encoded_seq[-index,]
#lm.fit <- glm(activity~., data=train)
#summary(lm.fit)
#pr.lm <- predict(lm.fit,test)
#MSE.lm <- sum((pr.lm - test$activity)^2)/nrow(test)
#predict(lm.fit, df_newpep_enc)

# to normalize (by column)
#maxs <- apply(df_encoded_seq, 2, max) 
#mins <- apply(df_encoded_seq, 2, min)
#scaled <- as.data.frame(scale(df_encoded_seq, center = mins, scale = maxs - mins))
#train_ <- scaled[index,]
#test_ <- scaled[-index,]

#library(neuralnet)
#n <- names(train_)
#f <- as.formula(paste("activity ~", paste(n[!n %in% "activity"], collapse = " + ")))
#nn <- neuralnet(f,data=train_,hidden=c(5,3),linear.output=T)

# if no normalize ()
#library(neuralnet)
#n <- names(train)
#f <- as.formula(paste("activity ~", paste(n[!n %in% "activity"], collapse = " + ")))
#nn <- neuralnet(f,data=train,hidden=c(5,3),linear.output=T)

# The problem here is that we have too many variables that a Neural Net is computationally demanding
# therefore -> reduce the dimensionality with the variables that matter (see in lm) and do ANN

# Keras

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
