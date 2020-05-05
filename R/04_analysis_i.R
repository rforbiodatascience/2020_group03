<<<<<<< HEAD
# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library("ggseqlogo")
library("protr")
library("tensorflow")
library("keras")

# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
aug_clean_data_set_2 <- read_tsv(file = "./data/03_aug_data_set_2.tsv")

# Wrangle data
# ------------------------------------------------------------------------------
df <- aug_clean_data_set_2 %>%
  filter(score>10)


bl62 = read.table(file = 'data/_raw/BLOSUM62_.txt')

X_train <- df %>%
  pull(sequence) %>%
  encode_peptide(m = bl62)

y_train = df %>%
  pull(score)



# Visualise data
# ---------------------------
# ht <- df %>%
#     ggplot(aes(x = score)) +
#     geom_histogram()

# Write data
# ------------------------------------------------------------------------------
#write_tsv(...)
#ggsave(...)
=======
## bego
# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library("devtools")
library(keras)
# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_functions.R")

# Load data
# ------------------------------------------------------------------------------
aug_clean_data_set_1 <- read_tsv(file = "./data/03_aug_data_set_1.tsv")

# Wrangle data
# ------------------------------------------------------------------------------
#clean_data <- read.table(file = './data/clean_test.txt')
data <- read_table2(file = "./data/clean_test.txt", col_names = FALSE)
# Remove non-important columns
data <- data %>% 
  rename(
    peptide = X1,
    activity = X2,
  ) %>%  select(c(1,2))

data$len = str_length(data$peptide)
data <- data %>% group_by(len) %>% filter(len == 15) 

# Add new columns (encoding)
encoded_seq = data %>%
  pull(peptide) %>% 
  encode_peptide(m = "blosum62")

# Model data
# ------------------------------------------------------------------------------
#my_data_clean_aug %>% ...
# The ouput of encode_peptide function is a matrix, we need to convert it to a dataframe
df_encoded_seq <- data.frame(matrix(unlist(encoded_seq), nrow=8910, byrow=T),stringsAsFactors=FALSE)
# We add the activity column to make the models 
df_encoded_seq$activity <- data$activity
# 
peptide_to_predict = 'LNVTSEDLGKTFSVG'
peptide_to_predict_data = df_encoded_seq[1,]
# 
df_encoded_seq <- df_encoded_seq[-1,]
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
newpeptide = 'LNVTSEDLGKTLEAK'
newpep_enc = encode_peptide(newpeptide, "blosum62")
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

# ANN + other regression model
# https://www.rdocumentation.org/packages/neuralnet/versions/1.44.2/topics/neuralnet

index <- sample(1:nrow(df_encoded_seq),round(0.75*nrow(df_encoded_seq)))
train <- df_encoded_seq[index,]
test <- df_encoded_seq[-index,]
lm.fit <- glm(activity~., data=train)
summary(lm.fit)
pr.lm <- predict(lm.fit,test)
MSE.lm <- sum((pr.lm - test$activity)^2)/nrow(test)
predict(lm.fit, df_newpep_enc)

# to normalize (by column)
maxs <- apply(df_encoded_seq, 2, max) 
mins <- apply(df_encoded_seq, 2, min)
scaled <- as.data.frame(scale(df_encoded_seq, center = mins, scale = maxs - mins))
train_ <- scaled[index,]
test_ <- scaled[-index,]

library(neuralnet)
n <- names(train_)
f <- as.formula(paste("activity ~", paste(n[!n %in% "activity"], collapse = " + ")))
nn <- neuralnet(f,data=train_,hidden=c(5,3),linear.output=T)

# if no normalize ()
#library(neuralnet)
#n <- names(train)
#f <- as.formula(paste("activity ~", paste(n[!n %in% "activity"], collapse = " + ")))
#nn <- neuralnet(f,data=train,hidden=c(5,3),linear.output=T)

# The problem here is that we have too many variables that a Neural Net is computationally demanding
# therefore -> reduce the dimensionality with the variables that matter (see in lm) and do ANN


# Visualise data
# ---------------------------
ht <- aug_clean_data_set_1 %>%
    ggplot(aes(x = E3_score)) +
    geom_histogram()

# Write data
# ------------------------------------------------------------------------------
#write_tsv(...)
#ggsave(...)
>>>>>>> master
