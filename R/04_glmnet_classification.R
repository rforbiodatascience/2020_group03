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
library("yardstick")

# Define functions
# ------------------------------------------------------------------------------
source(file = "./R/99_functions.R")


# Define functions
# ------------------------------------------------------------------------------
source(file = "./R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
#clean_data_set_1 <- read_tsv(file = "./data/02_clean_data_set_1.tsv")
#clean_data_set_2 <- read_tsv(file = "./data/02_clean_data_set_2.tsv")
#clean_data_set_3 <- read_tsv(file = "./data/02_clean_data_set_3.tsv")
clean_data_set_4 <- read_tsv(file = "./data/02_clean_data_set_4.tsv")

# Wrangle data
# ------------------------------------------------------------------------------
#aug_data_set_1 <- clean_data_set_1 %>%
#  mutate(variant = case_when(variant=="NA-NA" ~ "D2D",
#                             !variant=="NA-NA" ~ variant))

# augmenting with sequence
#aug_data_set_1 <- convert_variant_to_sequence(aug_data_set_1, "P38398")
#aug_data_set_2 <- convert_variant_to_sequence(clean_data_set_2, "P28482")
#aug_data_set_3 <- convert_variant_to_sequence(clean_data_set_3, "Q5SW96")
aug_data_set_4 <- convert_variant_to_sequence(clean_data_set_4, "P04147")

data <- aug_data_set_4 %>%
  drop_na(score) %>%
  mutate(class = case_when(score<=0 ~ 1,
                           score>0 & score<=0.5 ~ 2,
                           score>0.5 & score<=1 ~ 3,
                           score>1 & score<= 1.5 ~ 4,
                           score>1.5 & score<=2 ~ 5,
                           score>2 ~ 6))
#
# wrap in function
#

sequence_window <- data %>%
  mutate(mutation_position = as.integer(mutation_position)) %>%
  pull(mutation_position)

data <- data %>%
  mutate(sequence_to_model = str_sub(sequence, min(sequence_window), max(sequence_window)))

Y <- data  %>%
  select(class)
 
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
  select(class)

Y_train <- training %>%
  select(class)

Y_test <- unname(as.matrix(Y_test))
Y_train <- unname(as.matrix(Y_train))

X_train <- unname(as.matrix(X_train))

fit <- cv.glmnet(X_train, Y_train, family = "multinomial")

Y_pred<-predict(fit, newx = X_test, s = "lambda.min", type = "class")

table<-confusionMatrix(Y_test, Y_pred)


plotTable <- table %>%
  mutate(goodbad = ifelse(table$Prediction == table$Reference, "good", "bad")) %>%
  group_by(Reference) %>%
  mutate(prop = Freq/sum(Freq))

ggplot(data = plotTable, mapping = aes(x = Prediction, y = Reference, fill = freq)) +
  geom_tile() +
  geom_text(aes(label = Freq), vjust = .5, fontface  = "bold", alpha = 1) +
  scale_fill_manual(values = c(good = "green", bad = "red")) +
  theme_bw() +
  ylim(rev(levels(table$Reference)))