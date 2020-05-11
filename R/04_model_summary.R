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
RMSE_data_set_1_ANN <- read_tsv(file = "./data/ANN_RMSE_data_set_1.tsv")
RMSE_data_set_2_ANN <- read_tsv(file = "./data/ANN_RMSE_data_set_2.tsv")
RMSE_data_set_3_ANN <- read_tsv(file = "./data/ANN_RMSE_data_set_3.tsv")
RMSE_data_set_4_ANN <- read_tsv(file = "./data/ANN_RMSE_data_set_4.tsv")
RMSE_data_set_1_glmnet <- read_tsv(file = "./data/glmnet_RMSE_data_set_1.tsv")
RMSE_data_set_2_glmnet <- read_tsv(file = "./data/glmnet_RMSE_data_set_2.tsv")
RMSE_data_set_3_glmnet <- read_tsv(file = "./data/glmnet_RMSE_data_set_3.tsv")
RMSE_data_set_4_glmnet <- read_tsv(file = "./data/glmnet_RMSE_data_set_4.tsv")

RMSE_data_set_1_ANN <- RMSE_data_set_1_ANN %>%
  add_column(dataset = "data_set_1",
             model = "ANN")

RMSE_data_set_2_ANN <- RMSE_data_set_2_ANN %>%
  add_column(dataset = "data_set_2",
             model = "ANN")

RMSE_data_set_3_ANN <- RMSE_data_set_3_ANN %>%
  add_column(dataset = "data_set_3",
             model = "ANN")

RMSE_data_set_4_ANN <- RMSE_data_set_4_ANN %>%
  add_column(dataset = "data_set_4",
             model = "ANN")

RMSE_data_set_1_glmnet <- RMSE_data_set_1_glmnet %>%
  add_column(dataset = "data_set_1",
             model = "ElasticNet")

RMSE_data_set_2_glmnet <- RMSE_data_set_2_glmnet %>%
  add_column(dataset = "data_set_2",
             model = "ElasticNet")

RMSE_data_set_3_glmnet <- RMSE_data_set_3_glmnet %>%
  add_column(dataset = "data_set_3",
             model = "ElasticNet")

RMSE_data_set_4_glmnet <- RMSE_data_set_4_glmnet %>%
  add_column(dataset = "data_set_4",
             model = "ElasticNet")


RMSE_summary <- bind_rows(RMSE_data_set_1_ANN,
                          RMSE_data_set_2_ANN,
                          RMSE_data_set_3_ANN,
                          RMSE_data_set_4_ANN,
                          RMSE_data_set_1_glmnet,
                          RMSE_data_set_2_glmnet,
                          RMSE_data_set_3_glmnet,
                          RMSE_data_set_4_glmnet)


RMSE_plot <- ggplot(RMSE_summary, aes(x = dataset, y = .estimate, fill=model)) +
  geom_bar(stat="identity", position="dodge") +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black", size = 0),
    axis.text.x = element_text(face = "bold", color = "#000000"),
    axis.text.y = element_text(face = "bold", color = "#000000")
  )

ggsave(plot = RMSE_plot, "./doc/model_summaries/RMSE_blosum62.png")