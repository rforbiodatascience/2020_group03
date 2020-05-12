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
RMSE_data_set_1_ANN_z <- read_tsv(file = "./data/ANN_zs_RMSE_data_set_1.tsv")
RMSE_data_set_2_ANN_z <- read_tsv(file = "./data/ANN_zs_RMSE_data_set_2.tsv")
RMSE_data_set_3_ANN_z <- read_tsv(file = "./data/ANN_zs_RMSE_data_set_3.tsv")
RMSE_data_set_4_ANN_z <- read_tsv(file = "./data/ANN_zs_RMSE_data_set_4.tsv")
RMSE_data_set_1_glmnet_z <- read_tsv(file = "./data/glmnet_zs_RMSE_data_set_1.tsv")
RMSE_data_set_2_glmnet_z <- read_tsv(file = "./data/glmnet_zs_RMSE_data_set_2.tsv")
RMSE_data_set_3_glmnet_z <- read_tsv(file = "./data/glmnet_zs_RMSE_data_set_3.tsv")
RMSE_data_set_4_glmnet_z <- read_tsv(file = "./data/glmnet_zs_RMSE_data_set_4.tsv")

RMSE_data_set_1_ANN_b <- read_tsv(file = "./data/ANN_blosum62_RMSE_data_set_1.tsv")
RMSE_data_set_2_ANN_b <- read_tsv(file = "./data/ANN_blosum62_RMSE_data_set_2.tsv")
RMSE_data_set_3_ANN_b <- read_tsv(file = "./data/ANN_blosum62_RMSE_data_set_3.tsv")
RMSE_data_set_4_ANN_b <- read_tsv(file = "./data/ANN_blosum62_RMSE_data_set_4.tsv")
RMSE_data_set_1_glmnet_b <- read_tsv(file = "./data/glmnet_blosum62_RMSE_data_set_1.tsv")
RMSE_data_set_2_glmnet_b <- read_tsv(file = "./data/glmnet_blosum62_RMSE_data_set_2.tsv")
RMSE_data_set_3_glmnet_b <- read_tsv(file = "./data/glmnet_blosum62_RMSE_data_set_3.tsv")
RMSE_data_set_4_glmnet_b <- read_tsv(file = "./data/glmnet_blosum62_RMSE_data_set_4.tsv")



RMSE_data_set_1_ANN_z <- RMSE_data_set_1_ANN_z %>%
  add_column(dataset = "data_set_1",
             model = "ANN",
             scale = "z_scales")

RMSE_data_set_2_ANN_z <- RMSE_data_set_2_ANN_z %>%
  add_column(dataset = "data_set_2",
             model = "ANN",
             scale = "z_scales")

RMSE_data_set_3_ANN_z <- RMSE_data_set_3_ANN_z %>%
  add_column(dataset = "data_set_3",
             model = "ANN",
             scale = "z_scales")

RMSE_data_set_4_ANN_z <- RMSE_data_set_4_ANN_z %>%
  add_column(dataset = "data_set_4",
             model = "ANN",
             scale = "z_scales")

RMSE_data_set_1_glmnet_z <- RMSE_data_set_1_glmnet_z %>%
  add_column(dataset = "data_set_1",
             model = "ElasticNet",
             scale = "z_scales")

RMSE_data_set_2_glmnet_z <- RMSE_data_set_2_glmnet_z %>%
  add_column(dataset = "data_set_2",
             model = "ElasticNet",
             scale = "z_scales")

RMSE_data_set_3_glmnet_z <- RMSE_data_set_3_glmnet_z %>%
  add_column(dataset = "data_set_3",
             model = "ElasticNet",
             scale = "z_scales")

RMSE_data_set_4_glmnet_z <- RMSE_data_set_4_glmnet_z %>%
  add_column(dataset = "data_set_4",
             model = "ElasticNet",
             scale = "z_scales")




RMSE_data_set_1_ANN_b <- RMSE_data_set_1_ANN_b %>%
  add_column(dataset = "data_set_1",
             model = "ANN",
             scale = "blosum62")

RMSE_data_set_2_ANN_b <- RMSE_data_set_2_ANN_b %>%
  add_column(dataset = "data_set_2",
             model = "ANN",
             scale = "blosum62")

RMSE_data_set_3_ANN_b <- RMSE_data_set_3_ANN_b %>%
  add_column(dataset = "data_set_3",
             model = "ANN",
             scale = "blosum62")

RMSE_data_set_4_ANN_b <- RMSE_data_set_4_ANN_b %>%
  add_column(dataset = "data_set_4",
             model = "ANN",
             scale = "blosum62")

RMSE_data_set_1_glmnet_b <- RMSE_data_set_1_glmnet_b %>%
  add_column(dataset = "data_set_1",
             model = "ElasticNet",
             scale = "blosum62")

RMSE_data_set_2_glmnet_b <- RMSE_data_set_2_glmnet_b %>%
  add_column(dataset = "data_set_2",
             model = "ElasticNet",
             scale = "blosum62")

RMSE_data_set_3_glmnet_b <- RMSE_data_set_3_glmnet_b %>%
  add_column(dataset = "data_set_3",
             model = "ElasticNet",
             scale = "blosum62")

RMSE_data_set_4_glmnet_b <- RMSE_data_set_4_glmnet_b %>%
  add_column(dataset = "data_set_4",
             model = "ElasticNet",
             scale = "blosum62")



RMSE_summary <- bind_rows(RMSE_data_set_1_ANN_z,
                          RMSE_data_set_2_ANN_z,
                          RMSE_data_set_3_ANN_z,
                          RMSE_data_set_4_ANN_z,
                          RMSE_data_set_1_glmnet_z,
                          RMSE_data_set_2_glmnet_z,
                          RMSE_data_set_3_glmnet_z,
                          RMSE_data_set_4_glmnet_z,
                          RMSE_data_set_1_ANN_b,
                          RMSE_data_set_2_ANN_b,
                          RMSE_data_set_3_ANN_b,
                          RMSE_data_set_4_ANN_b,
                          RMSE_data_set_1_glmnet_b,
                          RMSE_data_set_2_glmnet_b,
                          RMSE_data_set_3_glmnet_b,
                          RMSE_data_set_4_glmnet_b)

RMSE_plot <- ggplot(RMSE_summary, aes(x = dataset, y = .estimate, fill = model)) +
  geom_bar(stat="identity", position="dodge") +
  theme_classic() +
  ylab("RMSE") +
  facet_grid(cols = vars(scale)) +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black", size = 0),
    axis.text.x = element_text(face = "bold", color = "#000000", angle = 90),
    axis.text.x.bottom = element_text(vjust = 0.5, angle = -90, face = "bold", color = "#000000"))


ggsave(plot = RMSE_plot, "./results/06_model_summaries/RMSE_blosum62_z_scale.png")