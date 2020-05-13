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
RMSE_data_set_1_ANN_z <- read_tsv(file = "./results/06_ann/06_ann_RMSE_z_scales_data_set_1.tsv")
RMSE_data_set_2_ANN_z <- read_tsv(file = "./results/06_ann/06_ann_RMSE_z_scales_data_set_2.tsv")
RMSE_data_set_3_ANN_z <- read_tsv(file = "./results/06_ann/06_ann_RMSE_z_scales_data_set_3.tsv")
RMSE_data_set_4_ANN_z <- read_tsv(file = "./results/06_ann/06_ann_RMSE_z_scales_data_set_4.tsv")
RMSE_data_set_1_glmnet_z <- read_tsv(file = "./results/06_glmnet/06_glmnet_RMSE_z_scales_data_set_1.tsv")
RMSE_data_set_2_glmnet_z <- read_tsv(file = "./results/06_glmnet/06_glmnet_RMSE_z_scales_data_set_2.tsv")
RMSE_data_set_3_glmnet_z <- read_tsv(file = "./results/06_glmnet/06_glmnet_RMSE_z_scales_data_set_3.tsv")
RMSE_data_set_4_glmnet_z <- read_tsv(file = "./results/06_glmnet/06_glmnet_RMSE_z_scales_data_set_4.tsv")
 
RMSE_data_set_1_ANN_b <- read_tsv(file = "./results/06_ann/06_ann_RMSE_blosum62_data_set_1.tsv")
RMSE_data_set_2_ANN_b <- read_tsv(file = "./results/06_ann/06_ann_RMSE_blosum62_data_set_2.tsv")
RMSE_data_set_3_ANN_b <- read_tsv(file = "./results/06_ann/06_ann_RMSE_blosum62_data_set_3.tsv")
RMSE_data_set_4_ANN_b <- read_tsv(file = "./results/06_ann/06_ann_RMSE_blosum62_data_set_4.tsv")
RMSE_data_set_1_glmnet_b <- read_tsv(file = "./results/06_glmnet/06_glmnet_RMSE_blosum62_data_set_1.tsv")
RMSE_data_set_2_glmnet_b <- read_tsv(file = "./results/06_glmnet/06_glmnet_RMSE_blosum62_data_set_2.tsv")
RMSE_data_set_3_glmnet_b <- read_tsv(file = "./results/06_glmnet/06_glmnet_RMSE_blosum62_data_set_3.tsv")
RMSE_data_set_4_glmnet_b <- read_tsv(file = "./results/06_glmnet/06_glmnet_RMSE_blosum62_data_set_4.tsv")


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