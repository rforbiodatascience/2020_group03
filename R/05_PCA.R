# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())



# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library("stats")
library("broom")
library("dplyr")
library("ggpubr")



# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_project_functions.R")



# Load data
# ------------------------------------------------------------------------------
df_encoded <- read.csv("encoded_seq.csv", header = TRUE)
seq_annotation <- read.csv("data_clean_ii.csv", header = TRUE)



# Wrangle data
# ------------------------------------------------------------------------------
# Make sure the length of annotations is the same as encoded output

genes_seq_annotation <- unique(seq_annotation$peptide)
unique_rows_annotation <- unique(seq_annotation)

genes_encoded <- unique(df_encoded$X)
df_encoded<- unique(df_encoded)


#Creation of two group labels based on activity measurements to use as a colour option later in plots
seq_annotation <- seq_annotation %>%
    mutate(Range_ = if_else(Mean_activity >= 0.5, 'High activity', 'Low activity'))



# Visualise data before PCA
# ------------------------------------------------------------------------------

#Activity distribution before and after the mean computation for common strings

dist1 <- ggplot(seq_annotation, aes(x=activity)) +
    geom_density() +
    geom_vline(data = seq_annotation, aes(xintercept=mean(activity)),
               linetype="dashed", size=1)+
  theme_classic() + labs(x = "Activity", y = "Density") + ggtitle("Activity distribution")+ 
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black", size = 0),
    axis.text.x = element_text(face = "bold", color = "#000000"),
    axis.text.y = element_text(face = "bold", color = "#000000")
  )

seq_annotation <- seq_annotation %>% 
   dplyr::group_by(seq_annotation$peptide) %>%
   dplyr::summarize(Mean_activity = mean(activity))



dist2 <- ggplot(seq_annotation, aes(x=Mean_activity)) +
    geom_density() +
    geom_vline(aes(xintercept=mean(seq_annotation$Mean_activity)),
               linetype="dashed", size=1) +
  theme_classic() + labs(x = "Mean activity", y = "Density") + ggtitle("Mean activity distribution")+
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black", size = 0),
    axis.text.x = element_text(face = "bold", color = "#000000"),
    axis.text.y = element_text(face = "bold", color = "#000000")
  )

png(filename="./data/04_activity_distribution_plot")
activity_distribution_plot <- ggarrange(dist1, dist2, ncol = 2, nrow = 1)
dev.off()




# PCA
# ------------------------------------------------------------------------------

# Compute PCA on encoded sequences -- zero centered values and unit variance

my_pca <- df_encoded %>% select(-X) %>% 
  
  
  prcomp(center = TRUE, scale. = TRUE)


my_pca %>% tidy("pcs")



# Scree plot
png(filename="./data/04_scree_plot_pca")
scree_plot_pca <- my_pca %>% tidy("pcs") %>% 
  ggplot(aes(x = PC, y = percent)) +
  geom_col() +  ggtitle("Variance explained for the principal components")+ 
  theme_classic() + labs(x = "Principal components", y = "Variance percentage") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black", size = 0),
    axis.text.x = element_text(face = "bold", color = "#000000"),
    axis.text.y = element_text(face = "bold", color = "#000000")
  )
dev.off()



# PC projections using the two groups to color

my_pca %>% tidy("samples")

my_pca_aug <- my_pca %>% augment(df_encoded)


png(filename="./data/04_PCA_PROJECTION_PLOT")
pca_projection_plot <- my_pca_aug %>% ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = seq_annotation$Range_)) +
  geom_point() + 
  theme_classic() + labs(x = "PC1", y = "PC2") +  ggtitle("Projection into the first two principal components")+ 
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank(),
    axis.line = element_line(colour = "black", size = 0),
    axis.text.x = element_text(face = "bold", color = "#000000"),
    axis.text.y = element_text(face = "bold", color = "#000000")
  )
dev.off()











