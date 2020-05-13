# Define project functions
# ------------------------------------------------------------------------------

#' @description
#' `encode_peptide` returns a dataframe encoding the sequence in the method `m' specified.
#' All the peptides have to be the same length.
#'
#' @param x peptide, string
#' @param m method, string, options: blosum45, blosum62, blosum80, blosum90, pam30, pam70 and pam250
#' 
#' @return if all the parameters are logical will return a dataframe with the row.names = peptide and each column the encoding. The number of columns is len(peptide) x 20 amino acids. 
#' @example
#' encode_peptide('ATRALPLTW', "blosum62")
#' 
encode_peptide = function(x, m){
  mat <- paste("./data/_raw/",m, ".txt", sep="")
  encoder <- read.table(file = mat, row.names = 1, header =TRUE)
  X_enc = x %>%
    str_split('') %>%
    lapply(function(x_i){
      encoder[x_i,] %>%
        as.matrix %>%
        t() %>%
        matrix(nrow = 1, byrow = TRUE) %>%
        return
    })
  X_enc = do.call(rbind, X_enc)
  rownames(X_enc) = x
  return(X_enc)
}


#' make_sequence
#' @description
#' `make_sequence` convert variant information and WT sequence to variant sequence
#' Supports only single mutations
#'
#' @param Mutated_residue, string residue to mutate
#' @param Mutation_position, numeric, residue number to mutate
#' @param Mutation, string mutation to which amino acid
#' @param uniprot_ID, string uniprot_ID
#' 
#' @return sequence, string protein sequence
#' 
#' @examples
#' make_sequence("D", 2, "K", "37636")
make_sequence = function(Mutated_residue, Mutation_position, Mutation, uniprot_ID)
{
  WT_sequence <- as_tibble(GetSequences(uniprot_ID)) %>%
    pull(Sequence)
  sequence <- WT_sequence
  str_sub(sequence, Mutation_position, Mutation_position) <- Mutation
  return(sequence)
}


#' @description
#' `get_mutated_residue` extracts inforationm about WT residue
#' Supports only single mutations
#'
#' @param variant, string format example "D23E"
#'
#' @return mutated_residue, string
#'
#' @examples
#' get_mutated_residue("D2K")
get_mutated_residue <- function(variant) {
  mutated_residue = str_extract(variant, "[A-Z]")
  return(mutated_residue)
}


#' @description
#' `get_mutation_position` extracts information about mutation position
#'
#' @param variant, string
#'
#' @return mutation_position, numeric
#'
#' @examples
#' get_mutation_position("D2K")
get_mutation_position <- function(variant) {
  mutation_position = str_extract(variant, "[0-9]+")
  return(mutation_position)
}


#' @description
#' `get_mutation` extracts information about mutation to which amino_acid
#'
#' @param variant, string
#'
#' @return mutation, string
#' @export
#'
#' @examples
#' get_mutation("D2K")
get_mutation <- function(variant) {
  mutation = str_extract(variant, "[A-Z,*]$")
  return(mutation)
}


#' @description
#' `convert_variant_to_sequence` converts a tibble of variant_IDs to a tibble with corrensponding
#' sequences based on a uniprot ID
#' 
#' @param df, tibble of variant_IDs column name = variant
#' @param protein, string uniprot ID
#'
#' @return tibble with sequences
#'
#' @examples convert_variant_to_sequence(df, "74742")
convert_variant_to_sequence <- function(df, protein)  {
  clean_df <- df %>%
  mutate(
    mutated_residue = get_mutated_residue(variant), 
    mutation_position = get_mutation_position(variant),
    mutation = get_mutation(variant),
    sequence = make_sequence(mutated_residue, mutation_position, mutation, protein))
  return(clean_df)
}

#' @description
#' `aminoacid_type` returns an aminoacid type
#' 
#' @param residue, a residue in 1 letter code
#'
#' @return tibble with sequences
#'
#' @examples aminoacid_type("A")
#' 
aminoacid_type <- function(residue) {
  type <- case_when(
    residue == "R" | residue == "H" | residue == "K" ~ "basic",
    residue == "D" | residue == "E" ~ "acid",
    residue == "S" | residue == "T" ~ "hydroxilic",
    residue == "Q" | residue == "N" ~ "amidic",
    residue == "A" | residue == "G" | residue == "I" | residue == "L" | residue == "P" | residue == "V" ~ "aliphatic",
    residue == "F" | residue == "W" | residue == "Y" ~ "aromatic",
    residue == "C" | residue == "M" ~ "sulfur-containing"
  )
return(type)
}

#' @description
#' `zscale2col` returns 5 columns one for each of the z-scales values stored in a list
#' 
#' @param zscales, tibble of zscales in each row, a list of 5 zscales
#' @param df, df where to cbind the zscales
#'
#' @return tibble with 5 rows
#'
zscale2col <- function(df, zscales) {
  zscales <- data.frame(zscales)
  zscales <- as_tibble(cbind(nms = names(zscales), t(zscales))) %>%
    select(.data$Z1, .data$Z2, .data$Z3, .data$Z4, .data$Z5)
  df <- bind_cols(df, zscales)
  df <- df %>%
    mutate(
      Z1 = as.numeric(Z1),
      Z2 = as.numeric(Z2),
      Z3 = as.numeric(Z3),
      Z4 = as.numeric(Z4),
      Z5 = as.numeric(Z5)
    )
  return(df)
}


#' Visualization functions
#' 
#' @description 
#' `heatmap` makes a ggplot2 tile plot based on a tibble
#' 
#'
#' @param df, tibble with columns score, mutation_position, mutation info
#' @param x,
#' @param y,
#' @param fill,
#'   
#' @return heatmap object
#'
#' @examples heatmap_score(df,x, y, fill)
heatmap <- function(df, x="mutation_position", y="mutation", fill="pI") {
  heatmap <- ggplot(df, aes_string(x, y, fill = fill), color = "black") +
    scale_x_continuous(expand=c(0,0)) +
    scale_fill_continuous(type ="viridis", oob= scales::squish) +
    geom_tile() +
    labs(x="Sequence position", y="Residue mutated",color="Score") +
    theme_classic() +
    theme(
      panel.background = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      axis.line = element_line(colour = "black", size = 0),
      axis.text.x = element_text(face = "bold", color = "#000000"),
      axis.text.y = element_text(face = "bold", color = "#000000")
    )
  return(heatmap)
}

#' Visualization functions
#' 
#' @description 
#' `score_distribution` makes a ggplot2 tile plot based on a tibble
#' 
#'
#' @param df, tibble with columns score, mutation_position, mutation info
#' @param threshold, cutoff of scoring
#'   
#' @return heatmap object
#'
#' @examples heatmap_score(df, threshold)
score_distribution <- function(df, threshold) {
  distribution <- df %>%
    group_by(mutation) %>%
    filter(score > threshold) %>%
    ggplot(aes(mutation_position, score, color = mutation)) +
    geom_point() +
    theme_classic() +
    labs(x = "Residue mutated", y = "Score", title = "Scores by position over a determined threshold") +
    theme(
      legend.position = "bottom",
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      axis.line = element_line(colour = "black", size = 0),
      axis.text.x = element_text(face = "bold", color = "#000000"),
      axis.text.y = element_text(face = "bold", color = "#000000")
    ) +
    guides(fill = guide_legend(title = "Residue Mutated", nrow = 4), color=guide_legend(title = "Residue Mutated",nrow = 4)
    ) +
    facet_grid(aminoacid_class ~ .)
    return(distribution)
}


#' Visualization functions
#' 
#' @description 
#' `heatmap_score` makes a ggplot2 tile plot based on a tibble
#' 
#'
#' @param df, tibble with columns score, mutation_position, mutation info
#' @param score_min, numeric lower cutoff
#' @param score_max, numeric higher cutoff
#'
#' @return heatmap object
#'
#' @examples heatmap_score(df, 0.2, 4.3)
heatmap_score <- function(df, score_min, score_max) {
  heatmap <- ggplot(df, aes(mutation_position, mutation, fill = score), color = "black") +
    scale_x_continuous(expand=c(0,0)) +
    scale_fill_continuous(limits=c(score_min, score_max),type ="viridis", oob= scales::squish) +
    geom_tile() +
    labs(x="Sequence position", y="Residue mutated",color="Score") +
    theme_classic() +
    theme(
      panel.background = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      axis.line = element_line(colour = "black", size = 0),
      axis.text.x = element_text(face = "bold", color = "#000000"),
      axis.text.y = element_text(face = "bold", color = "#000000")
    )
  return(heatmap)
}


#' @description 
#' `heatmap_g_score` makes a ggplot2 tile plot based on a tibble
#' 
#'
#' @param df, tibble with x and y columns
#' @param score_min, numeric lower cutoff
#' @param score_max, numeric higher cutoff
#'
#' @return heatmap object
#'
#' @examples heatmap_score(df, score_min=0.2, score_max=4.3)
heatmap_g_score <- function(df, score_min, score_max) {
  heatmap_plot<- ggplot(df, aes(x=mutation_position,y=aminoacid_class, fill=score), color = "black") +
    scale_x_continuous(expand=c(0,0)) +
    scale_fill_continuous(limits=c(score_min, score_max),type ="viridis", oob= scales::squish) +
    geom_tile() +
    labs(x="Sequence position", y="Amino acid classification",color="Score", title="Score heatmap") + 
    theme_classic() +
    theme(
      panel.background = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      axis.line = element_line(colour = "black", size = 0),
      axis.text.x = element_text(face = "bold", color = "#000000"),
      axis.text.y = element_text(face = "bold", color = "#000000")
    )
  return(heatmap_plot)
}


#' @description 
#' `density_plot` makes a ggplot2 density plot based on a tibble
#' 
#' 
#' @param df, tibble with score column 
#'
#' @return density plot object
#'
#' @examples density_plot(df)
#' 
#' 
density_plot <- function(df, title) {
  density <- ggplot(df, aes(score)) +
    geom_density(alpha=0.4, fill="lightblue") +
    ggtitle(title) +
    labs(x="Score", y="Density",color="Score") +
    theme_classic() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      axis.line = element_line(colour = "black", size = 0),
      axis.text.x = element_text(face = "bold", color = "#000000"),
      axis.text.y = element_text(face = "bold", color = "#000000")) +
    guides(fill = guide_legend(ncol = 5))
  return(density)
}

#' @description 
#' `density_plot_residue` makes a ggplot2 density plot based on a tibble
#' 
#' 
#' @param df, tibble with score column and mutation columns
#'
#' @return density plot object
#'
#' @examples density_plot_residue(df)
density_plot_residue <- function(df){
  density_residue <- ggplot(df, aes(score, color=mutation, fill=mutation)) +
    geom_density(alpha=0.1) +
    theme_classic() +
    labs(x="Score", y="Density",title="Score density plof scores per residue mutated") +
    theme(
      legend.position = "bottom",
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      axis.line = element_line(colour = "black", size = 0),
      axis.text.x = element_text(face = "bold", color = "#000000"),
      axis.text.y = element_text(face = "bold", color = "#000000")
      ) +
    guides(fill = guide_legend(title = "Residue Mutated", nrow = 4), color=guide_legend(title = "Residue Mutated",nrow = 4)
    )
  return(density_residue)
}

#' @description 
#' `quick_sar` makes a ggplot2 line plot based on a tibble
#' 
#' @param df, tibble with score column and mutation_position columns
#'
#' @return density plot object
#'
#' @examples quick_sar(df)
quick_sar <- function(df, cutoff) {
  df_eff <- df %>%
    mutate(effective = case_when(score <= cutoff ~ 0,
                                 score > cutoff ~ 1)) %>%
    group_by(mutation_position) %>%
    summarise(N_eff = sum(effective))
  
  neff_plot <- df_eff %>%
    ggplot(aes(x=mutation_position, y = N_eff)) +
    geom_line() +
    theme_classic() +
    labs(x="Mutation position", y="Number of residues",title="Residues allowed per position in the sequence") +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      axis.line = element_line(colour = "black", size = 0),
      axis.text.x = element_text(face = "bold", color = "#000000"),
      axis.text.y = element_text(face = "bold", color = "#000000")
    )
  return(neff_plot)
}


#' sequence_logo
#'
#' @param df tibble 
#' @param start_position start of sequence window
#' @param end_position start of sequence window
#'
#' @return ggseq sequence plot
#' @export
#'
#' @examples
sequence_logo <- function(df, start_position, end_position) {
  
  # drop na and normalize score
  data_set <- df %>%
    drop_na() %>%
    mutate(score = (score - min(score)) / (max((score)-min(score))))
  
  # sum up score per mutation_position
  score_summary_per_position <- data_set %>%
    group_by(mutation_position) %>%
    summarise(score_sum=sum(score))
  
  # join score_summary_per_position with data_set and calculate score/score_sum for each position
  # replace any occuring na with 0
  # pivot to wider
  aa_tibble <- data_set %>%
    full_join(score_summary_per_position, by = "mutation_position") %>%
    filter(mutation_position >= start_position & mutation_position <= end_position) %>%
    mutate(score = score / score_sum) %>%
    select(mutation_position, mutation, score) %>%
    replace(is.na(.), 0) %>%
    pivot_wider(names_from = mutation_position, values_from = score)
  
  # prepping for ggseglogo
  aa_matrix = as.matrix(aa_tibble) #convert to matrix
  aa_column = aa_matrix[,1] # extract aa names
  labels = colnames(aa_matrix) # get mutation positions
  colnames(aa_matrix) <- NULL # remove column names
  row.names(aa_matrix) <- aa_column  # name rows with aa
  aa_matrix <- aa_matrix[,-1]  # remove aa columns in matrix area
  
  labels <- labels[-1] # remove mutation word
  
  window_length = end_position - start_position + 1 # calculate break range for x_scale
  
  # Generate sequence logo
  logo <- ggseqlogo(aa_matrix, method='custom', seq_type='aa')
  ggseq <- logo + 
   
    scale_x_continuous(breaks = 1:window_length, labels = labels) +
    theme_classic() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      axis.line = element_line(colour = "black", size = 0),
      axis.text.y = element_text(face = "bold", color = "#000000"),
      axis.text.x.bottom = element_text(vjust = 0.5, angle = -90, face = "bold", color = "#000000")) +
    xlab("Mutation position") + 
    ylab('Fraction of scoring')
  
  return(ggseq)
}


#' ANN2 - Artificial Neural Network
#'
#' @param df datset
#' @param folder folder to store output
#' @param name name of dataset
#' @param epochs 
#' @param hidden_layers 
#' @param scale descriptor scale
#' @param train_size 
#' @param seed_value 
#'
#' @return
#' @export
#'
#' @examples
ANN2 <- function(df, folder, name, epochs, hidden_layers,
                 scale, train_size, seed_value) {
  
  print(name)
  print(scale)
  
  graph_filename = paste(folder,"_ann_corr_",scale,"_",name,".png", sep="")
  obj_filename = paste(folder,"_ann_",scale,"_",name,".RData", sep="")
  RMSE_filename = paste(folder,"_ann_RMSE_",scale,"_",name,".tsv", sep="")
  
  data <- df
  
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
  Y_pred <- as.double(pred$predictions)
  res <- tibble(Y_pred, Y_test)
  
  RMSE_ <- rmse(res, Y_test[,1], Y_pred) %>%
    add_column(dataset = name,
             model = "ANN",
             scale = scale)
  
  corr_plot <- ggplot(res, aes(x=Y_pred, y=Y_test)) +
    geom_point() +
    theme_classic() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      axis.line = element_line(colour = "black", size = 0),
      axis.text.x = element_text(face = "bold", color = "#000000"),
      axis.text.y = element_text(face = "bold", color = "#000000")
    )
  
  ggsave(plot=corr_plot, graph_filename)
  
  write_ANN(NN, obj_filename)
  
  write_tsv(x = RMSE_, path = RMSE_filename)
  
  return() # could return object
}


glmnet_reg <- function (df, folder, name, alpha_, s_,
                        scale, train_size, seed_value) {
  print(name)
  print(scale)
  
  graph_filename = paste(folder,"_glmnet_corr_",scale,"_",name,".png", sep="")
  obj_filename = paste(folder,"_glmnet_",scale,"_",name,".Rdata", sep="")
  RMSE_filename = paste(folder,"_glmnet_RMSE_",scale,"_",name,".tsv", sep="")
  
  data <- df

  
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
  
  
  # Fit elasticnet model
  # ------------------------------------------------------------------------------
  fit.elasticnet <- glmnet(X_train, Y_train, type.measure="mse", alpha = alpha_, s = s_,
                           family="gaussian",
                           standardize=TRUE)
  
  # Preduct using elasticnet model
  # ------------------------------------------------------------------------------
  Y_pred_ElasticNet <- predict(fit.elasticnet, newx = X_test, s = s_, alpha=alpha_)
  results = tibble(Y_test, Y_pred_ElasticNet)
  RMSE_table = results
  
  
  # Pivot results and test response variables
  # ------------------------------------------------------------------------------
  results <- results %>%
    mutate(score = Y_test) %>%
    full_join(test, by = "score") %>%
    select(variant, score, Y_test, Y_pred_ElasticNet) %>%
    pivot_longer(c(-variant, -score), names_to = "model", values_to = "prediction")
  
  # Plot correlation between score and prediction
  # ------------------------------------------------------------------------------
  corr_plot <- ggplot(results, aes(x=score, y=prediction, color=model)) +
    geom_point() +
    theme_classic() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      axis.line = element_line(colour = "black", size = 0),
      axis.text.x = element_text(face = "bold", color = "#000000"),
      axis.text.y = element_text(face = "bold", color = "#000000")
    )
  
  ggsave(plot=corr_plot, graph_filename)
  
  
  RMSE_table <- RMSE_table %>%
    select(Y_test, Y_pred_ElasticNet)
    
  
  RMSE_ = rmse(RMSE_table, Y_test[,1], Y_pred_ElasticNet[,"1"]) %>%
    add_column(dataset = name,
               model = "glmnet",
               scale = scale)
  
  write_tsv(x = RMSE_, path = RMSE_filename)
  
  save(fit.elasticnet, file=obj_filename)
  
  return() # could return something
  
}



glmnet_CV <- function (df, name, folder, scale, train_size, seed_value, alpha_) {
  print(name)
  print(scale)
  
  graph_filename = paste(folder,"_plots/glmnet_CV_",scale,"_",name,".png", sep="")
  CV_filename = paste(folder,"_data/glmnet_CV_",scale,"_",name,".tsv", sep="")
  
  data <- df

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
  
  
  # Fit elasticnet model with CV
  # ------------------------------------------------------------------------------
  fit.elasticnet.cv <- cv.glmnet(X_train, Y_train, type.measure="mse", alpha = alpha_,
                               family="gaussian",
                               standardize=TRUE)
  
  
  CV_results <- tidy(fit.elasticnet.cv) %>%
  mutate(log_lambda = log10(lambda))
  
  
  CV_plot <- ggplot(CV_results, aes(x=log_lambda, y=estimate)) +
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error)) +
  geom_point() +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(colour = "black", size = 0),
    axis.text.x = element_text(face = "bold", color = "#000000"),
    axis.text.y = element_text(face = "bold", color = "#000000")
  )
  
  ggsave(plot=CV_plot, graph_filename)
  write_tsv(x = CV_results, path = CV_filename)
  
  return() # clould return something
}