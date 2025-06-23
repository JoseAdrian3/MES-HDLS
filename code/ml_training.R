root_dir        <- "../data/loci"
prob_threshold  <- 0.5
model_type      <- "svm_linear"
rf_always_split <- c()
set.seed(123)

library(caret)
library(dplyr)
library(ranger)

# ML training function
ml_training <- function(region_dir, comparison_dir, root_dir,
                           threshold = 0.5, type = "svm_radial",
                           always_split = character()) {
  
  # Positive label
  positive_label <- sub("vsCTRL$", "", comparison_dir)

  # Train and test location
  train_file <- file.path(root_dir, region_dir, comparison_dir,
                          sprintf("%s_train_diff_exp_%s.csv", region_dir, comparison_dir))
  test_file  <- file.path(root_dir, region_dir, comparison_dir,
                          sprintf("%s_test_diff_exp_%s.csv",  region_dir, comparison_dir))

  # Load train and test files
  train_df <- tryCatch({
    df <- read.csv(train_file, row.names = 1, check.names = FALSE)
    if (ncol(df) < 2 || nrow(df) == 0) stop("Invalid train file.")
    df
  }, error = function(e) {
    message("[ERROR] Reading train file: ", train_file)
    return(NULL)
  })

  test_df <- tryCatch({
    df <- read.csv(test_file, row.names = 1, check.names = FALSE)
    if (ncol(df) < 2 || nrow(df) == 0) stop("Invalid test file.")
    df
  }, error = function(e) {
    message("[ERROR] Reading test file: ", test_file)
    return(NULL)
  })

  # Here we subsample the train and test datasets to train with the correct disease
  if (positive_label == "neuro") {
    train_df$Sample_Group <- ifelse(grepl("^CTRL", rownames(train_df)), "CTRL", "neuro")
    test_df$Sample_Group  <- ifelse(grepl("^CTRL", rownames(test_df)),  "CTRL", "neuro")
  } else {
    keep_regex <- paste0("^(CTRL|", positive_label, ")")
    train_df <- train_df[grepl(keep_regex, rownames(train_df)), , drop = FALSE]
    test_df  <- test_df [grepl(keep_regex, rownames(test_df )), , drop = FALSE]

    train_df$Sample_Group <- ifelse(grepl(paste0("^", positive_label), rownames(train_df)),
                                    positive_label, "CTRL")
    test_df$Sample_Group  <- ifelse(grepl(paste0("^", positive_label), rownames(test_df)),
                                    positive_label, "CTRL")
  }

  if (is.null(train_df) || is.null(test_df)) return(NULL)

  levels_vec <- c("CTRL", positive_label)
  train_df$Sample_Group <- factor(train_df$Sample_Group, levels = levels_vec)
  test_df$Sample_Group  <- factor(test_df$Sample_Group,  levels = levels_vec)

  # Training control
  tr_control <- trainControl(method = "repeatedcv",
                             number = 5, repeats = 3,
                             classProbs = TRUE,
                             savePredictions = "final",
                             summaryFunction = twoClassSummary)

  # Train model
  model_fit <- tryCatch({
    if (type == "svm_linear") {

      train(Sample_Group ~ .,
            data       = train_df,
            method     = "svmLinear",
            trControl  = tr_control,
            tuneLength = 5)

    } else if (type == "svm_radial") {

      train(Sample_Group ~ .,
            data       = train_df,
            method     = "svmRadial",
            trControl  = tr_control,
            tuneLength = 5)

    } else if (type == "random_forest") {

      valid_split_vars <- intersect(always_split, setdiff(colnames(train_df), "Sample_Group"))

      train(Sample_Group ~ .,
            data       = train_df,
            method     = "ranger",
            trControl  = tr_control,
            tuneLength = 5,
            always.split.variables = valid_split_vars)

    } else {
      stop("Unknown model type: ", type)
    }
  }, error = function(e) {
    message("[ERROR] Training failed in ", region_dir, "|", comparison_dir)
    return(NULL)
  })

  if (is.null(model_fit)) return(NULL)

  # Predictions
  pred_prob  <- predict(model_fit, newdata = test_df, type = "prob")
  pred_class <- ifelse(pred_prob[[positive_label]] > threshold, positive_label, "CTRL")
  pred_class <- factor(pred_class, levels = levels_vec)

  cm <- confusionMatrix(pred_class, test_df$Sample_Group, positive = positive_label)

  data.frame(
    comparison        = comparison_dir,
    region            = region_dir,
    p_value           = unname(cm$overall["AccuracyPValue"]),
    mcnemar_p_value   = unname(cm$overall["McnemarPValue"]),
    NIR               = unname(cm$overall["AccuracyNull"]),
    accuracy          = unname(cm$overall["Accuracy"]),
    balanced_accuracy = unname(cm$byClass["Balanced Accuracy"]),
    kappa             = unname(cm$overall["Kappa"]),
    stringsAsFactors  = FALSE
  )
}

# Results dataframe with all stadistics
results_df <- data.frame()

regions <- list.dirs(root_dir, full.names = FALSE, recursive = FALSE)

# Locis loop
for (region in regions) {
  comparison_dirs <- list.dirs(file.path(root_dir, region), full.names = FALSE, recursive = FALSE)
  # Disease loop
  for (comp in comparison_dirs) {
    if (!grepl("vsCTRL$", comp)) next
    message("Processing: ", region, ";", comp)

    row_res <- ml_training(region, comp, root_dir,
                              threshold    = prob_threshold,
                              type         = model_type,
                              always_split = rf_always_split)
    if (!is.null(row_res)) {
      results_df <- bind_rows(results_df, row_res)
    }
  }
}

outfile <- sprintf("%s_results_all_regions.csv", model_type)
write.csv(results_df, outfile, row.names = FALSE)
print(head(results_df))

print("\nML completed. Results saved in: ", outfile)
