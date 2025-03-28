library(dplyr)

cat("Performing differencial expression analysis for every loci...\n")

# Command-line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 ) {
  stop("Please provide the test dataset path and the train sheet as an argument.")
}
# Set path to test and sheet datasets and load them
test_path <- args[1]
sheet_path <- args[2]

test_df <- read.csv(test_path, row.names = 1)
test_df <- t(test_df)
sheet_df <- read.csv(sheet_path, row.names = 1)

# Main directory containing loci data.
data_dir <- "./data/loci"

# Load train datasets
promoter_train <- read.csv(file.path(data_dir, "promoter", "promoter_train.csv"), row.names = 1)
protein_coding_train <- read.csv(file.path(data_dir, "protein_coding", "protein_coding_train.csv"), row.names = 1)
tss200_train <- read.csv(file.path(data_dir, "tss200", "tss200_train.csv"), row.names = 1)
body_train <- read.csv(file.path(data_dir, "body", "body_train.csv"), row.names = 1)
stExon_train <- read.csv(file.path(data_dir, "stExon", "stExon_train.csv"), row.names = 1)
utr5_train <- read.csv(file.path(data_dir, "utr5", "utr5_train.csv"), row.names = 1)
tss1500_train <- read.csv(file.path(data_dir, "tss1500", "tss1500_train.csv"), row.names = 1)
utr3_train <- read.csv(file.path(data_dir, "utr3", "utr3_train.csv"), row.names = 1)
exonbnd_train <- read.csv(file.path(data_dir, "exonbnd", "exonbnd_train.csv"), row.names = 1)

# Function to perform differential expression analysis using the limma package.
limma_analysis <- function(train_df, sample_sheet_train, output_path) {
  
  # Convert sample information to appropriate data types
  f       <- factor(sample_sheet_train$Group)               # Sample groups as a factor.
  PMD     <- as.numeric(sample_sheet_train$PMD)             # Numeric covariate: PMD.
  Age     <- as.numeric(sample_sheet_train$Age)             # Numeric covariate: Age.
  Sex     <- as.factor(sample_sheet_train$Sex)              # Categorical covariate: Sex.
  NeuNP   <- as.numeric(sample_sheet_train$NeuNP)           # Numeric covariate: NeuNP.
  DoubleN <- as.numeric(sample_sheet_train$DoubleN)         # Numeric covariate: DoubleN.
  Sox10P  <- as.numeric(sample_sheet_train$Sox10P)          # Numeric covariate: Sox10P.
  Batch1  <- as.factor(sample_sheet_train$Sentrix_ID)       # Batch effect as a factor.
  Batch2  <- as.factor(sample_sheet_train$Sentrix_Position) # Secondary batch effect as a factor.
  
  # Construct the design matrix for the linear model
  design <- model.matrix(~ 0 + f + Age + Sex + NeuNP + DoubleN + PMD + Batch1 + Batch2)
  
  # Ensure all data in the train set is numeric
  train_df[] <- lapply(train_df, function(x) as.numeric(as.character(x)))
  
  # Transpose the train data to have probes as rows and samples as columns
  t_train_data_limma <- t(train_df)
  
  # Fit linear model to train data
  fit <- limma::lmFit(t_train_data_limma, design)
  
  # Define contrasts for the comparisons
  # This part MUST be modified manually if you change a dataset with other comparisons
  cm <- limma::makeContrasts(
    DLBvsCTRL    = fDLB - fCTRL,
    PDvsCTRL     = fPD - fCTRL,
    PDDvsCTRL    = fPDD - fCTRL,
    neurovsCTRL  = ((fDLB + fPD + fPDD) / 3) - fCTRL,
    levels = design
  )
  
  # Apply contrasts to fitted model.
  fit2 <- limma::contrasts.fit(fit, cm)
  
  # Compute empirical Bayes statistics to improve variance estimates.
  fit3 <- limma::eBayes(fit2)
  
  # Determine significant results based on the tests.
  resultsfit <- limma::decideTests(fit3)
  summary_results <- summary(resultsfit)
  
  # Return and save results
  result <- list(
    summary_fit = summary_results,
    final_fit   = fit3
  )
  
  saveRDS(result, output_path)
  
  return(result)
}

# Function to extract probes with significant differential expression.
get_significant_probes <- function(pvalue_vector, threshold = 2) {
  
  # Ngative log10 of p-values.
  neg_log10_pvals <- -log10(pvalue_vector)
  
  # Indices where transformed p-value exceeds threshold.
  sig_idx <- neg_log10_pvals > threshold
  
  result <- data.frame(
    probe = names(pvalue_vector)[sig_idx],
    neg_log10_pval = neg_log10_pvals[sig_idx],
    row.names = NULL
  )
  
  return(result)
}

# Perform differential expression analysis for each loci
result_promoter <- limma_analysis(promoter_train, sheet_df, file.path(data_dir, "promoter", "diff_exp_promoter.rds"))
result_protein_coding <- limma_analysis(protein_coding_train, sheet_df, file.path(data_dir, "protein_coding", "diff_exp_protein_coding.rds"))
result_tss200 <- limma_analysis(tss200_train, sheet_df, file.path(data_dir, "tss200", "diff_exp_tss200.rds"))
result_body <- limma_analysis(body_train, sheet_df, file.path(data_dir, "body", "diff_exp_body.rds"))
result_stExon <- limma_analysis(stExon_train, sheet_df, file.path(data_dir, "stExon", "diff_exp_stExon.rds"))
result_utr5 <- limma_analysis(utr5_train, sheet_df, file.path(data_dir, "utr5", "diff_exp_utr5.rds"))
result_tss1500 <- limma_analysis(tss1500_train, sheet_df, file.path(data_dir, "tss1500", "diff_exp_tss1500.rds"))
result_utr3 <- limma_analysis(utr3_train, sheet_df, file.path(data_dir, "utr3", "diff_exp_utr3.rds"))
result_exonbnd <- limma_analysis(exonbnd_train, sheet_df, file.path(data_dir, "exonbnd", "diff_exp_exonbnd.rds"))

# loci and constrats to analyze diff exp
loci <- c("promoter", "protein_coding", "tss200", "body", "stExon", "utr5", "tss1500", "utr3", "exonbnd")
contrasts <- c("DLBvsCTRL", "PDvsCTRL", "PDDvsCTRL", "neurovsCTRL")

# Iterate for every loci
for (locis in loci) {
  
  # Get loci train data
  train_data <- get(paste0(locis, "_train"))
  
  # Get result data
  result <- get(paste0("result_", locis))
  
  # Iterate for every contrast
  for (contrast in contrasts) {
    
    # Create contrast dir
    out_dir <- file.path(data_dir, locis, contrast)
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
    }
    
    # Get current significant probes for contrast and loci
    sig_probes <- get_significant_probes(result$final_fit$p.value[, contrast], threshold = 2)
    
    # Filter train and set data for significant probes
    train_diff <- train_data[, sig_probes$probe, drop = FALSE]
    test_diff  <- test_df[, sig_probes$probe, drop = FALSE]
    
    write.csv(train_diff,
              file = file.path(out_dir, paste0(locis, "_train_diff_exp_", contrast, ".csv")),
              row.names = TRUE)
    write.csv(test_diff,
              file = file.path(out_dir, paste0(locis, "_test_diff_exp_", contrast, ".csv")),
              row.names = TRUE)
  }
}

cat("Differencial expression analysis done succesfully!\n")