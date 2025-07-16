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
data_dir <- "../data/loci"

# Load train datasets from .rds files
promoter_train <- readRDS(file.path(data_dir, "promoter", "promoter_train.rds"))
protein_coding_train <- readRDS(file.path(data_dir, "protein_coding", "protein_coding_train.rds"))
tss200_train <- readRDS(file.path(data_dir, "tss200", "tss200_train.rds"))
body_train <- readRDS(file.path(data_dir, "body", "body_train.rds"))
stExon_train <- readRDS(file.path(data_dir, "stExon", "stExon_train.rds"))
utr5_train <- readRDS(file.path(data_dir, "utr5", "utr5_train.rds"))
tss1500_train <- readRDS(file.path(data_dir, "tss1500", "tss1500_train.rds"))
utr3_train <- readRDS(file.path(data_dir, "utr3", "utr3_train.rds"))
exonbnd_train <- readRDS(file.path(data_dir, "exonbnd", "exonbnd_train.rds"))


limma_analysis <- function(train_df, sample_sheet_train, contrast_name, output_path) {
  # Mapeo explícito de nombres de contraste a fórmulas válidas de limma
  contrast_map <- list(
    DLBvsCTRL   = "fDLB - fCTRL",
    PDvsCTRL    = "fPD - fCTRL",
    PDDvsCTRL   = "fPDD - fCTRL",
    neurovsCTRL = "(fDLB + fPD + fPDD)/3 - fCTRL"
  )
  
  # Verificar que el contraste exista
  if (!(contrast_name %in% names(contrast_map))) {
    stop(paste0("Contraste no reconocido: ", contrast_name))
  }
  
  # Convertir columnas
  f       <- factor(sample_sheet_train$Group)
  PMD     <- as.numeric(sample_sheet_train$PMD)
  Age     <- as.numeric(sample_sheet_train$Age)
  Sex     <- as.factor(sample_sheet_train$Sex)
  NeuNP   <- as.numeric(sample_sheet_train$NeuNP)
  DoubleN <- as.numeric(sample_sheet_train$DoubleN)
  Sox10P  <- as.numeric(sample_sheet_train$Sox10P)
  Batch1  <- as.factor(sample_sheet_train$Sentrix_ID)
  Batch2  <- as.factor(sample_sheet_train$Sentrix_Position)
  
  # Diseño del modelo
  design <- model.matrix(~ 0 + f + Age + Sex + NeuNP + DoubleN + PMD + Batch1 + Batch2)
  
  # Asegurar formato numérico
  train_df[] <- lapply(train_df, function(x) as.numeric(as.character(x)))
  t_train_data_limma <- t(train_df)
  
  # Ajustar modelo
  fit <- limma::lmFit(t_train_data_limma, design)
  
  # Crear contraste específico
  contrast_formula <- contrast_map[[contrast_name]]
  cm <- limma::makeContrasts(contrasts = contrast_formula, levels = design)
  
  fit2 <- limma::contrasts.fit(fit, cm)
  fit3 <- limma::eBayes(fit2)
  
  resultsfit <- limma::decideTests(fit3)
  summary_results <- summary(resultsfit)
  
  result <- list(
    summary_fit = summary_results,
    final_fit   = fit3,
    contrast_used = contrast_formula
  )
  
  saveRDS(result, output_path)
  return(result)
}

filter_by_contrast <- function(train_df, sample_sheet, test_df, contrast) {
  if (contrast == "neurovsCTRL") {
    # No filtrar: usar todo
    return(list(train = train_df, sheet = sample_sheet, test = test_df))
  }
  
  contrast_parts <- unlist(strsplit(contrast, "vs"))
  groups <- gsub("f", "", contrast_parts)
  
  filtered_sheet <- sample_sheet[sample_sheet$Group %in% groups, ]
  filtered_train <- train_df[rownames(filtered_sheet), , drop = FALSE]
  
  test_samples_to_keep <- rownames(test_df)[grepl(paste(groups, collapse = "|"), rownames(test_df))]
  filtered_test <- test_df[test_samples_to_keep, , drop = FALSE]
  
  return(list(train = filtered_train, sheet = filtered_sheet, test = filtered_test))
}

# Lista de datasets
loci_datasets <- list(
  promoter = promoter_train,
  protein_coding = protein_coding_train,
  tss200 = tss200_train,
  body = body_train,
  stExon = stExon_train,
  utr5 = utr5_train,
  tss1500 = tss1500_train,
  utr3 = utr3_train,
  exonbnd = exonbnd_train
)

loci <- names(loci_datasets)
contrasts <- c("DLBvsCTRL", "PDvsCTRL", "PDDvsCTRL", "neurovsCTRL")

for (locis in loci) {
  for (contrast in contrasts) {
    
    filtered <- filter_by_contrast(loci_datasets[[locis]], sheet_df, test_df, contrast)
    
    out_path <- file.path(data_dir, locis, paste0("diff_exp_", locis, "_", contrast, ".rds"))
    
    result <- limma_analysis(filtered$train, filtered$sheet, contrast, out_path)
    
    real_contrast <- result$contrast_used
    
    sig_probes <- get_significant_probes(result$final_fit$p.value[, real_contrast], threshold = 2)
    
    train_diff <- filtered$train[, sig_probes$probe, drop = FALSE]
    test_diff  <- filtered$test[, sig_probes$probe, drop = FALSE]
    
    out_dir <- file.path(data_dir, locis, contrast)
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
    }
    
    write.csv(train_diff, file = file.path(out_dir, paste0(locis, "_train_diff_exp_", contrast, ".csv")), row.names = TRUE)
    write.csv(test_diff,  file = file.path(out_dir, paste0(locis, "_test_diff_exp_", contrast, ".csv")), row.names = TRUE)
  }
}


cat("Differencial expression analysis done succesfully!\n")
