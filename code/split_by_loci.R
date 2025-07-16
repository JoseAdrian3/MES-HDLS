library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(biomaRt)

cat("Splitting by locis...\n")

# Reproducibility
set.seed(123)

# Line command arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("You must provide one the path to the train file")
}
train_path <- args[1]

# Data load
train_df <- read.csv(train_path, row.names = 1)

train_samples <- colnames(train_df)

# Transpose the training data so rows represent samples and columns represent probes
train_df <- as.data.frame(t(train_df))

# Probe illumina annotation
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Probe ids
probe_ids <- colnames(train_df)

# Columns with Illumina annotation for locis
probe_info <- anno[rownames(anno) %in% probe_ids, c("UCSC_RefGene_Group", "Regulatory_Feature_Group", "UCSC_RefGene_Name")]

# TSS filter
is_in_tss <- function(group) {
  terms <- unlist(strsplit(group, ";"))
  any(terms %in% c("TSS200", "TSS1500"))
}

tss_probes <- rownames(probe_info)[sapply(probe_info$UCSC_RefGene_Group, is_in_tss)]

# Promotor probes
regulatory_criteria <- c("Promoter_Associated", "Promoter_Associated_Cell_type_specific")
promoter_associated_probes <- rownames(probe_info)[grep(paste(regulatory_criteria, collapse = "|"), probe_info$Regulatory_Feature_Group)]

# Combine TSS probes and promoter probes to create a unique set of promoter probes
promoter_probes_names <- unique(c(tss_probes, promoter_associated_probes))
promoter_probes <- probe_info[rownames(probe_info) %in% promoter_probes_names, ]

# Protein coding probes
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
protein_coding_genes <- getBM(attributes = c("hgnc_symbol"),
                              filters = "biotype",
                              values = "protein_coding",
                              mart = ensembl)
protein_coding_genes_list <- protein_coding_genes$hgnc_symbol

is_protein_coding <- function(gene_list) {
  genes <- unlist(strsplit(gene_list, ";"))
  any(genes %in% protein_coding_genes_list)
}

protein_coding_probes <- probe_info[sapply(probe_info$UCSC_RefGene_Name, is_protein_coding), ]

# Probes by loci regions
unique_values <- unique(unlist(strsplit(as.character(probe_info$UCSC_RefGene_Group), ";")))

lista_subsets <- lapply(unique_values, function(valor) {
  indices_coinciden <- sapply(
    strsplit(as.character(probe_info$UCSC_RefGene_Group), ";"),
    function(x) valor %in% x
  )
  subset_actual <- probe_info[indices_coinciden, ]
  return(subset_actual)
})
names(lista_subsets) <- unique_values

# Directorio principal para loci
output_dir <- "../data/loci"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Dirs for every loci
loci_dirs <- list(
  promoter = file.path(output_dir, "promoter"),
  protein_coding = file.path(output_dir, "protein_coding"),
  tss200 = file.path(output_dir, "tss200"),
  body = file.path(output_dir, "body"),
  stExon = file.path(output_dir, "stExon"),
  utr5 = file.path(output_dir, "utr5"),
  tss1500 = file.path(output_dir, "tss1500"),
  utr3 = file.path(output_dir, "utr3"),
  exonbnd = file.path(output_dir, "exonbnd")
)

# Create every loci dir
for(dir in loci_dirs) {
  if(!dir.exists(dir)) {
    dir.create(dir)
  }
}

# Promoter
promoter_limma_df_train <- train_df[, colnames(train_df) %in% rownames(promoter_probes)]
saveRDS(promoter_limma_df_train, file = file.path(loci_dirs$promoter, "promoter_train.rds"))

# Protein coding
protein_coding_limma_df_train <- train_df[, colnames(train_df) %in% rownames(protein_coding_probes)]
saveRDS(protein_coding_limma_df_train, file = file.path(loci_dirs$protein_coding, "protein_coding_train.rds"))

# TSS200
tss200_limma_df_train <- train_df[, colnames(train_df) %in% rownames(lista_subsets$TSS200)]
saveRDS(tss200_limma_df_train, file = file.path(loci_dirs$tss200, "tss200_train.rds"))

# Body
body_limma_df_train <- train_df[, colnames(train_df) %in% rownames(lista_subsets$Body)]
saveRDS(body_limma_df_train, file = file.path(loci_dirs$body, "body_train.rds"))

# stExon
stexon_limma_df_train <- train_df[, colnames(train_df) %in% rownames(lista_subsets$`1stExon`)]
saveRDS(stexon_limma_df_train, file = file.path(loci_dirs$stExon, "stExon_train.rds"))

# 5'UTR
utr5_limma_df_train <- train_df[, colnames(train_df) %in% rownames(lista_subsets$`5'UTR`)]
saveRDS(utr5_limma_df_train, file = file.path(loci_dirs$utr5, "utr5_train.rds"))

# TSS1500
tss1500_limma_df_train <- train_df[, colnames(train_df) %in% rownames(lista_subsets$TSS1500)]
saveRDS(tss1500_limma_df_train, file = file.path(loci_dirs$tss1500, "tss1500_train.rds"))

# 3'UTR
utr3_limma_df_train <- train_df[, colnames(train_df) %in% rownames(lista_subsets$`3'UTR`)]
saveRDS(utr3_limma_df_train, file = file.path(loci_dirs$utr3, "utr3_train.rds"))

# ExonBnd
exonbnd_limma_df_train <- train_df[, colnames(train_df) %in% rownames(lista_subsets$ExonBnd)]
saveRDS(exonbnd_limma_df_train, file = file.path(loci_dirs$exonbnd, "exonbnd_train.rds"))

cat("Division by locis done successfully!\n")
