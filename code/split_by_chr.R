library(dplyr)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

cat("Splitting by chrs...\n")

# Command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript chr_sep.R <input_file_path> <output_directory>")
}

# Exec Arguments 
file_path <- normalizePath(args[1])
out_directory <- normalizePath(args[2])

print(paste("Reading file from:", file_path))
print(paste("Saving files to:", out_directory))

if (!file.exists(file_path)) {
  stop(paste("Error: File", file_path, "does not exist. Check the path."))
}

# Probes dataset
limma_df <- read.csv(file_path, row.names = 1)

# Probes are supposed to be, at first, in the rows
limma_df_var_lim_t <- as.data.frame(limma_df)

# Probe annotations
annotation <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Probes map to chromosomes
probe_to_chr <- annotation$chr
names(probe_to_chr) <- rownames(annotation)

# Probe names
probe_names <- rownames(limma_df_var_lim_t)

# Probe names to chromosome names
chromosome_names <- probe_to_chr[probe_names]

# Chromosome column
limma_df_var_lim_t$chr <- chromosome_names

# All probes,in theory if all types are present, there are from chr1 to chr22.
chrom_list <- unique(limma_df_var_lim_t$chr)

# Directories created and data saved by chromosome
for (chr in chrom_list) {
  dir_path <- file.path(out_directory, chr)
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  
  sub_df <- limma_df_var_lim_t[limma_df_var_lim_t$chr == chr, ]
  
  write.csv(sub_df, file = file.path(dir_path, paste0("limma_df_", chr, ".csv")), row.names = TRUE)
}

cat("Division by chrs done succesfully!\n")