library(data.table)
library(tidyr)
library(dplyr)
source("XCMAX4.R")

args <- commandArgs(trailingOnly = TRUE)

chrX_file_path <- args[1]
pheno_file_path <- args[2]
output_name <- args[3]

chrX <- fread(chrX_file_path, header=T, data.table = F)
pheno <- read.table(pheno_file_path, header=T)
pheno$gender <- ifelse(pheno$gender == "Male", 0, 1)

results_df <- data.frame(CHR = character(),
                         BP = numeric(),
                         A1 = character(),
                         A2 = character(),
                         Beta = numeric(),
                         P = numeric(),
                         SE = numeric(),
                         stringsAsFactors = FALSE)

# Loop through all variants in the chrX dataframe
for (i in 2:ncol(chrX)) {
  
  # Subset the 1st and ith column from df chrX
  chrX_subset <- chrX[, c(1, i)]
  
  # Merge the chrX_subset with the pheno dataframe by the first column (ID)
  merged_data <- merge(chrX_subset, pheno, by = "ID")
  
  # Rearrange the columns to "FA", genotype, and "gender"
  final_data <- merged_data[, c("FA", colnames(merged_data)[2], "gender", setdiff(names(merged_data), c("FA", colnames(merged_data)[2], "gender", "ID")))]
  
  # Run analysis
  res <- XCMAX4(final_data)
  
  split_names <- strsplit(colnames(merged_data)[2], ":")[[1]]

  results_df <- rbind(results_df, data.frame(
    CHR = split_names[1],
    BP = split_names[2],
    A1 = split_names[4], # alt allele
    A2 = split_names[3], # ref allele
    Beta = res$statictic,
    P = res$`p-value`
  ))
}

# Filter SNPs with missing statistics (results == NA) and sort using p-value
results_df <- results_df[!is.na(results_df$P), ] %>%
  arrange(P)

write.table(results_df, output_name, quote=F, row.names = F, sep = "\t")