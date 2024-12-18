library(dplyr)
library(data.table)

# Read the summary statistics and select columns (MarkerID, AF, BETA, SE and P-value)
female_stats <- fread('PATH_TO_FEMALE_SUMMARY_STATISTICS', sep = "\t", header = TRUE, select = c(3,7,9,10,13))
male_stats <- fread('PATH_TO_MALE_SUMMARY_STATISTICS', sep = "\t", header = TRUE, select = c(3,7,9,10,13))

# Rename columns
colnames(female_stats) <- c("MarkerID", "AF_female", "BETA_female", "SE_female", "P_female")
colnames(male_stats) <- c("MarkerID", "AF_male", "BETA_male", "SE_male", "P_male")

# Merge the data on the 'SNP' column
merged_stats <- merge(female_stats, male_stats, by = "MarkerID", all = TRUE)

# Remove individual summary statistics to save memory space
rm(female_stats, male_stats)

# Define a function to calculate the Z-score and p-value
compare_betas <- function(beta_female, se_female, beta_male, se_male) {
  # Return NA if either beta or SE is missing for one group
  if (is.na(beta_female) || is.na(beta_male) || is.na(se_female) || is.na(se_male)) {
    return(c(NA, NA))
  } else {
    # Z-score calculation
    z_score <- (beta_female - beta_male) / sqrt(se_female^2 + se_male^2)
    # Two-tailed p-value from Z-score
    p_value <- 2 * (1 - pnorm(abs(z_score)))
    return(c(z_score, p_value))
  }
}

# Apply the function to calculate Z-scores and p-values for each row in the merged data
results <- t(apply(merged_stats, 1, function(row) {
  compare_betas(as.numeric(row["BETA_female"]), as.numeric(row["SE_female"]),
                as.numeric(row["BETA_male"]), as.numeric(row["SE_male"]))
}))

# Split the results into Z-scores and p-values
merged_stats$Z_score <- results[, 1]  # First column is the Z-score
merged_stats$p_value <- results[, 2]  # Second column is the p-value

# Sort the merged_stats dataframe by the p_value column in ascending order
merged_stats <- merged_stats %>%
  arrange(p_value)

# FUMA format
merged_stats_FUMA <- merged_stats
split_marker <- strsplit(merged_stats_FUMA$MarkerID, ":")

# Convert list to a data frame
split_df <- do.call(rbind, split_marker)

# Combine with original data frame
merged_stats_FUMA <- cbind(merged_stats_FUMA[, c("p_value")], split_df)
colnames(merged_stats_FUMA) <- c("p_value", "CHR", "BP", "A2", "A1")

# Save the final output to a CSV file
write.csv(merged_stats, 'interaction_analysis_results.csv', row.names = FALSE, quote=F)
