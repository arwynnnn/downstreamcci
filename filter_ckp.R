library(dplyr)
library(purrr)

# Load the full interaction details
full_interactions <- readRDS("/home/projects2/kam_project/outputs/full_interactions.rds")

# Set your custom threshold
required_fraction <- 0.2  # 20%
ratio_cutoff <- 1         # Ratio must be >1 to count

# Function to compute fraction of ratios >1 for a list of details
compute_fraction_above_cutoff <- function(details_list, cutoff = 1) {
  if (nrow(details_list) == 0) {
    return(0)
  }
  mean(details_list$ratio > cutoff)
}

# Apply to source and target
full_interactions <- full_interactions %>%
  mutate(
    frac_source_above_1 = purrr::map_dbl(gene_set_results_source, compute_fraction_above_cutoff),
    frac_target_above_1 = purrr::map_dbl(gene_set_results_target, compute_fraction_above_cutoff)
  )

# Now apply your filtering
final_interactions_custom <- full_interactions %>%
  filter(frac_source_above_1 >= required_fraction,
         frac_target_above_1 >= required_fraction) %>%
  arrange(desc(composite_score))  # Sort if you want

# Save the newly filtered interactions
saveRDS(final_interactions_custom, file = "/home/projects2/kam_project/outputs/final_interactions_20q.rds")
