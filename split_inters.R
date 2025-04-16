# Load the RDS file
data <- readRDS("/home/projects2/kam_project/outputs/checkpoints_r5/full_interactions.rds")

# Filter rows where values in columns A and B are equal
no_self_inter <- subset(data, source_cell_type == target_cell_type)
only_self_inter <- subset(data, source_cell_type != target_cell_type)

# Save the filtered data to a new RDS file
saveRDS(no_self_inter, "/home/projects2/kam_project/outputs/checkpoints_r5/filter_intneractions-no_self_inter.rds")
saveRDS(only_self_inter, "/home/projects2/kam_project/outputs/checkpoints_r5/filter_interactions-only_self_inter.rds")
