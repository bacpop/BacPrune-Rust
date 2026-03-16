# Set seed for reproducibility
set.seed(123)

# Generate the initial matrix of 1s and 0s
n_rows <- 10
n_cols <- 100
matrix_data <- matrix(sample(c(0, 1), n_rows * n_cols, replace = TRUE), nrow = n_rows, ncol = n_cols)

# Create 3 sets of identical columns
set_size <- 5  # Number of columns in each set
n_sets <- 3  # Number of sets

# Randomly choose columns for the identical sets
chosen_columns <- sample(1:n_cols, n_sets * set_size)

# Generate the identical sets
for (i in 1:n_sets) {
  start_idx <- (i - 1) * set_size + 1
  end_idx <- i * set_size
  columns <- chosen_columns[start_idx:end_idx]
  matrix_data[, columns] <- matrix_data[, columns[1]]
}

# Assign column names
colnames(matrix_data) <- as.character(1:n_cols)

# Function to find groups of identical columns
group_identical_columns <- function(mat) {
  column_groups <- list()
  visited <- rep(FALSE, ncol(mat))
  
  group_index <- 1
  
  for (i in 1:ncol(mat)) {
    if (!visited[i]) {
      identical_cols <- which(apply(mat, 2, function(col) all(col == mat[, i])))
      if (length(identical_cols) > 1) {
        column_groups[[group_index]] <- colnames(mat)[identical_cols]
        visited[identical_cols] <- TRUE
        group_index <- group_index + 1
      }
    }
  }
  
  column_groups
}

duplicated.matrix(t(matrix_data))

matrix_data[!duplicated(as.list(matrix_data))]
colSums(matrix_data)/10

# Find groups of identical columns
identical_column_groups <- group_identical_columns(matrix_data)

# Find the maximum length of the groups
max_length <- max(sapply(identical_column_groups, length))

# Convert the list of groups to a more readable format with padding
grouped_matrix <- do.call(rbind, lapply(seq_along(identical_column_groups), function(idx) {
  c(paste("Group", idx), identical_column_groups[[idx]], rep(NA, max_length - length(identical_column_groups[[idx]])))
}))

# Set column names for the grouped matrix
colnames(grouped_matrix) <- c("Group", paste0("Column", seq(max_length)))

na.omit(grouped_matrix)
sapply(grouped_matrix, as.numeric)


# Output the result
print(grouped_matrix)

write.csv(x = matrix_data, file = "/nfs/research/jlees/jacqueline/gwas_code/ld_pruning/test_data_new/test_genotypes.csv", col.names = T, row.names = F)
