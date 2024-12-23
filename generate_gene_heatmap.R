
# Define the function for creating the heatmap
generate_gene_heatmap <- function(gene_lists) {
  # Clean the gene names by removing "up" and "down" suffixes
  cleaned_gene_lists <- lapply(gene_lists, function(x) gsub("up|down", "", x))

  # Create a vector of all unique genes across all lists (after cleaning)
  unique_genes <- unique(unlist(cleaned_gene_lists))

  # Initialize the matrix with empty strings, representing no gene presence
  gene_matrix <- matrix("", nrow = length(gene_lists), ncol = length(unique_genes),
                        dimnames = list(names(gene_lists), unique_genes))

  # Fill the matrix with "up" (blue), "down" (red), and "" (empty) for missing genes
  for (i in 1:length(gene_lists)) {
    for (gene in gene_lists[[i]]) {
      gene_name <- gsub("up|down", "", gene)  # Remove "up" or "down" to get the gene name
      if (grepl("up", gene)) {
        gene_matrix[i, gene_name] <- "up"
      } else if (grepl("down", gene)) {
        gene_matrix[i, gene_name] <- "down"
      }
    }
  }

  # Convert the matrix to a numeric matrix for pheatmap
  # "up" will be 1, "down" will be 2, and "" will be 0 (for missing genes)
  gene_matrix_numeric <- gene_matrix
  gene_matrix_numeric[gene_matrix == "up"] <- 1
  gene_matrix_numeric[gene_matrix == "down"] <- 2
  gene_matrix_numeric[gene_matrix == ""] <- 0

  # Convert to numeric for pheatmap
  gene_matrix_numeric <- as.numeric(gene_matrix_numeric)

  # Reshape into the same matrix format
  gene_matrix_numeric <- matrix(gene_matrix_numeric, nrow = length(gene_lists),
                                ncol = length(unique_genes),
                                dimnames = list(names(gene_lists), unique_genes))

  # Define the color palette for the heatmap
  color_palette <- c("white", "blue", "red")

  # Create the heatmap
  pheatmap(gene_matrix_numeric,
           color = color_palette,
           legend = FALSE,
           show_rownames = TRUE,
           show_colnames = TRUE,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           fontsize_col = 5, # Adjust the size of the gene names for readability
           fontsize_row = 10, # Adjust the size of the list names
           main = "Regulon Gene Activity Heatmap per TP53 mutation type")
}

