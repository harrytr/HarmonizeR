# Load necessary libraries
library(dplyr)

# Read the mutation types from TCGA_mutations.csv
mutation_types <- read.csv("TCGA_mutations.csv", stringsAsFactors = FALSE)$X1

# Read the GNN_labels.csv
GNN_labels <- read.csv("GNN_labels.csv", stringsAsFactors = FALSE)

# Function to find mutation type in each row
detect_mutation <- function(text, mutations) {
  match <- mutations[grepl(paste(mutations, collapse = "|"), text)]
  if (length(match) > 0) {
    return(match[1])  # Take the first match found
  } else {
    return(NA)
  }
}

# Apply the function to the first column of GNN_labels
GNN_labels$Detected_Mutation <- sapply(GNN_labels[[1]], detect_mutation, mutations = mutation_types)

# Save the updated dataframe
write.csv(GNN_labels, "GNN_labels_updated.csv", row.names = FALSE)
