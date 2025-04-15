update_labels <- function() {
  # Define the mutation types to look for
  mutation_types <- c(
    "Missense_Mutation", "Nonsense_Mutation", "In_Frame_Del", 
    "Splice_Site", "Frame_Shift_Ins", "Frame_Shift_Del", "In_Frame_Ins", "WT", "Fusion", "Splice_Region"
  )
  
  # Read the CSV file
  df <- read.csv("GNN_labels_TCGA_MT.csv", stringsAsFactors = FALSE)
  
  # Function to extract mutation type from filename
  extract_label <- function(filename) {
    for (mutation in mutation_types) {
      if (grepl(mutation, filename)) {
        return(mutation)
      }
    }
    return(NA)
  }
  
  # Apply the function to each row
  df$labels <- sapply(df$filename, extract_label)
  
  # Optionally write back to file
  df <- na.omit(df)
  write.csv(df, "GNN_labels_TCGA_updated.csv", row.names = FALSE)
  
  return(df)
}