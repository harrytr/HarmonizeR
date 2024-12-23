visualize_gene_sigs <- function(gene_lists) {


  library(ggplot2)
  library(dplyr)

  # Convert input list of lists into a data frame
  gene_data <- do.call(rbind, lapply(seq_along(gene_lists), function(i) {
    list_data <- gene_lists[[i]]
    data.frame(
      Gene = gsub("_up|_down", "", list_data),
      Status = ifelse(grepl("_up", list_data), "up", "down"),
      List = paste("List", i),
      stringsAsFactors = FALSE
    )
  }))

  # Assign colors based on "up" or "down" status
  gene_data <- gene_data %>%
    mutate(Color = ifelse(Status == "up", "blue", "red"))

  # Plot the gene data
  p <- ggplot(gene_data, aes(x = Gene, y = List, color = Color)) +
    geom_text(aes(label = Gene), size = 4, show.legend = FALSE) +
    scale_color_manual(values = c("red", "blue")) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    labs(title = "Gene Visualization by Status")

  # Display the plot
  print(p)
}
