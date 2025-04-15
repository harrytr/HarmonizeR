
isRegulon <-function(genes_file, TF){
  library(dplyr)

  regulon_df <- import_TFregulons_Interactions(select_organism = 9606) #keeping human regulons
  
  regulon_df <-  regulon_df[which((regulon_df$source_genesymbol %in% TF)), ]
  
  df <- matrix(data = , nrow = nrow(regulon_df), ncol = 5) #creating the regulon dataframe for the createRegulonList function
  
  df[, 1] = regulon_df$source_genesymbol
  df[, 3] = regulon_df$target_genesymbol
  df[, 4] = regulon_df$tfregulons_level
  df[, 5] = regulon_df$sources
  
  df[which(regulon_df$is_stimulation==1), 2] <- 1
  df[which(regulon_df$is_inhibition==1), 2] <- -1
  colnames(df) = c("Source", "Sign", "Target", "Confidence", "source")
  df <- as.data.frame(df)
  df$Source = as.character(df$Source)
  df$Sign = as.numeric(as.character(df$Sign))
  df$Target = as.character(df$Target)
  df$Confidence = as.character(df$Confidence)
  df$source = as.character(df$source)

  write.csv(df,"isRegulon.csv")
  return(df)
}