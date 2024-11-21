
isRegulon <-function(genes_file, TF){
  library(dplyr)
  #load(file = system.file("BEST_viperRegulon.rdata",package="CARNIVAL"))
  #gene_list_df <- as.data.frame(read.csv(genes_file), header=TRUE)
  
  # now per cancer type:
  #gene_list_df <- gene_list_df %>%  dplyr:: filter(all %in% c("up-reg", "down-reg","small change", "Others", "not expressed","not significant"))
 # gene_list_df <- gene_list_df %>%  dplyr:: filter(all %in% c("up-reg", "down-reg"))
 # gene_list <-gene_list_df$genes
  
  #regulon_df<- import_Omnipath_Interactions(select_organism = 9606)
  
  regulon_df <- import_TFregulons_Interactions(select_organism = 9606) #keeping human regulons
  #regulon_df <- regulon_df[which(regulon_df$tfregulons_level%in%c("A","B","C","D","E")), ] #keeping 'A', 'B' and 'C' levels
  #regulon_df <- regulon_df[which(regulon_df$is_directed==1), ] #keeping only directed
  #regulon_df <- regulon_df[which((regulon_df$is_stimulation+regulon_df$is_inhibition)==1), ] #keeping only regulons which are either activations/inhibitions
  #only for the TF now:
  regulon_df <-  regulon_df[which((regulon_df$source_genesymbol %in% TF)), ]
  # find only the ones that are targets of the TF from gene_list:
  #regulon_df <-  regulon_df[which((regulon_df$target_genesymbol %in% gene_list)), ]
  df <- matrix(data = , nrow = nrow(regulon_df), ncol = 5) #creating the regulon dataframe for the createRegulonList function
  #print(df)
  
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
  #regulons = createRegulonList(regulon_table = df)
  
  write.csv(df,"isRegulon.csv")
  return(df)
}