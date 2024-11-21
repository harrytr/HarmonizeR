
createRegulons <-function(results_dir){

load(file = system.file("BEST_viperRegulon.rdata",package="CARNIVAL"))

  
#regulon_df<- import_Omnipath_Interactions(select_organism = 9606)

regulon_df <- import_TFregulons_Interactions(select_organism = 9606) #keeping human regulons
#regulon_df <- regulon_df[which(regulon_df$tfregulons_level%in%c("A", "B", "C")), ] #keeping 'A', 'B' and 'C' levels
regulon_df <- regulon_df[which(regulon_df$tfregulons_level%in%c("A","B","C","D","E")), ] #keeping 'A', 'B' and 'C' levels

regulon_df <- regulon_df[which(regulon_df$is_directed==1), ] #keeping only directed
regulon_df <- regulon_df[which((regulon_df$is_stimulation+regulon_df$is_inhibition)==1), ] #keeping only regulons which are either activations/inhibitions

df <- matrix(data = , nrow = nrow(regulon_df), ncol = 3) #creating the regulon dataframe for the createRegulonList function
#print(df)

df[, 1] = regulon_df$source_genesymbol
df[, 3] = regulon_df$target_genesymbol
df[which(regulon_df$is_stimulation==1), 2] <- 1
df[which(regulon_df$is_inhibition==1), 2] <- -1
colnames(df) = c("Source", "Sign", "Target")
df <- as.data.frame(df)
df$Source = as.character(df$Source)
df$Sign = as.numeric(as.character(df$Sign))
df$Target = as.character(df$Target)

regulons = createRegulonList(regulon_table = df)

save(regulons, file = paste0(results_dir,"/regulons.RData"))

return(df)
}