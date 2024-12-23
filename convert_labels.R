
convert_labels <- function() {
  


df <- as.data.frame(read.csv("GNN_labels.csv"))



df <- df %>% mutate(labels = case_when((str_detect(filename, paste0("missense")) == TRUE & str_detect(filename, paste0("Metastatic")) == TRUE) ~ "missense_Metastatic",
                                      (str_detect(filename, paste0("WT")) == TRUE & str_detect(filename, paste0("Metastatic")) == TRUE) ~ "WT_Metastatic",
                                      (str_detect(filename, paste0("missense")) == TRUE & str_detect(filename, paste0("Primary")) == TRUE) ~ "missense_Primary",
                                      (str_detect(filename, paste0("WT")) == TRUE & str_detect(filename, paste0("Primary")) == TRUE) ~ "WT_Primary"))



# df <- df %>% mutate(labels = case_when((str_detect(filename, paste0("NOTisDeleterious")) == TRUE & str_detect(filename, paste0("Metastatic")) == TRUE) ~ "ND_Metastatic",
#                                       (str_detect(filename, paste0("_isDeleterious")) == TRUE & str_detect(filename, paste0("Metastatic")) == TRUE) ~ "D_Metastatic",
#                                       (str_detect(filename, paste0("NOTisDeleterious")) == TRUE & str_detect(filename, paste0("Primary")) == TRUE) ~ "ND_Primary",
#                                       (str_detect(filename, paste0("_isDeleterious")) == TRUE & str_detect(filename, paste0("Primary")) == TRUE) ~ "D_Primary"))


# df <- df %>% mutate(labels = case_when((str_detect(filename, paste0("_Metastatic")) == TRUE) ~ "Metastatic",
#                                        (str_detect(filename, paste0("_Primary"))== TRUE) ~ "Primary"))
                                      #(str_detect(filename, paste0("WT")) == TRUE) ~ "WT"))
 


#df <- df %>% mutate(labels = case_when((str_detect(filename, paste0("_NOTisDeleterious")) == TRUE & str_detect(filename, paste0("Metastatic")) == TRUE & str_detect(filename, paste0("missense")) == TRUE) ~ "M_ND_Metastatic",
#                                      (str_detect(filename, paste0("_isDeleterious")) == TRUE & str_detect(filename, paste0("Metastatic")) == TRUE & str_detect(filename, paste0("missense")) == TRUE) ~ "M_D_Metastatic",
#                                      (str_detect(filename, paste0("NOTisDeleterious")) == TRUE & str_detect(filename, paste0("Primary")) == TRUE & str_detect(filename, paste0("missense")) == TRUE) ~ "M_ND_Primary",
  #                                    (str_detect(filename, paste0("NOTisDeleterious")) == FALSE & str_detect(filename, paste0("Primary")) == TRUE & str_detect(filename, paste0("missense")) == TRUE) ~ "M_D_Primary"))




print(nrow(df))
#df <- df %>% dplyr::filter (labels %in% c("Missense_Primary","Missense_Metastatic","WT"))
print(nrow(df))
print(unique(df$labels))
print(table(df$labels))
write.csv(df,"GNN_labels_PM2.csv")


}