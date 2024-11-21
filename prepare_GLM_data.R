prepare_GLM_data <- function(RNAseq_matrix,
                             rgenes_e,
                             mapk_data,
                             r_expr_matrix,
                             load_other_GLM,
                             GLM_all,
                             reg_type,
                             eXML_dir,
                             key,
                             disease_name,
                             do_GLM,
                             inputs_dir,LM, folds,GLM_user_word)  {


  if (do_GLM == TRUE && load_other_GLM == FALSE)
  {

    print("    Using CCLE expression data for GLM...")
    #RNAseq_matrix[,-1] <- RNAseq_matrix[,-1] #the "-1" excludes column 1



    # choose data for GLM, RNA raw seq counts or log2 normalized?

    print("    Preparing the RNAseq matrix...")
    if (GLM_all == TRUE) {

      r_expr_matrix_RNA_seq <- RNAseq_matrix

    }

    else {
      r_expr_matrix_RNA_seq <- RNAseq_matrix %>% dplyr::select(matches(rgenes_e))
    }

    only_mutations <- mapk_data %>% dplyr::select(c("CELLLINE","Hugo_Symbol","Variant_Classification"))
    WT_violin_gene <- r_expr_matrix %>% dplyr::filter(!(CELLLINE %in% mapk_data$CELLLINE))
    WT_violin_gene <- WT_violin_gene %>% dplyr::select(c("CELLLINE"))
    mapk_data_RNA <- merge(only_mutations,WT_violin_gene, by = "CELLLINE" ,all = TRUE)
    mapk_data_RNA <- merge(mapk_data_RNA,r_expr_matrix_RNA_seq, by = "CELLLINE")
    # order by type of mutation in the output file
    mapk_data_RNA<-mapk_data_RNA[order(mapk_data_RNA[,'Variant_Classification']), ]
    mapk_data_RNA$Variant_Classification[is.na(mapk_data_RNA$Variant_Classification)] <- "WT"
    mapk_data_RNA <- mutate(mapk_data_RNA, Variant_Classification = as.character(Variant_Classification))
    mapk_data_RNA[is.na(mapk_data_RNA)] <- 0

    # remove all WT
    if (key != "WT") {
      mapk_data_RNA <- mapk_data_RNA[-(grep("WT", mapk_data_RNA$Variant_Classification)),]
    }

    print("    Preparing expression matrix for Generalized Linear Models regression...")

    saveRDS(mapk_data_RNA[,-(1:2)], file = "regression.rds")
    write.csv(mapk_data_RNA[,-(1:2)], "Renoir_CCLE.csv")

    resDir <- eXML_dir


      print("    Trying Generalized Linear Regression Modelling (Training) ...")
     glregression <- MNR("regression.rds", 50, disease_name,
                          reg_type, resDir, "alternative",key,LM, folds,GLM_user_word)




  }
  return(print("    GLM check completed."))
}
