prepare_Carnival <- function(mapk_data,
                             filtered_expression_matrix,
                             USER_EM2,
                             results_dir,
                             inputs_dir,
                             disease_filename,
                             j,
                             violin_gene,
                             cloned,
                             GAP,
                             cpu_threads,
                             network_similarity_threshold,
                             network_similarity_names,
                             top_user,
                             top_score,
                             load_other,
                             do_GLM,
                             genes_e,
                             carnival_flag,
                             GLM_all,
                             rgenes_e,
                             eXML_dir,
                             reg_type,
                             key,
                             genes_HGNC_bkp,
                             regulons_violin_gene,
                             FEM_user,
                             ccle_iterator,
                             key_opt,
                             key_opt_type,mutations,source_code_dir,epochs,lr,tts,update_it,bsu, hidden_dim, dropout, GNN_flag)

{
  print("Preparing Carnival data...")

  radar_plot_data <<- as.data.frame(matrix(data = , nrow = length(disease_filename), ncol = (length(network_similarity_threshold))*(length(network_similarity_names))))

  rownames(radar_plot_data)[j] <- disease_filename[j]
  f= 1
    for (j in 1:length(network_similarity_names)){
      for (k in 1:length(network_similarity_threshold)){
        colnames(radar_plot_data)[f] <- paste0(network_similarity_names[j],"_",network_similarity_threshold[k])
        f= f+ 1;
      }
    }

    `%!in%` = Negate(`%in%`)
    mapk_data_carnival <- merge(mapk_data,filtered_expression_matrix, by = "CELLLINE", all = TRUE)
    NA_cellines <- mapk_data_carnival %>% dplyr::filter(is.na(mapk_data_carnival$Variant_Classification))
    mapk_data_carnival_only_expression <- mapk_data_carnival[, -c(2)]
    USER_EM2 <- USER_EM2 %>% dplyr::filter(CELLLINE %in% NA_cellines$CELLLINE)
    mapk_data_carnival_small <- mapk_data_carnival %>% dplyr::select(c("CELLLINE", "Variant_Classification"))
    df_merged <- merge(mapk_data_carnival_small, USER_EM2, by = "CELLLINE", all = TRUE)
    
    df_merged$Variant_Classification <- ifelse(is.na(df_merged$Variant_Classification.x), df_merged$Variant_Classification.y, df_merged$Variant_Classification.x)
    
    df_merged <- df_merged[c("CELLLINE", "Variant_Classification")]
    write.csv(df_merged,"celllines_variant.csv")
    
    mapk_data_carnival <- NULL
    mapk_data_carnival <- merge(df_merged,filtered_expression_matrix, by = "CELLLINE", all = TRUE)
    mapk_data_carnival <- as.data.frame(mapk_data_carnival)
    #write.csv(mapk_data_carnival,"test.csv")
    setwd(results_dir)
    
    mapk_data_carnival$CELLLINE <-NULL
    colnames(mapk_data_carnival) <- sapply(strsplit(colnames(mapk_data_carnival), split='..', fixed=TRUE),function(x) (x[1]))
    carnival <- t(mapk_data_carnival)
    carnival <- as.data.frame(carnival)
    colnames(carnival) <- mapk_data_carnival[,1]
    carnival <- carnival[-1,]
    carnival[,1] <- rownames(carnival)
    rownames(carnival) <- NULL
    colnames(carnival)[1] <- "Gene"
    carnival <- carnival[,colSums(is.na(carnival))<nrow(carnival)]
    write.csv(carnival, "Carnival_EM.csv")
    df_EM<-as.data.frame(read.csv("Carnival_EM.csv", row.names = 'Gene', header = TRUE))

    # if (key_opt != "WT") {
    #   print("Removing WT samples...")
    # #df_EM <- df_EM[,-(grep("WT", colnames(df_EM)))]
    # df_EM <- select(df_EM, -contains("WT"))
    # }

    #write.csv(mapk_data_carnival,"Carnival_EM.csv")
    df_EM <- df_EM[,-1]
    write.csv(df_EM,"Carnival_input.csv")

    if (.Platform$OS.type == "unix") {
      
      Carnival_opt_res <- Carnival_opt(ccle_iterator,
                                           df_EM,
                                           results_dir,
                                           inputs_dir,
                                           disease_filename[j],
                                           violin_gene,
                                           cloned,
                                           GAP,
                                           cpu_threads,
                                           network_similarity_threshold,
                                           network_similarity_names,
                                           top_user,
                                           top_score,
                                           genes_HGNC_bkp,
                                           regulons_violin_gene,
                                           radar_plot_data,key_opt,key_opt_type, mutations,source_code_dir,epochs,lr,tts,update_it,bsu,hidden_dim, dropout, GNN_flag)
    }
    else if  (.Platform$OS.type == "windows") {
      Carnival_opt_res <- Carnival_opt_SPN(ccle_iterator,
                                           df_EM,
                                           results_dir,
                                           inputs_dir,
                                           disease_filename,
                                           violin_gene,
                                           cloned,
                                           GAP,
                                           cpu_threads,
                                           network_similarity_threshold,
                                           network_similarity_names,
                                           top_user,
                                           top_score,
                                           genes_HGNC_bkp,
                                           regulons_violin_gene,
                                           radar_plot_data,key_opt,key_opt_type, mutations,source_code_dir,epochs,lr,tts,update_it,bsu,hidden_dim, dropout, GNN_flag)
    }
    
      

    # CASE TCGA :




} # end of function
