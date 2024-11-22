prepare_Carnival <- function(mapk_data,
                             filtered_expression_matrix,
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
                             key_opt_type,mutations,source_code_dir,epochs,lr,tts,update_it)

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

  #network_similarity_threshold = c(0.5,0.75,0.9)
  #network_similarity_names <- c("All_networks","Mutation_Type","Deleterious","Hotspot")

    mapk_data_carnival <- merge(mapk_data,filtered_expression_matrix, by = "CELLLINE", all = TRUE)
    mapk_data_carnival$Variant_Classification[is.na(mapk_data_carnival$Variant_Classification)] <- "WT"
    setwd(results_dir)
    colnames(mapk_data_carnival) <- as.character(colnames(mapk_data_carnival))
    colnames(mapk_data_carnival) <- sapply(strsplit(colnames(mapk_data_carnival), split='..', fixed=TRUE),function(x) (x[1]))
    mapk_data_carnival <- mapk_data_carnival[,-(1:2),drop=F]
    mapk_data_carnival <- t(mapk_data_carnival)
    colnames(mapk_data_carnival) <- mapk_data_carnival[1,]
    mapk_data_carnival <- cbind(Gene = rownames(mapk_data_carnival), mapk_data_carnival)
    rownames(mapk_data_carnival) <- NULL
    mapk_data_carnival <- mapk_data_carnival[-1,,drop=F]
    mapk_data_carnival <- mapk_data_carnival[,colSums(is.na(mapk_data_carnival))<nrow(mapk_data_carnival)]




    write.csv(mapk_data_carnival,"Carnival_EM.csv")

    df_EM<-as.data.frame(read.csv("Carnival_EM.csv", row.names = 'Gene'), header = TRUE)


    if (key_opt != "WT") {
      print("Removing WT samples...")
    #df_EM <- df_EM[,-(grep("WT", colnames(df_EM)))]
    df_EM <- select(df_EM, -contains("WT"))
    }

    #write.csv(mapk_data_carnival,"Carnival_EM.csv")

    df_EM <- df_EM[,-1]
    write.csv(df_EM,"Carnival_input.csv")


      Carnival_opt_res <- Carnival_opt(ccle_iterator,
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
                                       radar_plot_data,key_opt,key_opt_type, mutations,source_code_dir,epochs,lr,tts,update_it)

    # CASE TCGA :




} # end of function
