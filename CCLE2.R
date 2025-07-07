CCLE2 <-function(disease_name,
                 PMA_user,
                 Age_user,
                 Sex_user,
                 OncotreeLineage_user,
                 target_gene,
                 do_GLM,
                 reg_type,
                 carnival_flag,
                 GLM_all,
                 new_version,
                 dataset_version,
                 GLM_user_word,
                 key,
                 load_other_GLM,
                 FEM_user,
                 LM,
                 folds,
                 key_opt,
                 key_opt_type,
                 molecular_features,
                 clinical_features,
                 epochs,lr,tts,update_it,bsu,hidden_dim,dropout,GNN_flag)



{
  options(warn=-1)
  GAP <- 0.05
  library(parallel)
  cpu_threads <- detectCores()
  top_user <- 0
  top_score <- 0.8
  source_code_dir <- getwd()

  if (  (molecular_features=="") & (clinical_features==""))
  {
    print("Warning: No selected Clinical or Molecular feature found. Defaulting to WT...")
    key_opt <- "WT"
    key_opt_type <- "Equal"
  }

  #options(warn=-1)
  wd = getwd()
  graphics.off()
 
  source(paste0(wd,'/analysis.R'))
  source(paste0(wd,'/runDoRothEA.R'))
  source(paste0(wd,'/createRegulonList.R'))
  source(paste0(wd,'/generateDataframeTF.R'))
  source(paste0(wd,'/Carnival_opt_SPN.R'))
  source(paste0(wd,'/Carnival_opt.R'))
  source(paste0(wd,'/createRegulons.R'))
  source(paste0(wd,'/createRegulonList.R'))
  source(paste0(wd,'/prepare_GLM_data.R'))
  source(paste0(wd,'/prepare_Carnival.R'))
  source(paste0(wd,"/assignPROGENyScores.R"))


  hotspots_bkp <- c("p.R175H","p.R248Q","p.R273H","p.R248W", "p.R273C", "p.R282W", "p.G245S")
  inputs_dir = paste0(wd,"//inputs")
  models <- as.data.frame(read.csv(paste0(inputs_dir,"/Model.csv"), header=TRUE))
  CCLE_models <- unique(models$DepmapModelType)
  results_dir = paste0(wd,"//CCLE_",target_gene)

  print(results_dir)
  dir.create(path = results_dir)
  results_dir_bkp <- results_dir

  if (new_version==TRUE) {

    new_version_folder <- as.character(dataset_version)
    print(paste0("Preparing input dataset for new CCLE version ", as.character(dataset_version[1])))
    print(paste0("Reading copy number matrix for ",as.character(dataset_version[1])))
    cn_csv <- as.data.frame(read.csv(paste0(inputs_dir,"/",new_version_folder,"/OmicsCNGene.csv"), header=TRUE))
    print(paste0("Reading expression matrix for ",as.character(dataset_version[1])))
    expr_matrix_csv <- as.data.frame(read.csv(paste0(inputs_dir,"/",new_version_folder,"/OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected.csv"), header=TRUE))
    print(paste0("Reading mutation profiles matrix for ",as.character(dataset_version[1])))
    mut_matrix_csv <- as.data.frame(read.csv(paste0(inputs_dir,"/",new_version_folder,"/OmicsSomaticMutations.csv"), header=TRUE))


    save(expr_matrix_csv,
         mut_matrix_csv,
         cn_csv,
         file=paste0(inputs_dir,"/",
                     as.character(dataset_version),".RData"))

  }


  print("Checking column names for newer version of CCLE mutation profiles...")
  try(colnames(mut_matrix_csv)[which(colnames(mut_matrix_csv)=="ModelID")] <- "Tumor_Sample_Barcode"
  )

  try(colnames(mut_matrix_csv)[which(colnames(mut_matrix_csv)=="VariantInfo")] <- "Variant_Classification"
  )


  mut_matrix_csv$Variant_Classification <- str_sub(mut_matrix_csv$Variant_Classification, start = 1, end=20)

  try(colnames(mut_matrix_csv)[which(colnames(mut_matrix_csv)=="HugoSymbol")] <- "Hugo_Symbol"
  )

  try(colnames(mut_matrix_csv)[which(colnames(mut_matrix_csv)=="ProteinChange")] <- "Protein_Change"
  )

  try(colnames(mut_matrix_csv)[which(colnames(mut_matrix_csv)=="Sift")] <- "isDeleterious"
  )

  try(colnames(mut_matrix_csv)[which(colnames(mut_matrix_csv)=="DNAChange")] <- "Codon_Change"
  )

  try(colnames(mut_matrix_csv)[which(colnames(mut_matrix_csv)=="Hotspot")] <- "isHotspot"
  )


  mut_matrix_csv$isDeleterious <- ifelse(str_detect(mut_matrix_csv$isDeleterious,"deleterious"), TRUE, FALSE)

  print("Reading the cell line names mapping and suppliers list...")
  names_dir = paste0(inputs_dir, "//", "Model.csv")
  cell_line_names <- as.data.frame(read.csv(names_dir, header = TRUE))

  try(colnames(cell_line_names)[which(colnames(cell_line_names)=="ModelID")] <- "BROAD_ID"
  )

  try(colnames(cell_line_names)[which(colnames(cell_line_names)=="CCLEName")] <- "CCLE_ID"
  )


  names_mat <-  cell_line_names  %>% dplyr::select(c("BROAD_ID","CCLE_ID"))
  colnames(names_mat) <- c("CELLLINE","CCLE_ID")



  plot_list <- c()

  lib_dir <- paste0(getwd(),"/libs")
  .libPaths( c( lib_dir , .libPaths() ) )

  print("Loading libraries required...")
  list.of.packages <- c("ggstatsplot","dplyr","ggplot2","ggrepel","ggpubr","viridis","tibble","stringr",
                        "corrplot","tidyverse","igraph","visNetwork", "data.table", "CARNIVAL",
                        "viper", "CellNOptR", "OmnipathR", "stringi","openxlsx",
                        "sna", "gplots","ggfortify","limma", "UpSetR","survival", "survminer","ggcorrplot","rstatix")# "edgeR",

  invisible(lapply(list.of.packages, library, character.only = TRUE))

  # user specifies the name of the cancer type to analyze (can input also "all" to analyze all of the supplied cell line files)
  disease_filename <- c()

  # measuring run-time
  start <- Sys.time()

  # TP53, CDC20, CENPA, KIF2C, PLK1 etc.
  violin_gene_name <- target_gene


  violin_gene <- strsplit(violin_gene_name,split='..', fixed=TRUE)
  violin_gene <- violin_gene[[1]][1]
  # rename some columns to tag cell lines in expression, mutation and CN matrices :

  colnames(expr_matrix_csv)[1] <- "CELLLINE"
  colnames(cn_csv)[1] <- "cell_lines"


  violin_gene_E_f <- paste0("^(", paste0("CELLLINE|",violin_gene_name), ")")
  violin_gene_E <- expr_matrix_csv %>% dplyr::select(matches(violin_gene_E_f))
  cell_line_dir = paste0(inputs_dir, "//", "CELL_LINES")


  if (disease_name == "all_custom") {
    disease_filename <- list.files(cell_line_dir)
  }
  else if (disease_name == "all_CCLE") {
    disease_filename <- c(disease_filename,CCLE_models)
  }
  else
  {
    disease_filename <- c(disease_filename,disease_name)

  }

  # csv file to save the results from CARNIVAL for the similarity comparison of the networks
  network_similarity_threshold = c(0.5,0.75,0.9)
  network_similarity_names <- c("All_networks","Mutation_Type","Deleterious","Hotspot",paste0(key_opt,"_",key_opt_type))

  print("Creating pdf results file...")
  setwd(results_dir)
  # create graph visual dir
  dynamic_graphs_dir = paste0(results_dir,"/dynamic_graphs")
  dir.create(path = dynamic_graphs_dir)

  eXML_dir <- paste0(results_dir,"/GLM")
  if (do_GLM == TRUE) {dir.create(path = eXML_dir)}


  all_cancers <- NULL
  all_cancers_all_genes <- NULL


  print(paste0("Reading the regulons of ", violin_gene," from DoRothEA (all confidence levels)"))
  regulons_dir = paste0(inputs_dir, "/", "regulons.csv")
  regulons_csv <- as.data.frame(read.csv(regulons_dir, header = TRUE))

  regulons <- regulons_csv  %>%  dplyr:: filter(source_genesymbol %in% violin_gene)
  regulons_violin_gene <- NULL

  regulons_violin_gene <- unique(as.character(regulons$target_genesymbol))
  #write.csv(regulons_violin_gene,paste0("regulon_",violin_gene,"_omnipath_downloaded.csv"))
  setwd(results_dir)
  # create violin expression data folder
  violin_expr_path = paste0(results_dir,"/", violin_gene, "_expression_data")
  violin_mut_path = paste0(results_dir,"/", violin_gene, "_mutation_data")
  cell_line_IDs_names = paste0(results_dir,"/", "_cell_line_IDs_names_data")
  pca_objects = paste0(results_dir,"/", "pca_data")
  dir.create(path = violin_expr_path)
  dir.create(path = violin_mut_path)
  dir.create(path = cell_line_IDs_names)
  dir.create(path = pca_objects)
  write.csv(regulons_violin_gene,paste0(violin_gene,"_regulons.csv"))
  setwd(inputs_dir)


  model_temp <- models

  if (PMA_user == "ALL") {

    #disease_cell_lines <- model_temp$ModelID
    if (Sex_user != "ALL") {
      model_temp <- model_temp %>% dplyr::filter(Sex %in% Sex_user)
    }
    if (Age_user != "ALL") {
      model_temp <- model_temp %>% dplyr::filter(AgeCategory %in% Age_user)
    }

  }
  else {

    model_temp <- model_temp %>% dplyr::filter(PrimaryOrMetastasis %in% PMA_user)
    #disease_cell_lines <- model_temp$ModelID
    if (Sex_user != "ALL") {
      model_temp <- model_temp %>% dplyr::filter(Sex %in% Sex_user)
    }
    if (Age_user != "ALL") {
      model_temp <- model_temp %>% dplyr::filter(AgeCategory %in% Age_user)
    }
  }

  genes_HGNC <- unique(c(regulons_violin_gene,violin_gene))
  library(stringi)
  disease_filename <- stri_remove_empty(disease_filename)
 
  ################################################### MAIN LOOP ##################################################
  for (j in 1: length(disease_filename)) {
    print(disease_filename)
    Sys.sleep(3)
    print("Main loop...")
    ccle_iterator <- j

    temp_dir <- paste0(results_dir_bkp ,"/cancers_",disease_filename[j])
    dir.create(path = temp_dir)
    results_dir <- temp_dir
    setwd(temp_dir)


    pdffile = paste0("Plots_cancers_",disease_filename[j],"_",dataset_version,".pdf")
    pdf(file = pdffile, height = 10, width = 18, onefile = TRUE)
    ###### separating title page for each cancer type #####
    a = paste0("Analyses for ", disease_filename[j])


    plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
         xaxt='n', yaxt='n', xlab='', ylab='')
    text(1,4,a, pos=4)

    if (disease_name == "all_cell_lines") {

      disease_cell_lines = model_temp$ModelID

    }
    else {


      if (disease_filename[j] %in% CCLE_models) {
        print(paste0("Reading single cancer type from predefined CCLE models...",disease_filename[j]))
        cell_lines <- model_temp %>% dplyr::filter(DepmapModelType %in% disease_filename[j])
        disease_cell_lines <- cell_lines %>% dplyr::select("ModelID")
        disease_cell_lines <- as.data.frame(disease_cell_lines)
        disease_cell_lines <- disease_cell_lines[,1]
      }
      else if (disease_filename[j] == "ALL_MODELS"){
        disease_filename[j] <- OncotreeLineage_user
        print(paste0("Reading single cancer type from predefined CCLE models (whole lineage)...",disease_filename[j]))
        cell_lines <- model_temp %>% dplyr::filter(OncotreeLineage %in% OncotreeLineage_user)
        disease_cell_lines <- cell_lines %>% dplyr::select("ModelID")
        disease_cell_lines <- as.data.frame(disease_cell_lines)
        disease_cell_lines <- disease_cell_lines[,1]
      }
      else{
        print(disease_filename[j])
        print(paste0("Reading custom single cancer type...",disease_filename[j]))
        command <- paste0(cell_line_dir,"/",disease_filename[j],".csv")
        command <- shQuote(command)
        command <- paste0("read.csv(",command)
        command <- paste0(command,",header = TRUE)")
        command <-paste0("as.data.frame(",command)
        command <- paste0(command,")")
        disease_csv <- eval(parse(text = command))
        disease_cell_lines <- disease_csv[, 1]
      }

    }
    print(disease_cell_lines)
    violin_column <- NULL

    ################################Processing names of genes#################
    exact_genes <- NULL
    exact_genes <- paste0("^",genes_HGNC, "$", collapse="|")

    genes <- paste0(genes_HGNC,"\\s*?\\.{2}",collapse="|")
    genes0 <- paste0("\\b",genes_HGNC,"\\s*?\\.{2}\\b", collapse="|")
    genes_e <- paste0("^(", paste0("CELLLINE|",genes), ")")
    genes_cn <- paste0("^(", paste0("cell_lines|",genes0), ")")
    tumor_samples <- paste0("^(", paste(disease_cell_lines, collapse="|"), ")")

    print(paste0("Filtering the expression profiles in ",disease_filename[j]))

    ############################################################################################################

    filtered_expression_matrix <- expr_matrix_csv  %>%  dplyr::filter(CELLLINE %in% disease_cell_lines)
    filtered_expression_matrix_disease <- filtered_expression_matrix
    # create a csv with the cell line names and BROAD IDs only for the specific disease
    names_mat <-  cell_line_names  %>% dplyr::select(c("BROAD_ID","CCLE_ID"))
    colnames(names_mat) <- c("CELLLINE","CCLE_ID")

    small_filtered_expression_matrix <- filtered_expression_matrix %>% dplyr::select(c("CELLLINE"))
    cell_line_map <- merge(names_mat,small_filtered_expression_matrix, by = "CELLLINE")
    write.csv(cell_line_map,paste0(cell_line_IDs_names,"/","_cell_line_mapping_",disease_filename[j]))
    colnames(names_mat) <- c("cell_lines","CCLE_ID")

    if (ncol(violin_gene_E)!=2) {
      print(paste0("No expression data was found for ", violin_gene, ". Exiting..."))
      stop()
    }
    else{
      setwd(results_dir)
      print("Processing expression data now...")
      # selected expression data only for the violin gene:
      violin_expression_data <- NULL

      if ((clinical_features != "") &  !(clinical_features %in% c("OncotreeLineage","AgeCategory","Sex","PrimaryOrMetastasis")) == TRUE ) {

        names_mat <-  cell_line_names  %>% dplyr::select(c("BROAD_ID","CCLE_ID","OncotreeLineage","AgeCategory","Sex","PrimaryOrMetastasis",clinical_features))
        colnames(names_mat) <- c("CELLLINE","CCLE_ID","OncotreeLineage","AgeCategory","Sex","PrimaryOrMetastasis", clinical_features)
      }
      else{
        names_mat <-  cell_line_names  %>% dplyr::select(c("BROAD_ID","CCLE_ID","OncotreeLineage","AgeCategory","Sex","PrimaryOrMetastasis"))
        colnames(names_mat) <- c("CELLLINE","CCLE_ID","OncotreeLineage","AgeCategory","Sex","PrimaryOrMetastasis")
      }

      violin_expression_data <- filtered_expression_matrix %>%  dplyr::select(matches(violin_gene_E_f))
      violin_expression_data <- merge(names_mat,violin_expression_data, by = "CELLLINE")
      write.csv(violin_expression_data,paste0(violin_expr_path,"/",violin_gene,"_expression_in_",disease_filename[j],".csv"))
    }

    filtered_expression_matrix <- filtered_expression_matrix %>%  dplyr::select(matches(genes_e))

    if (ncol(filtered_expression_matrix) <=2) {
      print(paste0("No expression data was found for your signature. Exiting..."))
      stop()
    }

    setwd(results_dir)

    ############################################################################################################
    print(paste0("Filtering the mutation profiles across only the GRN and Regulons in ",disease_filename[j]))
    filtered_mutation_matrix <- mut_matrix_csv %>% dplyr::filter(Tumor_Sample_Barcode %in% disease_cell_lines)

    skip_C = FALSE
    print(disease_cell_lines)
    # we now filter only for the genes that also belong to our GRN
    filtered_mutation_matrix <- filtered_mutation_matrix %>% dplyr::filter(stringr::str_detect(Hugo_Symbol, exact_genes))
    if (!(violin_gene %in% filtered_mutation_matrix$Hugo_Symbol) | length(disease_cell_lines) < 3) {
      print(paste0("No enough mutation data found for ", violin_gene, "!"))
      skip_C = TRUE
    }

    if (skip_C == FALSE) {
      setwd(results_dir)

      # selected expression data only for the violin gene:
      violin_mutation_data <- NULL

      if ((clinical_features != "") &  !(clinical_features %in% c("OncotreeLineage","AgeCategory","Sex","PrimaryOrMetastasis"))== TRUE) {
        names_mat <-  cell_line_names  %>% dplyr::select(c("BROAD_ID","CCLE_ID","OncotreeLineage","AgeCategory","Sex","PrimaryOrMetastasis",clinical_features))
        colnames(names_mat) <- c("Tumor_Sample_Barcode","CCLE_ID","OncotreeLineage","AgeCategory","Sex","PrimaryOrMetastasis", clinical_features)
      }

      else
      {
        names_mat <-  cell_line_names  %>% dplyr::select(c("BROAD_ID","CCLE_ID","OncotreeLineage","AgeCategory","Sex","PrimaryOrMetastasis"))
        colnames(names_mat) <- c("Tumor_Sample_Barcode","CCLE_ID","OncotreeLineage","AgeCategory","Sex","PrimaryOrMetastasis")
      }




      violin_mutation_data <- filtered_mutation_matrix %>%  dplyr::filter(Hugo_Symbol %in% violin_gene)
      violin_mutation_data <- merge(names_mat,violin_mutation_data, by = "Tumor_Sample_Barcode")
      write.csv(violin_mutation_data,paste0(violin_mut_path,"/",violin_gene,"_mutations_in_",disease_filename[j],".csv"))


      rgenes <- paste0(regulons_violin_gene,"\\s*?\\.{2}",collapse="|")


      rgenes_e <- paste0("^(", paste0("CELLLINE|",rgenes), ")")

      # #  regression data-set :

      mapk_data <- mut_matrix_csv %>% dplyr::filter(Tumor_Sample_Barcode %in% disease_cell_lines)
      mapk_data <- mapk_data %>% dplyr::filter(Hugo_Symbol %in% violin_gene)
      mapk_data_basis <- mapk_data

      # get whether a mutation is deleterious for the gene function or not so as to approp. change perturbation for CARNIVAL to -1 or 1
      # isDeleterious <-  mapk_data %>% dplyr::select("isDeleterious")

      if ((molecular_features != "") & !(molecular_features %in% c("Tumor_Sample_Barcode", "Hugo_Symbol",
                                                                   "Variant_Classification","Protein_Change","isDeleterious"))== TRUE) {
        mapk_data <- mapk_data %>% dplyr::select("Tumor_Sample_Barcode", "Hugo_Symbol",
                                                 "Variant_Classification","Protein_Change","isDeleterious",molecular_features)
        mapk_data <- mapk_data %>% group_by(Variant_Classification) %>% mutate(Unique_Mut_ID = cur_group_id())
        colnames(mapk_data) <- c("CELLLINE", "Hugo_Symbol", "Variant_Classification","Protein_Change","isDeleterious", "Mutation_ID",molecular_features)
      }
      else
      {
        mapk_data <- mapk_data %>% dplyr::select("Tumor_Sample_Barcode", "Hugo_Symbol",
                                                 "Variant_Classification","Protein_Change","isDeleterious")
        mapk_data <- mapk_data %>% group_by(Variant_Classification) %>% mutate(Unique_Mut_ID = cur_group_id())
        colnames(mapk_data) <- c("CELLLINE", "Hugo_Symbol", "Variant_Classification","Protein_Change","isDeleterious", "Mutation_ID")
      }


      if ((clinical_features != "") &  !(clinical_features %in% c("OncotreeLineage","AgeCategory","Sex","PrimaryOrMetastasis"))== TRUE) {

        names_mat <-  cell_line_names  %>% dplyr::select(c("BROAD_ID","CCLE_ID","OncotreeLineage","AgeCategory","Sex","PrimaryOrMetastasis", clinical_features))
        colnames(names_mat) <- c("CELLLINE","CCLE_ID","OncotreeLineage","AgeCategory","Sex","PrimaryOrMetastasis", clinical_features)
      }
      else{
        names_mat <-  cell_line_names  %>% dplyr::select(c("BROAD_ID","CCLE_ID","OncotreeLineage","AgeCategory","Sex","PrimaryOrMetastasis"))
        colnames(names_mat) <- c("CELLLINE","CCLE_ID","OncotreeLineage","AgeCategory","Sex","PrimaryOrMetastasis")
      }

      mapk_data <- merge(mapk_data,names_mat, by = "CELLLINE")

      #####################################################################################################################################################
      # insert start_end positions now in a copy :

      if ((molecular_features != "") & !(molecular_features %in% c("Tumor_Sample_Barcode", "Hugo_Symbol",
                                                                   "Variant_Classification","Protein_Change","isDeleterious"))== TRUE) {

        mapk_data_SE <- mapk_data_basis %>% dplyr::select("Tumor_Sample_Barcode", "Hugo_Symbol",
                                                          "Variant_Classification","Protein_Change","isDeleterious","isHotspot",molecular_features)
        mapk_data_SE <- mapk_data_SE %>% group_by(Variant_Classification) %>% mutate(Unique_Mut_ID = cur_group_id())
        colnames(mapk_data_SE) <- c("CELLLINE", "Hugo_Symbol", "Variant_Classification","Protein_Change","isDeleterious","isHotspot","Mutation_ID",molecular_features)

      }
      else{
        mapk_data_SE <- mapk_data_basis %>% dplyr::select("Tumor_Sample_Barcode", "Hugo_Symbol",
                                                          "Variant_Classification","Protein_Change","isDeleterious","isHotspot")
        mapk_data_SE <- mapk_data_SE %>% group_by(Variant_Classification) %>% mutate(Unique_Mut_ID = cur_group_id())
        colnames(mapk_data_SE) <- c("CELLLINE", "Hugo_Symbol", "Variant_Classification","Protein_Change","isDeleterious","isHotspot","Mutation_ID")
      }

      if ((clinical_features != "") &  !(clinical_features %in% c("OncotreeLineage","AgeCategory","Sex","PrimaryOrMetastasis"))== TRUE) {

        names_mat <-  cell_line_names  %>% dplyr::select(c("BROAD_ID","CCLE_ID","OncotreeLineage","AgeCategory","Sex","PrimaryOrMetastasis",clinical_features))
        colnames(names_mat) <- c("CELLLINE","CCLE_ID","OncotreeLineage","AgeCategory","Sex","PrimaryOrMetastasis",clinical_features)

      }
      else
      {
        names_mat <-  cell_line_names  %>% dplyr::select(c("BROAD_ID","CCLE_ID","OncotreeLineage","AgeCategory","Sex","PrimaryOrMetastasis"))
        colnames(names_mat) <- c("CELLLINE","CCLE_ID","OncotreeLineage","AgeCategory","Sex","PrimaryOrMetastasis")
      }

      #mapk_data_SE <- merge(mapk_data_SE,names_mat, by = "CELLLINE")
      mapk_data_SE <- merge(mapk_data_SE,names_mat, by = "CELLLINE")
      non_unite_data <- mapk_data_SE

      #####################################################################################################################################################

      cloned <- merge(non_unite_data,filtered_expression_matrix, by = "CELLLINE")


      k1 <- ncol(non_unite_data)+1
      k2 <-  ncol(cloned)

      #install_github("vqv/ggbiplot")
      cloned_t <- as.data.frame(cloned)

      pc_matrix <- cloned_t[,k1:k2]

      print("Calculating PCA...")
      pc <- prcomp(pc_matrix,
                   center = TRUE,
                   scale. = TRUE)

      PCA_plot <- ggbiplot::ggbiplot(pc,
                    obs.scale = 1,
                    var.scale = 1,
                    groups = cloned$Variant_Classification,
                    ellipse = TRUE,
                    circle = TRUE,
                    ellipse.prob = 0.68)
      PCA_plot <- PCA_plot + scale_color_discrete(name = '')
      PCA_plot <- PCA_plot + theme(legend.direction = 'horizontal',
                     legend.position = 'top')
      print(PCA_plot)




      write.csv(cloned,"cloned.csv")
      print("Uniting molecular features..")
      data_heatmap <-mapk_data



      # Transform all columns so there is no confusion from multiple TRUEs etc
      valid_values <- c("TRUE","FALSE")
      for (i in 1:ncol(mapk_data)) {
        if (all(mapk_data[,i] %in% valid_values)) {
          print(paste0("Column ", colnames(mapk_data)[i], " found to be only logical"))
          print(paste0("Converting to {",colnames(mapk_data)[i],",", paste0("NOT",colnames(mapk_data)[i]),"}"))
          #invisible(readline(prompt="Press [enter] to continue"))
          Sys.sleep(3)
          mapk_data[,i] <- ifelse(mapk_data[,i] == "TRUE",colnames(mapk_data)[i],paste0("NOT",colnames(mapk_data)[i]))
        }

      }

      mapk_data <-mapk_data %>% tidyr::unite(Variant_Classification,3:ncol(mapk_data),sep = "_", remove = TRUE)
      print("Done..")



      colnames(mapk_data) <- c("CELLLINE", "Hugo_Symbol", "Variant_Classification")


      # use expression log2 values for any other multinomial regression to output the data

      r_expr_matrix <- expr_matrix_csv %>% dplyr::select(matches(rgenes_e))
      filtered_expression_matrix_pca <- filtered_expression_matrix_disease  %>% dplyr::select(matches(rgenes_e))

      # prepare GLM data and run it

      if (do_GLM == TRUE) {
        prepare_GLM_data(filtered_expression_matrix_pca,
                         rgenes_e,
                         mapk_data,
                         r_expr_matrix,
                         load_other_GLM,
                         GLM_all,
                         reg_type,
                         eXML_dir,
                         key,
                         disease_filename[j],
                         do_GLM,
                         inputs_dir,LM,folds,GLM_user_word)
      }

      if (carnival_flag == TRUE || load_other_GLM == TRUE) {
        # prepare Carnival and run it
        if (FEM_user == TRUE) {
          USER_EM <- filtered_expression_matrix_disease
        }
        else{
          USER_EM <- filtered_expression_matrix
        }
        mutations <- unique(cloned$Variant_Classification)
        
        
        ############### carry features even if WT ##############
        `%!in%` = Negate(`%in%`)
        names_mat_WT <- names_mat %>% dplyr::filter(CELLLINE %in% expr_matrix_csv)
        names_mat_WT <- names_mat %>% dplyr::filter(CELLLINE %!in% mapk_data$CELLLINE)
        
        names_mat_WT <- names_mat_WT %>% tidyr::unite(CCLE_ID,2:ncol(names_mat_WT),sep = "_", remove = TRUE)
        colnames(names_mat_WT) <- c("CELLLINE", "Variant_Classification")

        
        names_mat_WT$Variant_Classification <- paste0("WT_", names_mat_WT$Variant_Classification)
        
        
        prepare_Carnival(mapk_data,
                         USER_EM,
                         names_mat_WT,
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
                         key_opt_type,
                         mutations,source_code_dir,
                         epochs,lr,tts,update_it, bsu, hidden_dim, dropout, GNN_flag)

      }

      cloned_heatmap<- as.data.frame(cloned[,c(1,3,13:ncol(cloned))])
      cloned_heatmap <- merge(filtered_expression_matrix,cloned_heatmap,  all = TRUE)
      cloned_heatmap <- cloned_heatmap %>% relocate(Variant_Classification)

      cloned_heatmap$CELLLINE <- make.unique(as.character(cloned_heatmap$CELLLINE), sep = "_")
      rownames(cloned_heatmap) <- as.factor(cloned_heatmap$CELLLINE)

      cloned_heatmap$CELLLINE <- NULL
      colnames(cloned_heatmap) <- sapply(strsplit(colnames(cloned_heatmap), split='..', fixed=TRUE),function(x) (x[1]))
      cloned_heatmap$Variant_Classification[is.na(cloned_heatmap$Variant_Classification)] <- "WT"

      cloned_heatmap<-cloned_heatmap[order(cloned_heatmap[,'Variant_Classification']), ]
      colnames(cloned_heatmap) <- sapply(strsplit(colnames(cloned_heatmap), split='..', fixed=TRUE),function(x) (x[1]))
      order_rows <- as.factor(cloned_heatmap$Variant_Classification)
      scheme<- unique(cloned_heatmap$Variant_Classification)
      print("Creating heatmap...")

      write.csv(cloned_heatmap,"heatmap_data.csv")
      dend <-  heatmap.2(data.matrix(cloned_heatmap[,-c(1,2)]), scale = "column",col=bluered(100), breaks=seq(-3, 3, length.out=101),
                         keysize = 1, Rowv=FALSE,
                         trace = "none", density.info = "none",
                         key.title = "Log2 expression",RowSideColors=as.character(as.numeric(order_rows)),
                         main = "Expression (Log2 Z-transformed) for the regulon VS mutation type",
                         xlab = paste0("Regulon of ", violin_gene, " in ", "CCLE", disease_filename[j]),
                         font.lab = 40,ylab = NULL, margins = c(8,20))

      legend("topright",
             legend =  scheme ,
             col = unique(as.numeric(order_rows)),
             lty= 1,
             lwd = 5,
             cex=0.7
      )

      print(dend)
      #dev.off()


      setwd(results_dir)
      filtered_mutation_matrix_DL <- mut_matrix_csv %>% dplyr::filter(Tumor_Sample_Barcode %in% disease_cell_lines &
                                                                        !(Variant_Classification %in% "Silent") &
                                                                        !(isDeleterious %in% "False"))

      filtered_mutation_matrix_DL <- filtered_mutation_matrix_DL %>% dplyr::filter(stringr::str_detect(Hugo_Symbol, exact_genes))


      grouped_mutations <-  filtered_mutation_matrix %>% dplyr::group_by(Hugo_Symbol, Variant_Classification)
      grouped_mutations_DL <-  filtered_mutation_matrix_DL %>% dplyr::group_by(Hugo_Symbol, Variant_Classification)



      GM <- grouped_mutations %>% dplyr::summarize(count = n())
      GM_DL <- grouped_mutations_DL %>% dplyr::summarize(count = n())


      # extract now the pairs of gene-cell line, that is for each gene only the cell line in which it was found mutated

      filtered_mutation_matrix <- filtered_mutation_matrix %>% dplyr::filter(stringr::str_detect(Hugo_Symbol, exact_genes))
      pairs_GC_CH <- filtered_mutation_matrix %>% dplyr::select(c("Hugo_Symbol","Chrom", "Tumor_Sample_Barcode",
                                                                  "Variant_Classification","Codon_Change","Protein_Change","isDeleterious","isHotspot"))

      pairs_GC <- filtered_mutation_matrix %>% dplyr::select(c("Hugo_Symbol","Tumor_Sample_Barcode", "Variant_Classification"))
      ############################################################################################################

      filtered_cn_matrix <- cn_csv %>% dplyr::filter(cell_lines %in% filtered_expression_matrix$CELLLINE)

      filtered_cn_matrix <- filtered_cn_matrix  %>% dplyr::select(matches(genes_cn))
      mutated_genes <- filtered_mutation_matrix  %>% dplyr::select(matches("Hugo_Symbol"))
      mutated_genes <- unique(mutated_genes[,1])

      #############################################################################################################################
      # calcuate the log fold change of the  median of mutations for each gene across the median on all cell lines in disease)

      all_means_mutants <- c()
      all_means_wildtype <- c()
      all_means <- c()
      all_cn_means <- c()

      temp <- data.matrix(filtered_expression_matrix[,2:ncol(filtered_expression_matrix)])
      means_E <- apply(temp,2,mean) # applies function 'mean' to 2nd dimension (columns)
      medians_E <- apply(temp,2,median) # applies function 'mean' to 2nd dimension (columns)
      means_medians <- means_E - medians_E

      log_fold_change <- c()
      e_val_st_vector <- c()


      title_cancer = paste0("% done - Analyzing expression profiles in ", disease_filename[j])
      temp_genes <- c() # all the genes appearing mutated but also in the expression matrix
      number_of_mutations <- c() # the number of times a gene was found mutant (maybe in same cell line)
      number_of_unique_mutations <- c() # number of cell lines in which each gene was found mutated (unique)

      ################################################ M A I N   L O O P ##################################################
      for (i in 1:length(mutated_genes))        {

        lfc <- 0
        means_mutants <- 0
        means_wildtype <- 0
        mutant_expression_cell_lines <- c()
        wildtype_expression_cell_lines <- c()
        temp_column <- c()
        temp_gene0 <-mutated_genes[i]

        # make sure a gene found mutated also exists in expression matrix:
        temp_genes_final <- paste0(temp_gene0,"..")

        temp_mutated_genes_s <- paste0("^",temp_genes_final)

        temp_gene0_test <-  filtered_expression_matrix %>% dplyr::select(matches(temp_mutated_genes_s))

        if (nrow(temp_gene0_test)  != 0 ) {

          temp_gene <- paste0("^(", paste0(temp_gene0, ".."), ")")

          temp_gene_all <- paste0("^(", paste0("CELLLINE|",temp_gene), ")")

          temp_column <- filtered_expression_matrix  %>% dplyr::select(matches(temp_gene_all))

          if (temp_gene0 == violin_gene) {
            violin_column <- temp_column[order(temp_column[,'CELLLINE']), ]
          }

          if (ncol(temp_column) == 2)
          {

            temp_genes <- c(temp_genes, toString(temp_gene0)) # all the genes found mutated but also exist in expression matrix
            mutant_cell_lines0 <-   pairs_GC  %>% dplyr::filter(Hugo_Symbol %in% temp_gene0)

            mutant_cell_lines1 <- mutant_cell_lines0 %>% dplyr::select(matches("Tumor_Sample_Barcode|Variant_Classification"))
            mutant_cell_lines <- as.list(mutant_cell_lines1$Tumor_Sample_Barcode)

            number_of_mutations <- c(number_of_mutations,length(mutant_cell_lines0$Tumor_Sample_Barcode))
            number_of_unique_mutations <- c(number_of_unique_mutations,length(unique(mutant_cell_lines0$Tumor_Sample_Barcode)))

            mutant_expression_cell_lines <- temp_column  %>% dplyr::filter(CELLLINE %in% mutant_cell_lines0[,2])
            wildtype_expression_cell_lines <- temp_column  %>% dplyr::filter(!(CELLLINE %in% mutant_cell_lines0[,2]))

            #print("----------------------------------------------------")

            # convert to numerical to get means and subtract for log fold change
            mutant_expression_cell_lines <-  as.numeric(na.omit(mutant_expression_cell_lines[,2]))
            wildtype_expression_cell_lines <- as.numeric(na.omit(wildtype_expression_cell_lines[,2]))

            means_mutants <- mean(mutant_expression_cell_lines)
            means_wildtype <- mean(wildtype_expression_cell_lines)


            all_means_mutants = c(all_means_mutants, means_mutants)
            all_means_wildtype = c(all_means_wildtype, means_wildtype)

            # this is the expression levels of each gene across all disease-specific  cell lines

            gene_expression_column <- filtered_expression_matrix  %>% dplyr::select(matches(temp_gene))

            cn_temp_gene <- paste0("^", temp_gene0)
            cn_temp_gene <- paste0(cn_temp_gene,"..{2}")
            gene_cn_column <- filtered_cn_matrix  %>% dplyr::select(matches(cn_temp_gene))

            # decide which matrix (Expression or CN) has the less number of disease cell lines to base the test on that subset
            # so as the test can be performed (as equal entries required from both vectors)

            l1 <- as.numeric(unlist(gene_expression_column))
            l2 <- as.numeric(unlist(gene_cn_column))

            if (ncol(gene_cn_column) == 1) {

              all_means <- c(all_means,median(l1))
              all_cn_means <- c(all_cn_means,median(l2))
            }

            ##############################################################################################################

            #}
          }
        }
      } # end of mutated genes iterator loop

      ttest_matrix = cbind(all_means_mutants,all_means_wildtype)
      p_val_st <- 0
      p_val <- 0
      e_val_st <- 0
      pval <- 0
      mean_of_differences <- 0

      tt <- NULL

      tt <- t.test(ttest_matrix[,1], ttest_matrix[,2], paired = TRUE)
      p_val <- round(tt$p.value, digits = 5)
      mean_of_differences <- round(tt$estimate, digits = 5)


      #print("Spearman test between the median of a Gene expression against the mean of its CN across all disease cell lines:")
      stest_matrix = cbind(all_means,all_cn_means)

      ST <- cor.test(stest_matrix[,1],stest_matrix[,2] , method = "spearman", exact=FALSE )

      p_val_st <- round(ST$p.value, digits = 5)
      e_val_st <- round(ST$estimate, digits = 5)

      coverage = round((length(filtered_expression_matrix$CELLLINE)/length(disease_cell_lines))*100, digits = 1)

      amplified_genes = length(all_cn_means)

      log2_fold_change = ttest_matrix[,1] -  ttest_matrix[,2]


      main1 = paste0("(Log2) Fold change (MT-WL) of MAPK/ERK pathway genes in found mutant VS WT cell lines in ",disease_filename[j],
                     ".\n ----------------------------------Statistics---------------------------------------- \n | Genes in Pathway: ", length(genes_HGNC), paste0(". Regulons of ", violin_gene, ":"),
                     length(regulons_violin_gene),"."
                     ,     "\n | P-value of pairwise T-test: ", p_val, ". Mean of differences : ",
                     mean_of_differences , "\n | P-value of Spearman test (expression VS CN): ",p_val_st ,
                     ". Estimate of Spearman correlation (rho): ",  e_val_st, "\n | Number of Cell-Lines found: ",length(filtered_expression_matrix$CELLLINE) ,
                     "/", length(disease_cell_lines) , ". Lineage Coverage : ",  coverage, "%", " .Genes found mutated/amplified: ", length(all_means_mutants) ,"/",
                     amplified_genes, ".",
                     "\n | (source data-sets: DepMap Public ", dataset_version, ")", "\n ---------------------------------------------------------------------------------")


      gene_labels <- paste0(as.list(temp_genes), "(", as.list(number_of_unique_mutations),"/",length(filtered_expression_matrix$CELLLINE),")")

      log2_fold_change_df <- data.frame(temp_genes,number_of_mutations,number_of_unique_mutations,
                                        ttest_matrix[,1],ttest_matrix[,2],log2_fold_change,gene_labels)
      names(log2_fold_change_df) <- c("Gene", "Number_of_Mutations", "Number_of_Unique_Mutations", "Mutant",
                                      "WT", "Log2_Fold_difference", "LABELS")

      log2_fold_change_df <- log2_fold_change_df %>% arrange(desc(number_of_mutations))

      Number_of_Mutations <- log2_fold_change_df$Number_of_Mutations

      ########################################### MUTATIONS OF ALL GENES GRAPH and T-test #################################################
      # just get the labels as : name of gene + (number_of_unique_mutations)


      sp2 <- ggplot(log2_fold_change_df,  aes(x=factor(Gene, levels=Gene), y=Log2_Fold_difference,label = LABELS,size = Number_of_Mutations)) +

        geom_point(color = dplyr::case_when(log2_fold_change_df$Log2_Fold_difference > 1 ~ "#FF0000",
                                            log2_fold_change_df$Log2_Fold_difference < -1 ~ "#FF0000",
                                            TRUE ~ "#00CC00"), alpha = 0.8) +
        geom_hline(yintercept = 1,linetype="dotted") +
        geom_hline(yintercept = -1,linetype="dotted") +

        geom_hline(yintercept = 2,linetype="dashed") +
        geom_hline(yintercept = -2,linetype="dashed") +
        geom_text_repel(label = ifelse(log2_fold_change_df$Number_of_Mutations >= 1 &
                                         log2_fold_change_df$Log2_Fold_difference < 1 &
                                         log2_fold_change_df$Log2_Fold_difference > -1, as.character(log2_fold_change_df$LABELS) , "" ), size = 2) +
        geom_text_repel( data          = subset(log2_fold_change_df, Log2_Fold_difference > 1),
                         nudge_y       = 16 - subset(log2_fold_change_df, Log2_Fold_difference > 1)$Log2_Fold_difference,
                         size          = 2,
                         box.padding   = 1.5,
                         point.padding = 0.5,
                         force         = 0.5,
                         segment.size  = 0.5,
                         segment.color = "grey50",
                         direction     = "y") +
        geom_label_repel(data         = subset(log2_fold_change_df, Log2_Fold_difference < -1),
                         nudge_y       = -16 - subset(log2_fold_change_df, Log2_Fold_difference < -1)$Log2_Fold_difference,
                         size          = 2,
                         box.padding   = 0.5,
                         point.padding = 0.5,
                         force         = 0.5,
                         segment.size  = 0.5,
                         segment.color = "grey50",
                         direction     = "y") +

        scale_x_discrete(expand = expand_scale(mult = c(0.005, .05))) +
        scale_y_continuous(expand = expand_scale(mult = c(0.005, .01)))  +
        theme(axis.text.x=element_text(size=5, angle=90,hjust=0.95,vjust=0.2),plot.title = element_text(size = 9)) + ggtitle(main1)




      #write.csv(pairs_GC_CH,"pairs_GC_CH.csv")
      pairs_GC_alp <-   pairs_GC_CH[order(pairs_GC_CH[,'Tumor_Sample_Barcode']), ]

      # for violin_gene
      pairs_GC_alp_bkp <-  pairs_GC_alp
      pairs_GC_alp <- pairs_GC_alp %>% dplyr::filter(Hugo_Symbol %in% violin_gene & Tumor_Sample_Barcode %in% violin_column$CELLLINE)


      pairs_GC_alp <- add_column(pairs_GC_alp, Expression = 0)
      pairs_GC_CH <- add_column(pairs_GC_CH, Expression = 0)

      pairs_GC_alp  <- pairs_GC_alp   %>% dplyr::mutate(Hugo_Symbol = disease_filename[j])

      # now pairs_GC_alph has all cellines, mutation types and corresponding expression levels for the violin gene

      if (!is.null(violin_column)) {

        violin_column <- violin_column %>% dplyr::filter(CELLLINE %in% pairs_GC_alp$Tumor_Sample_Barcode)
        # for-loop to fill-in the expression levels of all mutated cell lines for the "Violin-Gene" (usually TP53)
        j2 = 1
        for (i in 1: length(pairs_GC_alp$Tumor_Sample_Barcode)) {

          index = match(pairs_GC_alp$Tumor_Sample_Barcode[i],violin_column$CELLLINE)
          pairs_GC_alp$Expression[j2] = violin_column[index,2]

          j2 = j2 + 1
        }
      }
      # same but for all genes
      j2 = 1

      print("Analyzing mutation profiles...")

      for (i in 1: length(pairs_GC_CH$Tumor_Sample_Barcode)) {
        temp_gene <- pairs_GC_CH$Hugo_Symbol[i]

        temp_gene_all <- paste0("^(", paste0("CELLLINE|",temp_gene), ")")


        temp_column <- NULL
        temp_column <- filtered_expression_matrix  %>% dplyr::select(matches(temp_gene_all))
        if (ncol(temp_column) == 2)     {

          index = match(pairs_GC_CH$Tumor_Sample_Barcode[i],temp_column$CELLLINE)

          pairs_GC_CH$Expression[j2] = temp_column[index,2]
        }
        else{
          pairs_GC_CH$Expression[j2] = NA

        }

        j2 = j2 + 1
      }

      ###### (Box - Violin - Scatter Plot) on all mutated genes across cell lines and types of mutations versus expression level ##########
      dodge <- position_dodge(width = 0.4)


      dodge <- position_dodge(width = 0.4)

      main6_1 = paste0("Violin and boxplots of types of mutations across all mutated genes found in ",
                       disease_filename[j]," \n (source data-sets: DepMap Public  ", dataset_version, ")")
      main6_2 = paste0("Violin and boxplots of types of mutations (annotated with Gene Hugo Symbols) across all mutated genes found in ",
                       disease_filename[j]," \n (source data-sets: DepMap Public  ", dataset_version, ")")
      main6_3 = paste0("Violin and boxplots of types of mutations (annotated with cell line names) across all mutated genes found in ",
                       disease_filename[j]," \n (source data-sets: DepMap Public  ", dataset_version, ")")

      sp6_1 <- ggplot(data = pairs_GC_CH,aes(x = Variant_Classification, y = Expression , fill = Variant_Classification))+
        #scale_fill_viridis_d( option = "D")+

        geom_boxplot(width=.1,notch = FALSE,  outlier.size = 0, color="black",lwd=1.2, alpha = 0.7, position = dodge) +
        geom_point( shape = 21,size=2, position = dodge, color="black",alpha=1) +

        geom_violin(alpha=0.2,position = dodge,trim= FALSE) +
        ylab(  c("Expression (log2 values)")  )  +
        xlab(  c(paste0("Mutation variation in ", disease_filename[j])) ) +

        font("xylab",size=20)+
        font("xy",size=20)+
        font("xy.text", size = 20) +
        font("legend.text",size = 20) +
        theme(axis.text.x=element_text(size=25, angle=90,hjust=0.95,vjust=0.02)) +
        ggtitle(main6_1)

      sp6_2 <- ggplot(data = pairs_GC_CH,aes(x = Variant_Classification, y = Expression , fill = Variant_Classification))+
        #scale_fill_viridis_d( option = "D")+

        geom_boxplot(width=.1,notch = FALSE,  outlier.size = 0, color="black",lwd=1.2, alpha = 0.7, position = dodge) +
        geom_point( shape = 21,size=2, position = dodge, color="black",alpha=1) +
        geom_label_repel(aes(label=pairs_GC_CH$Hugo_Symbol),
                         box.padding   = 0.5,
                         point.padding = 0.005,
                         segment.color = 'grey50', size = 1.5) +
        geom_violin(alpha=0.2,position = dodge,trim= FALSE) +
        ylab(  c("Expression (log2 values)")  )  +
        xlab(  c(paste0("Mutation variation in ", disease_filename[j])) ) +

        font("xylab",size=20)+
        font("xy",size=20)+
        font("xy.text", size = 20) +
        font("legend.text",size = 20) +
        theme(axis.text.x=element_text(size=25, angle=90,hjust=0.95,vjust=0.02)) +
        ggtitle(main6_2)



      sp6_3 <- ggplot(data = pairs_GC_CH,aes(x = Variant_Classification, y = Expression , fill = Variant_Classification))+
        #scale_fill_viridis_d( option = "D")+

        geom_boxplot(width=.1,notch = FALSE,  outlier.size = 0, color="black",lwd=1.2, alpha = 0.7, position = dodge) +
        geom_point( shape = 21,size=2, position = dodge, color="black",alpha=1) +
        geom_label_repel(aes(label=pairs_GC_CH$Tumor_Sample),
                         box.padding   = 0.5,
                         point.padding = 0.005,
                         segment.color = 'grey50', size = 1) +
        geom_violin(alpha=0.2,position = dodge,trim= FALSE) +
        ylab(  c("Expression (log2 values)")  )  +
        xlab(  c(paste0("Mutation variation in ", disease_filename[j])) ) +

        font("xylab",size=20)+
        font("xy",size=20)+
        font("xy.text", size = 20) +
        font("legend.text",size = 20) +
        theme(axis.text.x=element_text(size=25, angle=90,hjust=0.95,vjust=0.02))+
        ggtitle(main6_3)



      ###### (Box - Violin - Scatter Plot) on violin_gene across cell lines and types of mutations versus expression level ##########
      genes <- paste0(violin_gene,"\\s*?\\.{2}",collapse="|")

      genes0 <- paste0("\\b",violin_gene,"\\s*?\\.{2}\\b", collapse="|")
      genes_e <- paste0("^(", paste0("CELLLINE|",genes), ")")
      genes_cn <- paste0("^(", paste0("cell_lines|",genes0), ")")

      test_data <- data.frame()

      test_data <- expr_matrix_csv  %>% dplyr::select(matches(genes_e))
      test_data$WT <- ifelse(expr_matrix_csv$CELLLINE %in% pairs_GC_alp$Tumor_Sample_Barcode, "MT", "WT")
      colnames( test_data) <- c("cell_lines", "Expression_log2", "State")
      #write.csv(test_data,"test_data.csv")
      #print(test_data)
      # copy number
      print("Analyzing CNV profiles...")
      cn_gene <- data.frame()
      cn_gene <- cn_csv  %>% dplyr::select(matches(genes_cn))

      if (ncol(cn_gene) == 1) # no CN data was found for the violin gene
      {
        cn_gene[,2] <- rep(0,nrow(cn_gene)) # filling then CN data with 0 zeroes
        print("no CNV data found for the selected gene")

      }
      else
      {
        print("CNV data found for the selected gene")
        # plot CN across Expression levels for violin gene only
        violin_gene_CN <- NULL
        violin_gene_E <- NULL

        violin_gene_CN <- cn_gene

        violin_gene_E <- expr_matrix_csv %>% dplyr::select(matches(genes_e))

        if (ncol(violin_gene_E)==2 & ncol(violin_gene_CN)==2) {
          colnames(violin_gene_E) <- c("CELLLINE",violin_gene)
          colnames(violin_gene_CN) <- c("CELLLINE",violin_gene)

          #violin_gene_CN_E <- NULL

          violin_gene_CN_E <- merge(violin_gene_CN,violin_gene_E, by = "CELLLINE")
          # print(violin_gene_CN_E)
          colnames(violin_gene_CN_E) <- c("CELLLINE","CN_log2","Expression_log2")
          violin_gene_CN_E <- violin_gene_CN_E %>% dplyr::select("CN_log2","Expression_log2")


          violin_gene_CN_E_plot <- ggplot(violin_gene_CN_E, aes(x=CN_log2, y=Expression_log2))+
            scale_x_continuous(breaks=seq(0,80,1)) +
            geom_point()+
            ggtitle(paste0("Copy number (log2) versus expression (log2) for ", violin_gene," in all cell lines"))
          print(violin_gene_CN_E_plot)
        }

      }
      print("Processing final graphs...")
      cn_gene <-  cn_gene  %>% dplyr::select(c(1,2))
      colnames(cn_gene) <- c("cell_lines", "CN")
      cn_gene$CN_S <- ifelse(cn_gene$CN >= 1, "Amplification", "Deletion")
      colnames(cn_gene) <- c("cell_lines", "CN", "GYSTIC_CALLS")
      cn_gene <-  cn_gene  %>% dplyr::filter(cell_lines %in% disease_cell_lines)
      write.csv(cn_gene,"cn_gene.csv")

      # merging dataframes to one
      total_o <- merge(test_data, cn_gene, by = "cell_lines")

      pairs_GC_alp_temp <- NULL
      pairs_GC_alp_temp <- pairs_GC_alp
      colnames(pairs_GC_alp_temp)[3] <- "cell_lines"

      pairs_GC_alp_temp<- pairs_GC_alp_temp %>% dplyr::select("cell_lines","Variant_Classification","Codon_Change","Protein_Change","isDeleterious", "isHotspot")

      total <- merge(total_o, pairs_GC_alp_temp, by = "cell_lines", all = TRUE)

      total <- total[!is.na(total$State),]

      total  <- mutate(total, Variant_Classification = ifelse(is.na(Variant_Classification), "WT", as.character(Variant_Classification)))

      write.csv(total,"total.csv")



      if (!is.null(violin_column)) {
        print("Printing violin plots now...")
        ####### is Hotspot ######
        dodge <- position_dodge(width = .4)
        main7 = paste0("Violin and boxplots of WT vesus Mutant expression levels for ", violin_gene, " in all cell lines",
                       " \n (source data-sets: DepMap Public  ", dataset_version, ")")

        total2 <- total[!is.na(total$GYSTIC_CALLS),]

        sp777 <- ggplot(data = total2,aes(x = State, y = Expression_log2, fill = isHotspot))+
          #scale_fill_viridis_d( option = "D")+
          geom_boxplot(width=.1,notch = FALSE,  outlier.size = 0, color="black",lwd=1.2, alpha = 0.7, position = dodge) +
          geom_point(shape = 21,size=2, position = dodge, color="black",alpha=1) +
          #geom_label_repel(aes(label=total$cell_lines),
          #                 box.padding   = 0.5,
          #                 point.padding = 0.005, size = 1.8) +
          geom_violin(alpha=0.2,position = dodge,trim= FALSE) +
          ylab(  c("Expression (log2 values)")  )  +
          xlab(  c(paste0(violin_gene," WT/MT Classification in ", disease_filename[j])) ) +
          font("xylab",size=20)+
          font("xy",size=20)+
          font("xy.text", size = 20) +
          font("legend.text",size = 20) +
          theme(axis.text.x=element_text(size=25, angle=90,hjust=0.95,vjust=0.02))+
          ggtitle(main7) +
          #stat_compare_means(method = "anova", label.y = 20, size = 5)  +
          #stat_compare_means(label.y = 21, size = 5) +
          stat_compare_means(method = "t.test",label.y = 12, size = 7) +
          stat_compare_means(method = "wilcox.test",label.y = 14, size = 7)

        print(sp777)
        plot_list <- c(plot_list,sp777)
        ggsave(filename="CCLE_WT_MT_ishotspot.png", plot=sp777)
        ####### is deleterious ######
        dodge <- position_dodge(width = .4)
        main7 = paste0("Violin and boxplots of WT vesus Mutant expression levels for ", violin_gene, " in all cell lines",
                       " \n (source data-sets: DepMap Public  ", dataset_version, ")")



        sp777 <- ggplot(data = total2,aes(x = State, y = Expression_log2, fill = isDeleterious))+
          #scale_fill_viridis_d( option = "D")+
          geom_boxplot(width=.1,notch = FALSE,  outlier.size = 0, color="black",lwd=1.2, alpha = 0.7, position = dodge) +
          geom_point(shape = 21,size=2, position = dodge, color="black",alpha=1) +
          #geom_label_repel(aes(label=total$cell_lines),
          #                 box.padding   = 0.5,
          #                 point.padding = 0.005, size = 1.8) +
          geom_violin(alpha=0.2,position = dodge,trim= FALSE) +
          ylab(  c("Expression (log2 values)")  )  +
          xlab(  c(paste0(violin_gene," WT/MT Classification in ", disease_filename[j])) ) +
          font("xylab",size=20)+
          font("xy",size=20)+
          font("xy.text", size = 20) +
          font("legend.text",size = 25) +
          theme(axis.text.x=element_text(size=25, angle=90,hjust=0.95,vjust=0.02))+
          ggtitle(main7) +
          #stat_compare_means(method = "anova", label.y = 20, size = 5)  +
          #stat_compare_means(label.y = 21, size = 5) +
          stat_compare_means(method = "t.test",label.y = 12, size = 7) +
          stat_compare_means(method = "wilcox.test",label.y = 14, size = 7)

        print(sp777)
        plot_list <- c(plot_list,sp777)
        ggsave(filename="CCLE_WT_MT_isdel.png", plot=sp777)

        dodge <- position_dodge(width = .4)
        main7 = paste0("Violin and boxplots of WT vesus Mutant expression levels for ", violin_gene, " in all cell lines",
                       " \n (source data-sets: DepMap Public  ", dataset_version, ")")

        total2 <- total[!is.na(total$GYSTIC_CALLS),]

        sp777 <- ggplot(data = total2,aes(x = State,y = Expression_log2, fill = GYSTIC_CALLS))+
          #scale_fill_viridis_d( option = "D")+
          geom_boxplot(width=.1,notch = TRUE,  outlier.size = 0, color="black",lwd=1.2, alpha = 0.7, position = dodge) +
          geom_point(shape = 21,size=2, position = dodge, color="black",alpha=1) +
          #geom_label_repel(aes(label=total$cell_lines),
          #                 box.padding   = 0.5,
          #                 point.padding = 0.005, size = 1.8) +
          geom_violin(alpha=0.2,position = dodge,trim= FALSE) +
          ylab(  c("Expression (log2 values)")  )  +
          xlab(  c(paste0(violin_gene," WT/MT Classification in ", disease_filename[j])) ) +
          font("xylab",size=20)+
          font("xy",size=20)+
          font("xy.text", size = 20) +
          font("legend.text",size = 20) +
          theme(axis.text.x=element_text(size=25, angle=90,hjust=0.95,vjust=0.02))+
          ggtitle(main7) +
          stat_compare_means(method = "anova", label.y = 12, size = 8)  +
          stat_compare_means(label.y = 15, size = 8)


        print(sp777)
        plot_list <- c(plot_list,sp777)
        ggsave(filename="CCLE_WT_MT_gystic.png", plot=sp777)
        ####### is deleterious ######
        dodge <- position_dodge(width = .4)
        main7 = paste0("Violin and boxplots of WT vesus Mutant expression levels for ", violin_gene, " in all cell lines",
                       " \n (source data-sets: DepMap Public  ", dataset_version, ")")





      }
      if (!is.null(violin_column)) {
        print("Printing violin plots II...")
        dodge <- position_dodge(width = .4)
        main7 = paste0("Violin and boxplots of WT vesus Mutant expression levels and types for ", violin_gene, " in all cell lines",
                       " \n (source data-sets: DepMap Public  ", dataset_version, ")")
        #my_comparisons <- list(  c("Missense","Nonsense","In_Frame_Del","Splice_Site","Frame_Shift_Ins","Frame_Shift_Del","In_Frame_Ins", "WT"))
        #scheme <- c("Missense","Nonsense","In_Frame_Del","Splice_Site","Frame_Shift_Ins","Frame_Shift_Del","In_Frame_Ins", "WT")
        total2 <- total[!is.na(total$GYSTIC_CALLS),]
        sp77 <- ggplot(data = total2,aes(x = Variant_Classification, y = Expression_log2, fill = Variant_Classification))+
          #scale_fill_viridis_d( option = "D")+
          geom_boxplot(width=.1,notch = FALSE,  outlier.size = 0, color="black",lwd=1.2, alpha = 0.7, position = dodge) +
          #geom_point(shape = ifelse(total$CN_S == "Amplification", 8, 21),size=ifelse(total$CN_S == "Amplification", 3, 1),
          #           position = dodge, color= ifelse(total$CN_S == "Amplification", "red", "yellow"),alpha=1) +
          #geom_label_repel(aes(label=total$cell_lines),
          #                 box.padding   = 0.5,
          #                 point.padding = 0.005, size = 1.8) +
          geom_violin(alpha=0.2,position = dodge,trim= FALSE) +
          ylab(  c("Expression (log2 values)")  )  +
          xlab(  c(paste0(violin_gene," WT/MT Classification in ", disease_filename[j])) ) +
          font("xylab",size=20)+
          font("xy",size=20)+
          font("xy.text", size = 20) +
          font("legend.text",size = 20) +
          theme(axis.text.x=element_text(size=25, angle=90,hjust=0.95,vjust=0.02))+
          ggtitle(main7) +
          stat_compare_means(method = "anova", label.y = 11,label.x = 2, size = 8)  +
          stat_compare_means(label.y = 10, label.x = 2, size = 8)# +
        #stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.", label.y=20, size = 5)

        print(sp77)
        plot_list <- c(plot_list,sp77)
        ggsave(filename="CCLE.png", plot=sp77)

      }

      if (!is.null(violin_column)) {
        print("Printing violin plots III...")
        dodge <- position_dodge(width = .4)
        main7 = paste0("Violin and boxplots of types of mutations for ", violin_gene, " in ",
                       disease_filename[j]," \n (source data-sets: DepMap Public  ", dataset_version, ")")


        colnames(names_mat)[1] <- "Tumor_Sample_Barcode"
        pairs_GC_alp <- merge(names_mat, pairs_GC_alp, by = "Tumor_Sample_Barcode")

        sp7 <- ggplot(data = pairs_GC_alp,aes(x = Variant_Classification, y = Expression , fill = Variant_Classification))+
          #scale_fill_viridis_d( option = "D")+
          geom_boxplot(width=.1,notch = FALSE,  outlier.size = 0, color="black",lwd=1.2, alpha = 0.7, position = dodge) +
          geom_point(shape = 21,size=2, position = dodge, color="black",alpha=1) +
          geom_label_repel(aes(label=pairs_GC_alp$CCLE_ID),
                           box.padding   = 0.5,
                           point.padding = 0.005, size = 3) +
          geom_violin(alpha=0.2,position = dodge,trim= FALSE) +
          ylab(  c("Expression (log2 values)")  )  +
          xlab(  c(paste0(violin_gene," Variant Classification in ", disease_filename[j])) ) +
          font("xylab",size=20)+
          font("xy",size=20)+
          font("xy.text", size = 20) +
          font("legend.text",size = 20) +
          theme(axis.text.x=element_text(size=25, angle=90,hjust=0.95,vjust=0.02))+
          ggtitle(main7)
        print(sp7)
        plot_list <- c(plot_list,sp7)

      }
      all_cancers <- rbind(all_cancers,pairs_GC_alp)

      print("Finalising...")
      #pairs_GC_CH_bkp <- NULL
      pairs_GC_CH_bkp <- pairs_GC_CH

      # map broad ids to cell line ids
      colnames(names_mat)[1] <- "Tumor_Sample_Barcode"
      pairs_GC_CH_bkp <- merge(names_mat, pairs_GC_CH_bkp, by = "Tumor_Sample_Barcode")
      pairs_GC_CH_bkp <- pairs_GC_CH_bkp %>% dplyr::select(c("Hugo_Symbol","Chrom", "CCLE_ID",
                                                             "Variant_Classification", "Codon_Change" , "Protein_Change" , "Expression"))
      colnames(pairs_GC_CH_bkp)[3] <- "Tumor_Sample_Barcode"

      all_cancers_all_genes <- rbind(all_cancers_all_genes,pairs_GC_CH)

      print(sp2)

      print(sp6_1)
      plot_list <- c(plot_list,sp2,sp6_1)

      pairs_GC_CH  <-  pairs_GC_CH %>% dplyr::group_by(Hugo_Symbol, Variant_Classification)
      pairs_GC_CH <- pairs_GC_CH %>% dplyr::summarize(count = n())

      pairs_GC_CH <- pairs_GC_CH  %>% spread(key=Variant_Classification, value=count)
      #pairs_GC_CH[is.na(pairs_GC_CH)] <- 0
      pairs_GC_CH$Hugo_Symbol <- NULL
      pairs_GC_CH <- pairs_GC_CH[2:nrow(pairs_GC_CH),]



      dev.off()
    } # end of cancer type iterator loop

    graphics.off()

  }

  print("Finished scan.")
  end <- Sys.time()
  time <- end - start
  print(time)
  print("---------------END OF ANALYSIS--------------")

  # end of function
}
