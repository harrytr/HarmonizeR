#CARNIVAL module

Carnival_opt <-function(iterator_index,
                        df_EM,
                        results_dir,
                        inputs_dir,
                        disease_filename_j,
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
                        radar_plot_data,key_opt,key_opt_type, mutations, source_code_dir,
                        epochs,lr,tts,update_it,bsu, hidden_dim, dropout, GNN_flag) {
  
  library(progress)
  library(dorothea)
  library(progeny)
  
  mutations <- c(mutations, "WT")
  
  if (.Platform$OS.type == "unix") {
    cplex_path = "/Applications/CPLEX_Studio221/cplex/bin/x86-64_osx/cplex"
  }
  else if  (.Platform$OS.type == "windows") {
    cplex_path <- "C:/Program Files/IBM/ILOG/CPLEX_Studio221/cplex/bin/x64_win64/"
  }
  #  PREPARE CARNIVAL INPUT FILES:
  print("Preparing the viper regulon list...")
  
  print(paste0("Standarizing Expression Matrix for the regulons of ", violin_gene, " ..."))
  
  df_EM <- scale(df_EM)
  df_EM <- as.data.frame(df_EM)
  
  # CREATE THE INPUT REGULON FOR CARNIVAL
  print("Creating the VIPER regulon object for CARNIVAL...")
  
  
  #=======================================================================
  #load(file = paste0(inputs_dir,"/omni.RData"))
  
  print("Getting data from OmniPathR...")
  regulon_df <- import_dorothea_interactions(organism=9606)
  save(regulon_df, file = "dorothea_omnipath_interactions.RData")
  #regulon_df <- omnipath_interactions(organism = 9606)
  regulon_df <- regulon_df[which(regulon_df$is_directed==1), ] #keeping only directed
  regulon_df <- regulon_df[which((regulon_df$is_stimulation+regulon_df$is_inhibition)==1), ] #keeping only regulons which are either activations/inhibitions
  
  df <- matrix(data = , nrow = nrow(regulon_df), ncol = 3) #creating the regulon dataframe for the createRegulonList function
  
  df[, 1] = regulon_df$source_genesymbol
  df[, 3] = regulon_df$target_genesymbol
  
  df[which(regulon_df$is_stimulation==1), 2] <- "1" #assigning correct sign
  df[which(regulon_df$is_inhibition==1), 2] <- "-1"
  
  print("Creating the network structure...")
  colnames(df) = c("Source", "Sign", "Target")
  df <- as.data.frame(unique(df))
  df$Source = as.character(df$Source)
  df$Sign = as.numeric(as.character(df$Sign))
  df$Target = as.character(df$Target)
  
  network <- df
  print("Saving PKN...")
  write.csv(network,paste0(results_dir,"/interactions.csv"))
  #save(network, file="network.RData")
  
  # regulon_object = createRegulonList(regulon_table = df)
  #
  # save(regulon_object, file = "regulon_object.RData")
  setwd(results_dir)
  #load(file = paste0(inputs_dir,"/regulon_A_B_C_D_E.RData"))
  #load(file = paste0(inputs_dir,"/regulon_A_B_C.RData"))
  
  #regulon_A_B_C <- regulon_A_B_C[which(regulon_A_B_C$tfregulons_level%in%c("A")), ]
  #############################################################################
  dir.create(path = paste0(results_dir,"/measurements"))
  
  
  carnival_path_1 <- paste0(results_dir)
  
  print(carnival_path_1)
  #dir.create(path = carnival_path_1)
  carnival_path <- paste0(carnival_path_1,"/opt/")
  dir.create(path = carnival_path)
  
  
  setwd(paste0(results_dir,"/measurements"))
  save(network, file="network.RData")
  tp_input_i <- -1
  tp_input_a <- 1
  input_df_i  <- as.data.frame(tp_input_i)
  input_df_a  <- as.data.frame(tp_input_a)
  colnames(input_df_i) <- violin_gene
  colnames(input_df_a) <- violin_gene
  save(input_df_i, file="inputs_i.RData")
  save(input_df_a, file="inputs_a.RData")
  
  
  load(file = paste0(inputs_dir,"/dorothea_hs_pancancer.RData"))
  dorothea_hs_pancancer <- dorothea_hs_pancancer %>% dplyr::filter(confidence %in% c("A","B","C","D","E"))
  print("Calculating TF activities...")
  TF_activities  = as.data.frame(dorothea::run_viper(df_EM, dorothea_hs_pancancer,
                                                     options = list(method = "none", minsize = 1,  nes = T,
                                                                    eset.filter = FALSE, cores = 1,
                                                                    verbose = T)))
  
  print("Saving the regulon activities from Viper/DoRothEA...")
  save(TF_activities, file="measurements.RData")
  
  tfList <- generateDataframeTF(df = TF_activities, top = top_user, 1:ncol(df_EM))
  
  resList <- list()
  list_of_graphs <- c()
  graph_heatmap <- matrix(data = , nrow = length(tfList), ncol = length(tfList))
  
  tfList_names <- names(tfList)
  
  graph_heatmap_ID <- matrix(data = , nrow = length(tfList), ncol = length(tfList)) # compare same type of mutation
  graph_heatmap_DEL <- matrix(data = , nrow = length(tfList), ncol = length(tfList)) # compare deleterious or not
  graph_heatmap_hotspot <- matrix(data = , nrow = length(tfList), ncol = length(tfList)) # compare hotspot or not
  graph_heatmap_KEY <- matrix(data = , nrow = length(tfList), ncol = length(tfList)) # user comparison
  
  top_similar_1 <- c()
  top_similar_2 <- c()
  top_bc_network <-  c()
  top_bc_vertex <- c()
  top_similar_score <- c()
  total_Bar <- length(tfList)
  pb <- progress_bar$new(total = total_Bar)
  
  print("optimising networks with CARNIVAL MILP...")
  
  
  for (ii in 1: length(tfList)){
    pb$tick()
    Sys.sleep(1 / total_Bar)
    temp_dir <- paste0(carnival_path, "//", names(tfList)[ii])
    print(temp_dir)
    dir.create(path = temp_dir)
    
    
    
    final_input <- NULL
    final_input <- ifelse(str_detect(names(tfList)[ii],"NOTisDeleterious|WT"),FALSE,TRUE)
    
    setwd(temp_dir)
    print(temp_dir)
    dot_file <- paste0(temp_dir,"//network_solution.dot")
    # !file.exists(dot_file)
    
    if (length(list.files(temp_dir)) == 0 ){
      if (!(is.null(final_input) )){
        try(
          # res <- runCARNIVAL(inputObj = NULL, measObj = tfList[[ii]], netObj = network,weightObj = NULL,
          #                                           solverPath = cplex_path , solver = "cplex",
          #                                          dir_name = temp_dir, mipGAP = GAP, threads = cpu_threads,limitPop = 10, poolCap = 10)
        )
        if (final_input == TRUE) {
          print(paste0(violin_gene," knockdown..."))
          try(
            res <- runCARNIVAL(inputObj = input_df_i, measObj = tfList[[ii]], netObj = network,weightObj = NULL,
                               solverPath = cplex_path , solver = "cplex",
                               dir_name = temp_dir, mipGAP = GAP, threads = cpu_threads,limitPop = 10, poolCap = 10)
          )
        }
        else {
          print(paste0(violin_gene," stimulation..."))
          try(
            res <- runCARNIVAL(inputObj = input_df_a, measObj = tfList[[ii]], netObj = network,weightObj = NULL,
                               solverPath = cplex_path , solver = "cplex",
                               dir_name = temp_dir, mipGAP = GAP, threads = cpu_threads,limitPop = 10, poolCap = 10)
          )
        }
      }
      
      
    }
    else {print("Network has already been optimized ; folder not empty! Continuing...")}
  }
  
  all_networks <- list()
  pb <- progress_bar$new(total = total_Bar)
  
  range_opt <- length(tfList)
  
  #### RUN PYTHON SCRIPT TO CONVERT .dot to .graphml so that igraph can read it
  setwd(source_code_dir)
  print("Running now Python script to convert .DOT to igraph...")
  
  if (.Platform$OS.type == "unix") {
    try(
      command_1 <- paste0("python3 dot_2_igraph.py ",'"',carnival_path,'"'," ",'"',violin_gene,'"')
    )
  }
  
  else if  (.Platform$OS.type == "windows") {
    try(
      command_1 <- paste0("python dot_2_igraph.py ",'"',carnival_path,'"'," ",'"',violin_gene,'"')
    )
  }
  
  print(command_1)
  
  system(command_1)
  print("Done!")
  
  #############################################################################
  
  labels_csv <-  matrix(data = , nrow = length(tfList_names), ncol = 3)
  labels_csv <- as.data.frame(labels_csv)
  colnames(labels_csv) <- c("filename","mutation","labels")
  print("Read all possible networks as graphs...")
  for (i in 1:  range_opt) {
    
    pb$tick()
    Sys.sleep(1 / total_Bar)
    net_base <- NULL
    temp_dir <- paste0(carnival_path, "//", names(tfList)[i])
    base_file2 <- paste0(temp_dir,"//network_solution.dot")
    base_file <- paste0(temp_dir,"//network_solution.graphml")
    
    new_name_gml <- paste0(temp_dir,"//", as.character(names(tfList)[i]),".graphml")
    if (file.exists(new_name_gml)){
      
      print("Removing old graphml file")
      try(
        file.remove(new_name_gml)
      )
    }
    
    if (file.exists(base_file)){
      
      print("Renaming graphml file")
      try(
        file.rename(base_file, new_name_gml)
      )
      
      
      labels_csv[i,1] <- paste0(as.character(names(tfList)[i]),".graphml")
      labels_csv[i,2] <- names(tfList)[i]
      net_base <- NULL
      print("Reading graphml file:")
      net_base <- igraph::read_graph(new_name_gml, format = "graphml")
      base <- net_base
      #my_data <- read.delim(base_file2)
      if (is.null(net_base) == FALSE) {
        all_networks[[i]] <- net_base
      }
      else{
        print("Graphml file could not be read through Python!")
        Sys.sleep(3)
        labels_csv[i,1] <- paste0(as.character(i),".graphml")
        labels_csv[i,2] <- names(tfList)[i]
        net_base <- NULL
      }
    }
    else
    {
      print("Graphml file not found!")
      Sys.sleep(3)
      labels_csv[i,1] <- paste0(as.character(i),".graphml")
      labels_csv[i,2] <- names(tfList)[i]
      net_base <- NULL
    }
    
    # ================ COMMUNITY DETECTION =======================
    
    if (is.null(net_base) == FALSE) {
      

      print("Converting to visNetwork object")
      GRN <- toVisNetworkData(net_base)
      GRN$nodes$name <- GRN$nodes$label
      GRN$nodes$color <- GRN$nodes$fillcolor
      
      GRN$nodes$label <- ifelse(GRN$nodes$fillcolor == "lavender", paste0(GRN$nodes$label,"up"),paste0(GRN$nodes$label,"down"))
      GRN$nodes$shape <- ifelse(GRN$nodes$label %in% c(paste0(violin_gene,"up"),paste0(violin_gene,"down")), "star", "circle")
      GRN$edges$id <- NULL
      print("Visualising ...")
      print(paste0(carnival_path, "/", names(tfList)[i],"/GRAPH",".html"))

      mainP = paste0("Optimized network for:: ",names(tfList)[i],"::" , violin_gene)
      sp_g_6 <- visNetwork(GRN$nodes,GRN$edges, main = mainP, height = "700px", width = "100%") %>%
        visEdges(labelHighlightBold= "TRUE",arrows = "to") %>%
        visInteraction(zoomView = TRUE) %>%
        visOptions(highlightNearest = TRUE,nodesIdSelection = TRUE) %>%
        visPhysics(stabilization = FALSE)  %>%
        visEdges(smooth = FALSE) %>% visLayout(hierarchical = TRUE)
      try(
        visSave(sp_g_6, file = paste0(carnival_path, "/", names(tfList)[i],"/GRAPH",".html"), selfcontained = TRUE, background = "white")
      )
      ####################################################
      g <- as.directed(net_base)
      sp_g <- cluster_edge_betweenness(g)
      V(g)$community <- sp_g$membership
      V(g)$name <- V(g)$id
      V(g)$label <- V(g)$name
      
      GRN <- toVisNetworkData(net_base)
      GRN$nodes$name<- GRN$nodes$label
      GRN$nodes$color <- NULL
      
      GRN$nodes$label <- ifelse(GRN$nodes$fillcolor == "lavender", paste0(GRN$nodes$label,"up"),paste0(GRN$nodes$label,"down"))
      GRN$nodes$shape <- ifelse(GRN$nodes$label %in% c(paste0(violin_gene,"up"),paste0(violin_gene,"down")), "star", "circle")
      GRN$edges$id <- NULL
      
      gg <- g
      V(gg)$name <- V(gg)$label
      sub_graphs <- c()
      sub_objects <- c()
      unique_communities <- unique(V(gg)$community)
      temp_vertex <- c()
      print("Calculating centralities per community...")
      for (ie in 1: length(unique_communities)){
        OV <- which(V(gg)$community == unique_communities[ie])
        g_subgraph_temp <- induced_subgraph(gg, OV)
        temp_vertex0 <- which.max(igraph::betweenness(g_subgraph_temp))
        print("Best centrality:")
        print(temp_vertex0)
        temp_vertex <- c(temp_vertex,names(temp_vertex0))
        
      }
      signature <- paste(temp_vertex, collapse = '_')
      
      top_bc_network <-  c(top_bc_network, paste0(names(tfList)[i]))
      top_bc_vertex <- c(top_bc_vertex,signature)
      
      
      GRN$nodes$group <- V(g)$community
  
      mainP = paste0("Communities of optimized network for:: ",names(tfList)[i],"::" , violin_gene)
      print(paste0("Saving visualized graph for ", print(names(tfList)[i])))
      sp_g_5 <- visNetwork(GRN$nodes,GRN$edges, main = mainP, height = "700px", width = "100%") %>%
        visEdges(labelHighlightBold= "TRUE",arrows = "to") %>%
        visInteraction(zoomView = TRUE) %>%
        visOptions(highlightNearest = TRUE,nodesIdSelection = TRUE,selectedBy = "group") %>%
        visPhysics(stabilization = FALSE)  %>% visIgraphLayout(layout = "layout_nicely") %>%
        visEdges(smooth = FALSE)
      
      
      nodes <- data.frame(id =  V(g)$name,
                          label = V(g)$label,
                          shape = ifelse(V(g)$label %in% c(paste0(violin_gene,"up"),paste0(violin_gene,"down")), "star", "circle"))
      
      
      print("Saving Community detection graph...")
      try(
        visSave(sp_g_5, file = paste0(carnival_path, "/", names(tfList)[i],"/CM",".html"), selfcontained = TRUE, background = "white")
      )
      print("Done")
      
    }
  }
  print("Creating the labels for GNN...")
  
  # labels_csv$labels <- sapply(strsplit(labels_csv$mutation, split='_', fixed=TRUE),function(x) paste0(tail(x,1)))
  
  labels_csv$labels <- sapply(strsplit(labels_csv$mutation, split='_', fixed=TRUE),function(x) paste0(x[1],"_",x[2]))
  labels_csv <- labels_csv %>% dplyr::mutate(labels = if_else(str_detect(labels,"WT"),str_sub(labels, start = 1, end = 2) , labels_csv$labels))
  #labels_csv  <- labels_csv %>% group_by(mutation) %>% mutate(label = cur_group_id())
  #$label <- as.numeric( labels_csv$label) - 1
  #labels_csv$labels <- as.character( labels_csv$labels)
  labels_csv$mutation <- NULL
  write.csv(labels_csv,paste0(carnival_path, "/","GNN_labels.csv"))
  
  print("Printing upSetR graph...")
  no_mutations <- length(mutations)
  print(no_mutations)
  print(mutations)
  
  all_bc <- cbind(top_bc_network,top_bc_vertex)
  write.csv(all_bc,paste0(carnival_path,"/Betweenness_vertex_",disease_filename_j,".csv"))
  signatures <- list()
  index = 1;
  for (i in 1:no_mutations) {
    
    print(mutations[i])
    print("-------")
    keyword <- mutations[i]
    df <- as.data.frame(all_bc)
    df <- df %>% dplyr::filter(stringr::str_detect(top_bc_network,keyword))
    signature<- unique(paste(df$top_bc_vertex, collapse = "_"))
    temp <- as.list(strsplit(signature,"_")[[1]])
    signature <- unique(unlist(temp))
    print(signature)
    if (!is.null(signature)) {
      signatures[[index]] <- signature
      names(signatures)[index] <- keyword
      index = index + 1
    }
  }
  print("All signatures:")
  print(signatures)
  save(signatures, file="signatures.RData")
  UR <- upset(fromList(signatures), nsets = no_mutations,
              order.by = "freq",decreasing = T, empty.intersections = "on",sets.bar.color = "#56B4E9")
  print(UR)
  
  URD <- upset(fromList(signatures), nsets = no_mutations,
               order.by = "degree",decreasing = T, empty.intersections = "on",sets.bar.color = "#56B4E9")
  print(URD)
  
  
  pb <- progress_bar$new(total = total_Bar)
  print("Starting optimized network comparisons...")
  if (!file.exists(paste0(carnival_path,"//graph_heatmap_ID.csv"))) {
    for (i in 1:  length(tfList)) {
      pb$tick()
      Sys.sleep(1 / total_Bar)
      net_base <- all_networks[[i]]
      
      if (is.null(net_base) == FALSE) {
        
        for (jj in 1: length(tfList)) {
          score <- 0
          
          graph_heatmap[i,jj] <- 0
          graph_heatmap_ID[i,jj] <- 0
          graph_heatmap_DEL[i,jj] <- 0
          graph_heatmap_hotspot[i,jj]  <- 0
          graph_heatmap_KEY[i,jj]  <- 0
          
          if (i>jj) # only lower symmetric is enough
          {
            
            net_temp <- all_networks[[jj]]
            
            if (is.null(net_temp) == FALSE) {
              
              g_sim <- igraph::graph.intersection(net_base, net_temp, byname = "auto", keep.all.vertices = FALSE)
              
              if (gsize(net_base) > gsize(net_temp)) {
                score <- gsize(g_sim)/gsize(net_base)
                
              }
              else {
                
                score <- gsize(g_sim)/gsize(net_temp)
              }
              
              graph_heatmap[i,jj] <- score
              
              if (score >= top_score) {
                top_similar_1 <-  c(top_similar_1, paste0(names(tfList)[i]))
                top_similar_2 <-  c(top_similar_2, paste0(names(tfList)[jj]))
                top_similar_score <- c(top_similar_score,score)
                
              }
              
              graph_heatmap_ID[i,jj] <- ifelse(cloned$Variant_Classification[i] == cloned$Variant_Classification[jj],1,0)
              graph_heatmap_DEL[i,jj] <- ifelse(cloned$isDeleterious[i] == cloned$isDeleterious[jj],1,0)
              graph_heatmap_hotspot[i,jj] <- ifelse(cloned$isHotspot[i] == cloned$isHotspot[jj],1,0)
              if (key_opt_type == "Equal")
              {
                graph_heatmap_KEY[i,jj] <- ifelse(str_detect(names(tfList[i]),key_opt)==TRUE & str_detect(names(tfList[jj]),key_opt)==TRUE , 1,0)
              }
              else
              {
                graph_heatmap_KEY[i,jj] <- ifelse(str_detect(names(tfList[i]),key_opt)==TRUE & str_detect(names(tfList[jj]),key_opt)==FALSE |
                                                    str_detect(names(tfList[i]),key_opt)==FALSE & str_detect(names(tfList[jj]),key_opt)==TRUE, 1,0)
              }
            }
            else{score <- 0}
          }
        }
      }
      else{score <- 0} # if base network iterating is already null
    }
    print("Finished network comparisons...")
    
    colnames(graph_heatmap) <- names(tfList)
    row.names(graph_heatmap) <- names(tfList)
    graph_heatmap[is.na(graph_heatmap)] <- 0
    graph_heatmap_ID[is.na(graph_heatmap_ID)] <- 0
    graph_heatmap_DEL[is.na(graph_heatmap_DEL)] <- 0
    graph_heatmap_hotspot[is.na(graph_heatmap_hotspot)] <- 0
    graph_heatmap_KEY[is.na(graph_heatmap_KEY)] <- 0
    
    
    #print(graph_heatmap)
    setwd(results_dir)
    total_perc <- c() # similar
    total_perc_filtered <- c()  # similar and same type of mutation
    total_perc_filtered_del <- c() # similar and same flag for deleterious
    total_perc_filtered_del_S <- c() # similar and same type of mutation and deleterious flag
    total_perc_filtered_hotspot <- c() # similar and from same hotspot
    total_perc_filtered_KEY <- c() # user given
    
    write.csv(graph_heatmap,paste0(carnival_path,"/Networks_Similarity_Scores_",disease_filename_j,".csv"))
    print("Finishing matrices init...")
    # now transform heatamp matrix to only account for similarity scores that come from same mutation types
    # this is done through filtering with graph_heatmap_ID as follows:
    graph_heatmap_filtered <- graph_heatmap*as.vector(graph_heatmap_ID) # this zeroes target based on zeroes in heatmap
    graph_heatmap_filtered_del <- graph_heatmap*as.vector(graph_heatmap_DEL)
    graph_heatmap_filtered_hotspot <- graph_heatmap*as.vector(graph_heatmap_hotspot)
    graph_heatmap_filtered_KEY <- graph_heatmap*as.vector(graph_heatmap_KEY)
    
    
    colnames(graph_heatmap_filtered) <- names(tfList)
    row.names(graph_heatmap_filtered) <- names(tfList)
    
    colnames(graph_heatmap_filtered_del) <- names(tfList)
    row.names(graph_heatmap_filtered_del) <- names(tfList)
    
    colnames(graph_heatmap_filtered_hotspot) <- names(tfList)
    row.names(graph_heatmap_filtered_hotspot) <- names(tfList)
    
    colnames(graph_heatmap_filtered_KEY) <- names(tfList)
    row.names(graph_heatmap_filtered_KEY) <- names(tfList)
    
    write.csv(graph_heatmap_filtered, "graph_heatmap_filtered.csv")
    write.csv(graph_heatmap_filtered_del, "graph_heatmap_filtered_del.csv")
    write.csv(graph_heatmap_filtered_hotspot, "graph_heatmap_filtered_hotspot.csv")
    write.csv(graph_heatmap_filtered_KEY, "graph_heatmap_filtered_KEY.csv")
    
    
    total_no_nets <- nrow(which( graph_heatmap >0, arr.ind=TRUE))
    total_no_nets_ID <- nrow(which( graph_heatmap_filtered >0, arr.ind=TRUE))
    total_n_nets_DEL <- nrow(which( graph_heatmap_filtered_del >0, arr.ind=TRUE))
    total_n_nets_hotspot <- nrow(which( graph_heatmap_filtered_hotspot >0, arr.ind=TRUE))
    total_n_nets_KEY <- nrow(which( graph_heatmap_filtered_KEY >0, arr.ind=TRUE))
    
    
    print("Finishing loops...")
    for (i in 1:length(network_similarity_threshold)){
      temp_score <- sum(graph_heatmap >= network_similarity_threshold[i])
      temp_score <- round(temp_score/total_no_nets, digits = 2)
      total_perc <- c(total_perc,temp_score*100)
      
    }
    
    # now the same but only for same type of mutation:
    for (i in 1:length(network_similarity_threshold)){
      #temp_score <- sum(sapply(graph_heatmap, function(x) sum(x>network_similarity_threshold[i])))
      temp_score <- sum(graph_heatmap_filtered >= network_similarity_threshold[i])
      temp_score <- round(temp_score/total_no_nets_ID, digits = 2)
      total_perc_filtered <- c(total_perc_filtered,temp_score*100)
      
    }
    
    
    # now the same but only for same type of mutation:
    for (i in 1:length(network_similarity_threshold)){
      #temp_score <- sum(sapply(graph_heatmap, function(x) sum(x>network_similarity_threshold[i])))
      temp_score <- sum(graph_heatmap_filtered_del >= network_similarity_threshold[i])
      temp_score <- round(temp_score/total_n_nets_DEL, digits = 2)
      total_perc_filtered_del <- c(total_perc_filtered_del,temp_score*100)
      
    }
    
    # now the same but only for hotspots of TP53:
    for (i in 1:length(network_similarity_threshold)){
      #temp_score <- sum(sapply(graph_heatmap, function(x) sum(x>network_similarity_threshold[i])))
      temp_score <- sum(graph_heatmap_filtered_hotspot >= network_similarity_threshold[i])
      temp_score <- round(temp_score/total_n_nets_hotspot, digits = 2)
      total_perc_filtered_hotspot <- c(total_perc_filtered_hotspot,temp_score*100)
      
    }
    # now for the user key
    for (i in 1:length(network_similarity_threshold)){
      #temp_score <- sum(sapply(graph_heatmap, function(x) sum(x>network_similarity_threshold[i])))
      temp_score <- sum(graph_heatmap_filtered_KEY >= network_similarity_threshold[i])
      temp_score <- round(temp_score/total_n_nets_KEY, digits = 2)
      total_perc_filtered_KEY <- c(total_perc_filtered_KEY,temp_score*100)
      
    }
    
    print("Writing files...")
    
    top_all <- cbind(top_similar_1,top_similar_2,top_similar_score)
    write.csv(top_all,paste0(carnival_path,"/top_similar_nets.csv"))
    write.csv(graph_heatmap_ID,paste0(carnival_path,"/graph_heatmap_ID.csv"))
    write.csv(graph_heatmap_DEL,paste0(carnival_path,"/graph_heatmap_DEL.csv"))
    write.csv(graph_heatmap_hotspot,paste0(carnival_path,"/graph_heatmap_hotspot.csv"))
    print("Finished..")
    write.csv(graph_heatmap_KEY,paste0(carnival_path,paste0("/graph_heatmap_",key_opt,".csv")))
    
    
    ylimit <- length(tfList)
    xlimit <- ylimit
    if (xlimit > 100){
      xlimit <- 100
      ylimit <- xlimit
    }
    print("Printing corrplots...")
    
    title2 = paste0("Percentage of networks with similarity greater or equal than ",
                    paste(network_similarity_threshold*100, collapse = ', '), "%"," : ",
                    paste(paste0(total_perc,"%"), collapse = ', '), " in ", disease_filename_j)
    sp_CARNIVAL_1 <- corrplot(graph_heatmap[1:xlimit,1:ylimit], type = 'lower', order = 'hclust', tl.col = 'black',
                              cl.ratio = 0.2, tl.srt = 45, col = COL2('PuOr', 10),title = title2,mar=c(0,0,1,0),is.corr = FALSE,
                              col.lim = c(0,1))
    
    #######################################################
    
    title2 = paste0("Percentage of networks with similarity greater or equal than ",
                    paste(network_similarity_threshold*100, collapse = ', '), "% and same type of mutation"," : ",
                    paste(paste0(total_perc_filtered,"%"), collapse = ', '), " in ", disease_filename_j)
    sp_CARNIVAL_2 <-  corrplot(graph_heatmap_filtered[1:xlimit,1:ylimit], type = 'lower', order = 'hclust', tl.col = 'black',
                               cl.ratio = 0.2, tl.srt = 45, col = COL2('PuOr', 10),title = title2,mar=c(0,0,1,0),is.corr = FALSE,
                               col.lim = c(0,1))
    
    #######################################################
    
    title2 = paste0("Percentage of networks with similarity greater or equal than ",
                    paste(network_similarity_threshold*100, collapse = ', '), "% and same deleterious flag"," : ",
                    paste(paste0(total_perc_filtered_del,"%"), collapse = ', '), " in ", disease_filename_j)
    
    sp_CARNIVAL_3 <-  corrplot(graph_heatmap_filtered_del[1:xlimit,1:ylimit], type = 'lower', order = 'hclust', tl.col = 'black',
                               cl.ratio = 0.2, tl.srt = 45, col = COL2('PuOr', 10),title = title2,mar=c(0,0,1,0),is.corr = FALSE,
                               col.lim = c(0,1))
    #######################################################
    title2 = paste0("Percentage of networks with similarity greater or equal than ",
                    paste(network_similarity_threshold*100, collapse = ', '), "% and same hotspot"," : ",
                    paste(paste0(total_perc_filtered_hotspot,"%"), collapse = ', '), " in ", disease_filename_j)
    sp_CARNIVAL_4 <- corrplot(graph_heatmap_filtered_hotspot[1:xlimit,1:ylimit], type = 'lower', order = 'hclust', tl.col = 'black',
                              cl.ratio = 0.2, tl.srt = 45, col = COL2('PuOr', 10),title = title2,mar=c(0,0,1,0),is.corr = FALSE,
                              col.lim = c(0,1))
    
    #######################################################
    
    title2 = paste0("Percentage of networks with similarity greater or equal than ",
                    paste(network_similarity_threshold*100, collapse = ', '), paste("% and", key_opt,"", key_opt_type)," : ",
                    paste(paste0(total_perc_filtered_KEY,"%"), collapse = ', '), " in ", disease_filename_j)
    sp_CARNIVAL_5 <- corrplot(graph_heatmap_filtered_KEY[1:xlimit,1:ylimit], type = 'lower', order = 'hclust', tl.col = 'black',
                              cl.ratio = 0.2, tl.srt = 45, col = COL2('PuOr', 10),title = title2,mar=c(0,0,1,0),is.corr = FALSE, #this is the parameter you were looking for
                              col.lim = c(0,1),tl.cex = 1)
    
    print(sp_CARNIVAL_1)
    print(sp_CARNIVAL_2)
    print(sp_CARNIVAL_3)
    print(sp_CARNIVAL_4)
    print(sp_CARNIVAL_4)
    
    sp_CARNIVAL <- c(sp_CARNIVAL_1,sp_CARNIVAL_2,sp_CARNIVAL_3,sp_CARNIVAL_4,sp_CARNIVAL_5)
    
    print("Assigning radar plot data...")
    k = 0
    radar_plot_data[iterator_index,1] <- disease_filename_j
    for (i in 1:length(network_similarity_threshold)){
      
      radar_plot_data[iterator_index,k+1] <- total_perc[i]
      k = k + 1;
    }
    
    for (i in 1:length(network_similarity_threshold)){
      
      radar_plot_data[iterator_index,k+1] <- total_perc_filtered[i]
      k = k + 1;
    }
    
    for (i in 1:length(network_similarity_threshold)){
      
      radar_plot_data[iterator_index,k+1] <- total_perc_filtered_del[i]
      k = k + 1;
    }
    
    for (i in 1:length(network_similarity_threshold)){
      
      radar_plot_data[iterator_index,k+1] <- total_perc_filtered_hotspot[i]
      k = k + 1;
    }
    
    for (i in 1:length(network_similarity_threshold)){
      
      radar_plot_data[iterator_index,k+1] <- total_perc_filtered_KEY[i]
      k = k + 1;
    }
    
    
    
    print(radar_plot_data)
    write.csv(cloned,paste0(disease_filename_j,".csv"))
    write.csv(radar_plot_data,"radar_plot_data_temp.csv")
    return_list<- list("radar_plot_data" = radar_plot_data, sp_CARNIVAL,total_perc,total_perc_filtered,total_perc_filtered_del,total_perc_filtered_hotspot)
  }
  else{print("Networks have been compared already - file 'graph_heatmap_ID.csv' already exists! Delete if you wish to rerun the comparisons")}
  
  #### RUN PYTHON SCRIPT TO CONVERT .dot to .graphml so that igraph can read it
  setwd(source_code_dir)
  if (GNN_flag == "TRUE") {
    print("Running now Python script to classify optimized networks at graph level using Graph Neural Networks...")
    if (.Platform$OS.type == "unix") {
      command <- paste0("python3 GNN3.py ",'"',carnival_path,'"'," ",'"',paste0(carnival_path, "//","GNN_labels.csv"),'"'," ", epochs," ",lr," ",tts," ",update_it, " ", bsu, " ", hidden_dim, " ", dropout )
    }
    else if  (.Platform$OS.type == "windows") {
      command <- paste0("python GNN3.py ",'"',carnival_path,'"'," ",'"',paste0(carnival_path, "//","GNN_labels.csv"),'"'," ",epochs," ",lr," ",tts," ",update_it, " ", bsu, " ", hidden_dim, " ", dropout)
    }
    
    
    
    print(command)
    system(command)
    
    
    Sys.sleep(10)
    
  }
  #stop()
  print("Done!")
  
  return()
  
  
  #############################################################################
  
  
}
