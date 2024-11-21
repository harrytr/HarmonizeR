runPROGENy <- function(df,weight_matrix,k = 10000, z_scores = T, get_nulldist = F)
{
  # library(tidyverse)
  library(dplyr)
  resList <- list()
  if(get_nulldist)
  {
    nullDist_list <- list()
  }
  
  for(i in 2:length(df[1,]))
  {
    current_df <- df[,c(1,i)]
    current_df <- current_df[complete.cases(current_df),]
    t_values <- current_df[,2]
    
    current_weights <- weight_matrix
    
    names(current_df)[1] <- "ID"
    names(current_weights)[1] <- "ID"
    
    common_ids <- merge(current_df, current_weights, by = "ID")
    common_ids <- common_ids$ID
    common_ids <- as.character(common_ids)
    
    row.names(current_df) <- current_df$ID
    current_df <- as.data.frame(current_df[common_ids,-1])
    
    row.names(current_weights) <- current_weights$ID
    current_weights <- as.data.frame(current_weights[common_ids,-1])
    
    current_mat <- as.matrix(current_df)
    current_weights <- t(current_weights)
    
    scores <- as.data.frame(current_weights %*% current_mat)
    
    null_dist_t <- replicate(k, sample(t_values,length(current_mat[,1]), replace = F))
    
    null_dist_scores <- current_weights %*% null_dist_t
    
    if(get_nulldist)
    {
      nullDist_list[[i-1]] <- null_dist_scores
    }
    
    if(z_scores)
    {
      scores$mean <- apply(null_dist_scores,1,mean)
      scores$sd <- apply(null_dist_scores,1,sd)
      resListCurrent <- (scores[,1]-scores[,2])/scores[,3]
      names(resListCurrent) <- names(weight_matrix[,-1])
      resList[[i-1]] <- resListCurrent
    }
    else
    {
      for(j in 1:length(weight_matrix[,-1]))
      {
        ecdf_function <- ecdf(null_dist_scores[j,])
        scores[j,1] <- ecdf_function(scores[j,1])
      }
      score_probas <- scores*2-1
      
      resListCurrent <- score_probas[,1]
      names(resListCurrent) <- names(weight_matrix[,-1])
      resList[[i-1]] <- resListCurrent
    }
  }
  names(resList) <- names(df[,-1])
  resDf <- as.data.frame(resList)
  colnames(resDf)=colnames(df)[-1]
  if(get_nulldist)
  {
    names(nullDist_list) <- names(df[,-1])
    return(list(resDf, nullDist_list))
  }
  else
  {
    return(resDf)
  }
  
}