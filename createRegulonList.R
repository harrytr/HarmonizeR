#'\code{createRegulonList}
#'
#' Creating regulon list similar to vipper from an interaction data-frame
#' 
#' @param regulon_table the regulon data-frame with columns the source, sign and
#' target of the TF interaction.
#'
#'Enio Gjerga, 2020

createRegulonList <- function(regulon_table = regulon_table){
  
  regulon = list()
  sourceList = unique(regulon_table[, 1])
  
  for(ii in 1:length(sourceList)){
    
    currRegulon <- list()
    
    currRegulon[[length(currRegulon)+1]] <- 
      regulon_table[, 2][which(regulon_table[, 1]==sourceList[ii])]
    
    names(currRegulon[[1]]) <- 
      regulon_table[, 3][which(regulon_table[, 1]==sourceList[ii])]
    
    currRegulon[[length(currRegulon)+1]] <- 
      rep(1, length(which(regulon_table[, 1]==sourceList[ii])))
    
    names(currRegulon) <- c("tfmode", "likelihood")
    
    regulon[[length(regulon)+1]] <- currRegulon
    
  }
  
  names(regulon) <- sourceList
  
  return(regulon)
  
}