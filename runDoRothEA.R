runDoRothEA<-function(df, regulon, confidence_level=c('A','B','C'), write2file = NULL){
  # library(tidyverse)
  library(dplyr)
  library(purrr)
  library(viper)
  library(tibble)
  library(tidyr)
  # names(regulon) <- sapply(strsplit(names(viper_regulon), split = ' - '), head, 1)
  names(regulon) <- sapply(strsplit(names(regulon), split = ' - '), head, 1)
  filtered_regulon <- regulon %>%
    map_df(.f = function(i) {
      tf_target = i$tfmode %>%
        enframe(name = "target", value="mor") %>%
        mutate(likelihood = i$likelihood)
    },.id = "tf")  %>%
    separate(tf, into=c("tf", "conf"), sep="_") %>%
    filter(conf %in% confidence_level) %>%
    arrange(tf)%>%
    split(.$tf) %>%
    map(function(dat) {
      tf = dat %>% distinct(tf) %>% pull()
      targets = setNames(dat$mor, dat$target)
      likelihood = dat$likelihood
      list(tfmode =targets, likelihood = likelihood)})
  
  TF_activities = as.data.frame(viper::viper(eset = df, regulon = filtered_regulon, nes = T, method = 'none', minsize = 4, eset.filter = F))
  if(!is.null(write2file)){write.csv2(TF_activities, file = write2file)}
  
  return(TF_activities)
}