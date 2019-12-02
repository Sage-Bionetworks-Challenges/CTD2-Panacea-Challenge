library(tidyverse)

template <- read_csv('template.csv')

null_model <- lapply(1:1000, function(x){
  
 foo <- template %>% 
    gather(cmpd_id, confidence ,-target) %>% 
    group_by(cmpd_id) %>% 
    arrange(-confidence, target) %>% 
    sample_n(10, replace = F) %>% 
   select(-confidence) %>%
   set_names(c('target_random', 'cmpd_id')) %>% 
   nest()
  
})

write_rds(null_model, 'null_model.rds')
