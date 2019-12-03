library(tidyverse)

template <- read_csv('template.csv')

null_model_sc1 <- lapply(1:1000, function(x){
  
  foo <- template %>% 
    gather(cmpd_id, confidence ,-target) %>% 
    group_by(cmpd_id) %>% 
    sample_n(10, replace = F) %>% 
    arrange(-confidence, target) %>% 
    select(-confidence) %>%
    set_names(c('target_random', 'cmpd_id')) %>% 
    nest()
  
})

null_model_sc2 <- lapply(1:1000, function(x){
  
  foo <- template %>% 
    gather(cmpd_id, confidence ,-target) %>% 
    group_by(cmpd_id) %>% 
    arrange(-confidence, target) %>% 
    sample_frac(1, replace = F) %>% 
    select(-confidence) %>%
    set_names(c('target_random', 'cmpd_id')) %>% 
    nest()
  
})

write_rds(null_model_sc1, 'null_model_sc1.rds')
write_rds(null_model_sc2, 'null_model_sc2.rds')

write_rds(null_model_sc1[1:10], 'null_model_sc1_leaderboard.rds')
write_rds(null_model_sc2[1:10], 'null_model_sc2_leaderboard.rds')
