library(tidyverse)
library(synapser)
synLogin()

template <- synGet('syn21321426')$path %>% read_csv

##create null model that considers all targets in template

## FULL TEMPLATE null model

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

write_rds(null_model_sc1, 'analysis/null_model_sc1_full_template.rds')
write_rds(null_model_sc2, 'analysis/null_model_sc2_full_template.rds')

synStore(File('analysis/null_model_sc1_full_template.rds', parentId = "syn22332108"))
synStore(File('analysis/null_model_sc2_full_template.rds', parentId = "syn22332108"))


# 
# 
# ##chembl target universe
# targs <- chembl$target %>% unique()
# sum(targs %in% template$target, na.rm = TRUE)
# 
# null_model_sc1_subset_chembl <- lapply(1:1000, function(x){
# 
#   foo <- template %>%
#     gather(cmpd_id, confidence ,-target) %>%
#     group_by(cmpd_id) %>%
#     sample_n(10, replace = F) %>%
#     arrange(-confidence, target) %>%
#     select(-confidence) %>%
#     set_names(c('target_random', 'cmpd_id')) %>%
#     filter(target_random %in% targs) %>%
#     nest()
# 
# })
# 
# null_model_sc2_subset_chembl <- lapply(1:1000, function(x){
# 
#   foo <- template %>%
#     gather(cmpd_id, confidence ,-target) %>%
#     group_by(cmpd_id) %>%
#     arrange(-confidence, target) %>%
#     sample_frac(1, replace = F) %>%
#     select(-confidence) %>%
#     set_names(c('target_random', 'cmpd_id')) %>%
#     filter(target_random %in% targs) %>%
#     nest()
# 
# })
# 
# write_rds(null_model_sc1_subset_chembl, 'analysis/null_model_sc1_subset_chembl.rds')
# write_rds(null_model_sc2_subset_chembl, 'analysis/null_model_sc2_subset_chembl.rds')
# 
# synStore(File('analysis/null_model_sc1_subset_chembl.rds', parentId = "syn22332108"))
# synStore(File('analysis/null_model_sc2_subset_chembl.rds', parentId = "syn22332108"))
# 
# file.remove('analysis/null_model_sc1_subset_chembl.rds')
# file.remove('analysis/null_model_sc2_subset_chembl.rds')


# ###drugbank
# 
# targs <- drugbank$target %>% unique()
# sum(targs %in% template$target, na.rm = TRUE)
# 
# null_model_sc1_subset_drugbank <- lapply(1:1000, function(x){
# 
#   foo <- template %>%
#     gather(cmpd_id, confidence ,-target) %>%
#     group_by(cmpd_id) %>%
#     sample_n(10, replace = F) %>%
#     arrange(-confidence, target) %>%
#     select(-confidence) %>%
#     set_names(c('target_random', 'cmpd_id')) %>%
#     filter(target_random %in% targs) %>%
#     nest()
# 
# })
# 
# null_model_sc2_subset_drugbank <- lapply(1:1000, function(x){
# 
#   foo <- template %>%
#     gather(cmpd_id, confidence ,-target) %>%
#     group_by(cmpd_id) %>%
#     arrange(-confidence, target) %>%
#     sample_frac(1, replace = F) %>%
#     select(-confidence) %>%
#     set_names(c('target_random', 'cmpd_id')) %>%
#     filter(target_random %in% targs) %>%
#     nest()
# 
# })
# 
# write_rds(null_model_sc1_subset_drugbank, 'analysis/null_model_sc1_subset_drugbank.rds')
# write_rds(null_model_sc2_subset_drugbank, 'analysis/null_model_sc2_subset_drugbank.rds')
# 
# synStore(File('analysis/null_model_sc1_subset_drugbank.rds', parentId = "syn22332108"))
# synStore(File('analysis/null_model_sc2_subset_drugbank.rds', parentId = "syn22332108"))
# 
# file.remove('analysis/null_model_sc1_subset_drugbank.rds')
# file.remove('analysis/null_model_sc2_subset_drugbank.rds')
# 
# ##drugtargetcommons
# 
# targs <- dtc$target %>% unique()
# sum(targs %in% template$target, na.rm = TRUE)
# 
# null_model_sc1_subset_dtc <- lapply(1:1000, function(x){
# 
#   foo <- template %>%
#     gather(cmpd_id, confidence ,-target) %>%
#     group_by(cmpd_id) %>%
#     sample_n(10, replace = F) %>%
#     arrange(-confidence, target) %>%
#     select(-confidence) %>%
#     set_names(c('target_random', 'cmpd_id')) %>%
#     filter(target_random %in% targs) %>%
#     nest()
# 
# })
# 
# null_model_sc2_subset_dtc <- lapply(1:1000, function(x){
# 
#   foo <- template %>%
#     gather(cmpd_id, confidence ,-target) %>%
#     group_by(cmpd_id) %>%
#     arrange(-confidence, target) %>%
#     sample_frac(1, replace = F) %>%
#     select(-confidence) %>%
#     set_names(c('target_random', 'cmpd_id')) %>%
#     filter(target_random %in% targs) %>%
#     nest()
# 
# })
# 
# write_rds(null_model_sc1_subset_dtc, 'analysis/null_model_sc1_subset_dtc.rds')
# write_rds(null_model_sc2_subset_dtc, 'analysis/null_model_sc2_subset_dtc.rds')
# 
# synStore(File('analysis/null_model_sc1_subset_dtc.rds', parentId = "syn22332108"))
# synStore(File('analysis/null_model_sc2_subset_dtc.rds', parentId = "syn22332108"))
# 
# file.remove('analysis/null_model_sc1_subset_dtc.rds')
# file.remove('analysis/null_model_sc1_subset_dtc.rds')
