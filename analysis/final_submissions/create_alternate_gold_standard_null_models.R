library(tidyverse)

template <- read_csv('template.csv')

## FULL TEMPLATE null model 

# null_model_sc1 <- lapply(1:1000, function(x){
#   
#   foo <- template %>% 
#     gather(cmpd_id, confidence ,-target) %>% 
#     group_by(cmpd_id) %>% 
#     sample_n(10, replace = F) %>% 
#     arrange(-confidence, target) %>% 
#     select(-confidence) %>%
#     set_names(c('target_random', 'cmpd_id')) %>% 
#     nest()
#   
# })
# 
# null_model_sc2 <- lapply(1:1000, function(x){
#   
#   foo <- template %>% 
#     gather(cmpd_id, confidence ,-target) %>% 
#     group_by(cmpd_id) %>% 
#     arrange(-confidence, target) %>% 
#     sample_frac(1, replace = F) %>% 
#     select(-confidence) %>%
#     set_names(c('target_random', 'cmpd_id')) %>% 
#     nest()
#   
# })
# 
# write_rds(null_model_sc1, 'null_model_sc1.rds')
# write_rds(null_model_sc2, 'null_model_sc2.rds')
# 
# write_rds(null_model_sc1[1:100], 'null_model_sc1_leaderboard.rds')
# write_rds(null_model_sc2[1:100], 'null_model_sc2_leaderboard.rds')

##gold standard only null model
##This is what will be the final null model for the challenge

targs <- c("EGFR","CSNK2A2","BMP2K","AAK1","Q6ZSR9","GAK","PRKD2","PRKD3","SIK2","SIK3","PLK4","AURKA","AURKB","PTK2",
           "FER","PTK2B","TNK1","TNK2","MAP2K5","MAP2K1","MAP2K2","PDGFRB","STK10","SLK","MAP4K2","MAP4K5","MAP4K3",
           "MET","CDK2","CDK5","CDK16","CDK9","GSK3A","GSK3B","CDK7","ERCC2","CDK1","CDK12","CDK13","CLK1","DYRK1A",
           "CDK17","CDK4","CDK6","FGFR1","MAPK10","MAPK8","MAPK9","ADCK1","STK16","NEK9","FECH","CSNK1A1","CSNK1E",
           "CSNK1D","MAP3K4","WEE1","MST1R","PIP4K2C","IRAK1","TAOK3","PIM1","TAOK1","CIT","PRKX","PRKY","CDC42BPA",
           "CDC42BPB","PRKACA","PRKACB","PRKCQ","CLK4","PAK6","PRKACG","RPS6KB1","MAPK7","PRKG1","PRKCI","CDKL5","CDK18",
           "ADRBK1","STK11","EPHB3","BRAF","NLK","TGFBR2","CSK","MYLK","FGR","TESK1","KIT","ULK1","FES","IGF1R","INSR",
           "JAK2","STK26","CAMK4","TAOK2","EIF2AK1","ICK","CSNK1G2","CSNK1G1","CSNK1G3","CLK2","TP53RK","MINK1","MAPK15",
           "ERN1","MAP2K6","MAP3K5","MAP3K6","MAP2K3","ERN2","MAPK1","MAPK3","AKT3","AKT1","AKT2","RPS6KA5","ACVRL1",
           "MAPKAPK5","NEK1","ARAF","ACAD10","ADCK3","ACVR2B","PAK2","CLK3","MAPKAPK3","PLK1","PDXK","DCK","ADK","PKMYT1",
           "NUAK2","HSP90AB2P","STK24","CMPK1","GARS","ACOX3","ACAD11","DCTPP1","NEK7","SMC2","TOP2B","MYH10","PIM2",
           "CAMK1G","PIP4K2A","PIP4K2B","AK2","NME2","CDC42BPG","AP1G1","DHCR24","SLC25A5","ACTR3","SNRNP200",
           "CARS","DNAJA1","PFKP","MAT2A","IPO7","STK38L","AIMP1","CDC7","PRKAR2A","EPHA1","EPHA7","BMPR1B","BMPR1A",
           "TYK2","BMPR2","MAP3K2","MAP3K3","SYK","NEK3","BUB1","LATS1","MAP3K11","IKBKE","TBK1","STK4","STK3","CDK3",
           "MARK4","MELK","IRAK4","RPS6KA4","RPS6KA1","RPS6KA3","MARK2","MARK3","RPS6KA6","PRKAG2","PRKAA1","PRKAG1","PAK4",
           "CHEK1","NTRK1","CAMK2D","CAMK2G","CAMKK2","PHKG2","MAP4K4","TNIK","MYLK3","MAP4K1","JAK1","ULK3","PRKCD",
           "PRKCA","PRKCB","PKN1","PKN2","ROCK2","ROCK1","NQO2","ACVR1","IRAK3","FLT3","RET","DDR1","DDR2","ABL2","BCR",
           "ABL1","RIPK2","MAPKAPK2","MAPK11","MAPK14","ZAK","RIPK3","YES1","LCK","SRC","FYN","HCK","LYN","FRK","EPHA2",
           "EPHA5","EPHB2","EPHB4","EPHA4","BTK","TEC","MAP3K1","LIMK1","LIMK2","PTK6","EPHB6","ACVR1B","TGFBR1")

null_model_sc1_subset <- lapply(1:1000, function(x){
  
  foo <- template %>% 
    gather(cmpd_id, confidence ,-target) %>% 
    group_by(cmpd_id) %>% 
    sample_n(10, replace = F) %>% 
    arrange(-confidence, target) %>% 
    select(-confidence) %>%
    set_names(c('target_random', 'cmpd_id')) %>% 
    filter(target_random %in% targs) %>% 
    nest()
  
})

null_model_sc2_subset <- lapply(1:1000, function(x){
  
  foo <- template %>% 
    gather(cmpd_id, confidence ,-target) %>% 
    group_by(cmpd_id) %>% 
    arrange(-confidence, target) %>% 
    sample_frac(1, replace = F) %>% 
    select(-confidence) %>%
    set_names(c('target_random', 'cmpd_id')) %>% 
    filter(target_random %in% targs) %>% 
    nest()
  
})

write_rds(null_model_sc1_subset, 'null_model_sc1_subset.rds')
write_rds(null_model_sc2_subset, 'null_model_sc2_subset.rds')

write_rds(null_model_sc1_subset[1:100], 'null_model_sc1_leaderboard_subset.rds')
write_rds(null_model_sc2_subset[1:100], 'null_model_sc2_leaderboard_subset.rds')




