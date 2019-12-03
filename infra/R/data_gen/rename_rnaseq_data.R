library(tidyverse)
library(synapser)
# library(random)
synLogin()

rnaseq_zip <- synGet("syn21036601")$path
rnaseq_files <- unzip(rnaseq_zip, list = T)$Name
rnaseq_paths <- unzip(rnaseq_zip, files = rnaseq_files[2:12])

rnaseq <- purrr::map(rnaseq_paths, readr::read_csv)

## code to generate random IDs. Do not rerun unless ids need to be recreated. this uses the random.org api to get random 2-letter values for each drug
# cmpd_names <- synGet("syn21036464")$path %>% readr::read_csv() %>%
#   rename(cmpd = X1) %>%
#   mutate(cmpd_id = random::randomStrings(n = nrow(.), len = 2, upperalpha = T, loweralpha = F, digits = F)[,1]) %>%
#   mutate(cmpd_id = glue::glue('cmpd_{cmpd_id}')) %>%
#   select(cmpd_id, cmpd) %>% 
#   write_csv('random_cmpd_map.csv')
# 
# synStore(File("random_cmpd_map.csv", parentId = "syn21036376"))

cmpd_ids <- synGet("syn21197825")$path %>% read_csv

id_map <- cmpd_ids %>% 
  add_row(cmpd_id = "cmpd_dmso", cmpd =  "DMSO") %>% 
  add_row(cmpd_id = "cmpd_untreated", cmpd = "UNTREATED")

id_map_vec <- id_map$cmpd_id
names(id_map_vec) <- id_map$cmpd

rnaseq_anon <- lapply(rnaseq, function(x){ #replace drug names with randomized ids from map above
  x <- x %>% 
    rename(hgnc_symbol = X1) 
  
  ##append _0 if no replicate id was added on read in
  colnames(x)[is.na(stringr::str_match(colnames(x), "_\\d+$"))][-1] <- 
    paste0(colnames(x)[is.na(stringr::str_match(colnames(x), "_\\d+$"))][-1], "_0")
  
  ##map to new cmpd names
  colnames(x) <- stringr::str_replace_all(colnames(x), id_map_vec)
  
  x
})

paths <- sapply(rnaseq_paths, function(x){ 
  
  path<-str_extract(x, "[A-Za-z0-9-]+.csv") %>% print()
  
  write_csv(rnaseq_anon[[which(rnaseq_paths == x)]], path)
  
  path
  
})

zipped_path<-zip('rnaseq_concealed.zip',paths)

synStore(File("rnaseq_concealed.zip", parentId = "syn21036376"))




