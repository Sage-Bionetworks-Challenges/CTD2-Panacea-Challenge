library(tidyverse)
library(synapser)
library(random)
synLogin()


dr_zip <- synGet("syn21036458")$path
dr_files <- unzip(dr_zip, list = T)$Name
dr_files <- dr_files[grep("dose-responses/[A-Za-z0-9-]+.csv",dr_files)]

dr_paths <- unzip(dr_zip, files = dr_files)

dr <- purrr::map(dr_paths, readr::read_csv)

cmpd_ids <- synGet("syn21197825")$path %>% read_csv

id_map <- cmpd_ids %>% 
  add_row(cmpd_id = "cmpd_dmso", cmpd =  "DMSO") %>% 
  add_row(cmpd_id = "cmpd_untreated", cmpd = "UNTREATED")

id_map_vec <- id_map$cmpd_id
names(id_map_vec) <- id_map$cmpd


dr_anon <- lapply(dr, function(x){
  x <- x %>% 
    rename(dose_log10_uM = X1)
  
  colnames(x) <- stringr::str_replace_all(colnames(x), id_map_vec)
  
  x
})

paths <- sapply(dr_paths, function(x){
  
  path<-str_extract(x, "[A-Za-z0-9-]+.csv") %>% print()
  
  write_csv(dr_anon[[which(dr_paths == x)]], path)
  
  path
  
})

zipped_path<-zip('dose_response_concealed.zip',paths)

synStore(File("dose_response_concealed.zip", parentId = "syn21036376"))




