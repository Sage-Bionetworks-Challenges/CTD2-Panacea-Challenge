library(dplyr)
library(magrittr)
library(readr)

###made up data until i get the real stuff
pred <- dplyr::tribble(
  ~compound_id, ~hgnc_id, ~confidence,
  "drugA", 'AURdfKA', 0.9,
  "drugA", 'AURKB', 0.1
)

gold <- dplyr::tribble(
  ~compound_id, ~target, ~rank,
  "drugA", 'AURKA', 10,
  "drugB", 'AURKB', 9
)

validate <- function(prediction_path, gold_path, hgnc_ids){
  
  pred <- readr::read_csv(prediction_path)
  gold <- readr::read_csv(gold_path)
  
  ###configure validation
  ncol_req <- 3
  colnames_req <-  c('compound_id', 'hgnc_id', 'confidence')
  compound_ids <- unique(gold$compound_id)
  
  ###
  hgnc_ids <- readr::read_tsv("hgnc_ids.txt") %>%  
    dplyr::filter(Status == "Approved") %>% 
    dplyr::pull('Approved symbol') 
  
  errs <- list()
  
  trim_vec <- function(vec, trim = 10){
    if(length(vec) > trim){
      vec <- vec[1:trim]
      vec <- as.character(vec)
      vec[11] <- '...'
    }else{
      vec
    }
  }

  if(ncol(pred)<ncol_req){
    errs["ncol_short"] <- paste0("Prediction file is missing rows. Only ", ncol(pred), " rows detected.")
  }
  
  if(ncol(pred)<ncol_req){
    errs["ncol_long"] <- paste0("Prediction file is missing rows. ", ncol(pred), " rows detected.")
  }
  
  if(isTRUE(colnames(pred) %in% colnames_req)){
    errs["colnames"] <- paste0("Column names are not correct. Column names must be ", cat(colnames_req))
  }
  
  if(!identical(pred$id, gold$id)){
    errs["id_mismatch"] <- paste0("Prediction file id column does not have matching values to gold standard file.")
  }
  
  if(!is.numeric(pred$confidence)){
    errs["non_numeric"] <- paste0("Predictions are not all numeric values.")
  }
  
  if(any(!pred$hgnc_id %in% hgnc_ids)){
    invalid <- unique(pred$hgnc_id[!pred$hgnc_id %in% hgnc_ids]) %>% trim_vec()
    errs["non_hgnc"] <- paste0("Invalid HGNC identifiers were included in your prediction file (up to 10 displayed): ", invalid)
  }
  
  if(all(gold$compound_id %in% pred$compound_id)){
    missing <- unique(gold$compound_id[!gold$compound_id %in% pred$compound_id]) %>% trim_vec()
    errs["non_hgnc"] <- paste0("Missing compound_ids (up to 10 displayed): ", missing)
  }
}


bootstrapped_score <- function(prediction_path, gold_path){
  
  ###placeholder
  
}

score <- function(prediction_path, gold_path){
  
  pred <- readr::read_csv(prediction_path)
  gold <- readr::read_csv(gold_path)
  
  df <- dplyr::full_join(pred, gold)
  
  score <- function(x){} ###TBD
  
  return(score)
}

