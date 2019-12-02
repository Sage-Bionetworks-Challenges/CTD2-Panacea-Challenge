library(dplyr)
library(magrittr)
library(readr)
library(challengescoring)

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

# for testing
#
# template <- read_csv('template.csv')
# 
# challenge_targets <- template$target
# 
# pred <- template %>% 
#   mutate_at(vars(-target), ~ rnorm(length(.), mean = 0.5, sd = 0.2)) 


frac_overlap <- function(gold, pred){
  sum(gold %in% pred)/length(gold)
}

mean_rank <- function(gold, pred){
  sapply(gold, function(x){
    
  })
}

score <- function(prediction_path,
                  gold_path,
                  null_model_path,
                  gold_targets_path){
  
  targets <- read_csv(gold_targets_path)
  gold <- read_csv(gold_path) 
  pred <- read_csv(prediction_path) %>% 
    filter(target %in% targets)
  
  ###SC1 
  
  null_model <- read_rds('null_model.rds')
  
  gold_df <- gold %>% 
    select(-cmpd) %>% 
    group_by(cmpd_id) %>% 
    nest()
  
  pred_df <- pred %>% 
    gather(cmpd_id, confidence ,-target) %>% 
    group_by(cmpd_id) %>% 
    arrange(-confidence, target) %>% 
    slice(1:10) %>% 
    nest() ##instead of top n. We eliminate ties alphabetically!
  
  vals <- sapply(null_model, function(x){
    
    join <- inner_join(gold_df, pred_df, by = 'cmpd_id', suffix = c("_gold", "_pred")) %>% 
      inner_join(x, by = 'cmpd_id') %>% 
      mutate(gold_pred = map2(data_gold, data_pred, function(x,y){
        frac_overlap(x$target, y$target)
      })) %>% 
      mutate(gold_null = map2(data_gold, data, function(x,y){
        frac_overlap(x$target, y$target_random)
      })) %>% 
      select(cmpd_id, gold_pred, gold_null) %>% 
      unnest(c(gold_pred, gold_null)) %>% 
      ungroup() %>% 
      summarize(pval = t.test(x = gold_pred, y= gold_null, paired = T)$p.value) %>% 
      purrr::pluck('pval')
  
  if(is.nan(join)){
    join <- 1
  }
  join
  
  })
  
  sc1 <- mean(vals) %>% signif(3)

  ##SC2
  
  pred_sc2 <- pred %>% 
    gather(cmpd_id, confidence ,-target) %>% 
    group_by(cmpd_id) %>% 
    arrange(-confidence, target) %>% 
    nest()
  
  gold_sc2 <- gold %>% 
    select(-cmpd) %>% 
    group_by(cmpd_id) %>% 
    nest()
  
  sc2 <- gold_sc2 %>% 
    inner_join(pred_sc2, by = 'cmpd_id', suffix = c("_gold", "_pred")) %>% 
    mutate(gold_pred = map2(data_gold, data_pred, function(x,y){
      frac_overlap(x$target, y$target)
    })) 
  
  score <- c("sc1" = sc1, 
             "sc2" = sc2)
  
  return(score)
}


