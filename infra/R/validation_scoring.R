suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(challengescoring))

trim_vec <- function(vec, trim = 10){
  if(length(vec) > trim){
    vec <- vec[1:trim]
    vec <- as.character(vec)
    vec[trim+1] <- '...'
  }else{
    vec
  }
}

validate <- function(prediction_path, template_path){
  
  pred <- readr::read_csv(prediction_path)
  temp <- readr::read_csv(template_path)
  
  ###configure validation
  ncol_req <- ncol(temp)
  nrow_req <- nrow(temp)
  colnames_req <- colnames(temp)
  target_ids <- unique(temp$target)
  

  errs <- list()

  if(ncol(pred)<ncol_req){
    errs["ncol_short"] <- paste0("Prediction file is missing cols. Only ", ncol(pred), " cols detected.")
  }
  
  if(ncol(pred)>ncol_req){
    errs["ncol_long"] <- paste0("Prediction file has extra  cols ", ncol(pred), " cols detected.")
  }
  
  if(nrow(pred)<nrow_req){
    errs["nrow_short"] <- paste0("Prediction file is missing rows Only ", nrow(pred), " rows detected.")
  }
  
  if(nrow(pred)>nrow_req){
    errs["nrow_long"] <- paste0("Prediction file has extra  rows ", nrow(pred), " rows detected.")
  }
  
  if(isTRUE(colnames(pred) %in% template)){
    errs["colnames"] <- paste0("Column names are not correct. Column names must be ", cat(colnames_req))
  }
  
  if(!all(pred[-1] > 0) | !all(pred[-1] < 1)){
    errs["wrong_range"] <- paste0("Confidence values are not between 0 and 1.")
  }
  
  if(!all(apply(pred[-1], 1:2, is.numeric))){
    errs["non_numeric"] <- paste0("Predictions are not all numeric values.")
  }
  
  if(any(!pred$target %in% target_ids)){
    invalid <- unique(pred$target[!pred$target %in% target_ids]) %>% trim_vec()
    errs["non_target"] <- paste0("Invalid target identifiers were included in your prediction file (up to 10 displayed): ", invalid)
  }
}

 
frac_overlap <- function(gold, pred){
  sum(gold %in% pred)/length(gold)
}

score <- function(prediction_path,
                  gold_path,
                  null_model_path_sc1, 
                  null_model_path_sc2,
                  round = c("leaderboard", "final")){
  
  if(round == "leaderboard"){
    random_sample <- 5
  }else if(round == "final"){
    random_sample <- expr(length(null_model)) ##to be evaluated later
  }

  gold <- suppressMessages(read_csv(gold_path))
  pred <- suppressMessages(read_csv(prediction_path))
  
  ###SC1 
  
  null_model <- read_rds(null_model_path_sc1)
  null_model <- sample(null_model, eval(random_sample))
  
  gold_df <- gold %>% 
    select(-cmpd) %>% 
    group_by(cmpd_id) %>% 
    nest() %>% 
    arrange(cmpd_id)
  
  pred_df <- pred %>% 
    gather(cmpd_id, confidence ,-target) %>% 
    arrange(cmpd_id) %>% 
    group_by(cmpd_id) %>% 
    arrange(-confidence, target) %>% 
    slice(1:10) %>% 
    nest() %>% 
    arrange(cmpd_id) ##instead of top n. We eliminate ties alphabetically!
  
  sc1_vals <- sapply(null_model, function(x){
    
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
  
  sc1 <- mean(-log10(sc1_vals)) %>% signif(5)

  ##SC2
  
  null_model <- read_rds(null_model_path_sc2) 
  null_model <- sample(null_model, eval(random_sample))  
  
  pred_sc2 <- pred %>% 
    gather(cmpd_id, confidence ,-target) %>% 
    group_by(cmpd_id) %>% 
    arrange(-confidence, target) %>% 
    nest() %>% 
    arrange(cmpd_id)
  
  gold_sc2 <- gold %>% 
    select(-cmpd) %>% 
    group_by(cmpd_id) %>% 
    nest() %>% 
    arrange(cmpd_id) 
  
  sc2_vals <- sapply(null_model, function(x){
    
    join <- inner_join(gold_sc2, pred_sc2, by = 'cmpd_id', suffix = c("_gold", "_pred")) %>% 
      inner_join(x, by = 'cmpd_id') %>% 
      mutate(gold_pred = map2(data_gold, data_pred, function(x,y){
        match(x$target, y$target) %>% unique()
      })) %>% 
      mutate(gold_null = map2(data_gold, data, function(x,y){
        match(x$target, y$target_random) %>% unique()
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
  
  sc2 <- mean(-log10(sc2_vals)) %>% signif(5)
  
  score <- c("sc1" = sc1, 
             "sc2" = sc2)
  
  return(score)
}



