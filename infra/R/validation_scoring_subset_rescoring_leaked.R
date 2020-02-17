suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(purrr))

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
    errs["ncol_long"] <- paste0("Prediction file has extra cols: ", ncol(pred), " cols detected.")
  }
  
  if(nrow(pred)<nrow_req){
    errs["nrow_short"] <- paste0("Prediction file is missing rows. Only ", nrow(pred), " rows detected.")
  }
  
  if(nrow(pred)>nrow_req){
    errs["nrow_long"] <- paste0("Prediction file has extra rows: ", nrow(pred), " rows detected.")
  }
  
  if(!identical(colnames(pred), colnames_req)){
    errs["colnames"] <- paste0("Column names are not correct. Column names must be ", paste(c(colnames_req), collapse=", "))
  }
  

  tryCatch({
    if(!all(pred[-1] >= 0) | !all(pred[-1] <= 1)){
      errs["wrong_range"] <- paste0("Confidence values are not between 0 and 1.")
    }
  }, error = function(e) {
     return(NA)
  })
  
  if(any(apply(pred[-1], 1:2, is.na))){
    errs["is_na"] <- paste0("One or more predictions are missing.")
  }
  
  if(!all(apply(pred[-1], 1:2, is.numeric))){
    errs["non_numeric"] <- paste0("Predictions are not all numeric values.")
  }
  
  if(any(!pred$target %in% target_ids)){
    invalid <- unique(pred$target[!pred$target %in% target_ids]) %>% trim_vec()
    errs["non_target"] <- paste0("Invalid target identifiers were included in your prediction file (up to 10 displayed): ", invalid)
  }
  
  return(errs)
}

 
frac_overlap <- function(gold, pred){
  sum(gold %in% pred)/length(gold)
}

 
score <- function(prediction_path,
                  gold_path,
                  null_model_path_sc1, 
                  null_model_path_sc2,
                  round = c("leaderboard", "final"),
                  no_rand = 10){
  
  #most obnoxiously long vector of targets here
  
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
  
  if(round == "leaderboard"){
    random_sample <- no_rand
  }else if(round == "final"){
    random_sample <- expr(length(null_model)) ##to be evaluated later
  }
  
  # ground truth leaked in biorxiv preprint
  leaked <- c("cmpd_CB","cmpd_DC","cmpd_DF","cmpd_EY","cmpd_FC","cmpd_FI","cmpd_IV","cmpd_LA",
              "cmpd_PQ","cmpd_QJ","cmpd_RH","cmpd_WB","cmpd_VQ","cmpd_XY")
  "%ni%" <- Negate("%in%")
  
  # remove the leaked compounds from datasets
  gold <- suppressMessages(read_csv(gold_path)) %>% filter(target %in% targs) %>% filter(cmpd_id %ni% leaked)
  pred <- suppressMessages(read_csv(prediction_path)) %>% filter(target %in% targs) %>% select(-leaked)
  
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
    group_by(cmpd_id) %>% 
    arrange(-confidence, target) %>% 
    slice(1:10) %>%  ##instead of top n. We eliminate ties alphabetically!
    nest() %>% 
    arrange(cmpd_id)
  
  sc1_vals <- sapply(null_model, function(x){
    x <- x %>% filter(cmpd_id %ni% leaked)
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
      summarize(pval = suppressWarnings(wilcox.test(x = gold_pred, y= gold_null, paired = T, exact = NULL)$p.value)) %>% 
      purrr::pluck('pval')
  

  if(is.nan(join)){
    join <- 1
  }
  join
  
  })
  
  sc1 <- mean(-log2(sc1_vals)) %>% signif(5)

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
    x <- x %>% filter(cmpd_id %ni% leaked)
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
      summarize(pval = suppressWarnings(wilcox.test(x = gold_pred, y= gold_null, paired = T, exact = NULL)$p.value)) %>% 
      purrr::pluck('pval')
    
    if(is.nan(join)){
      join <- 1
    }
    join
    
  })
  
  sc2 <- mean(-log2(sc2_vals)) %>% signif(5)
  
  score <- c("sc1" = sc1, 
             "sc2" = sc2)
  
  return(score)
}



