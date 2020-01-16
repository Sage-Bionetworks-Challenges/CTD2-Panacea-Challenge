library(argparse)
library(rjson)

parser = ArgumentParser()

parser$add_argument(
    "--inputfile",
    type = "character",
    required = TRUE
)
parser$add_argument(
    "--goldstandard",
    type = "character",
    required = TRUE
)
parser$add_argument(
    "--round",
    type = "character",
    required = TRUE
)
parser$add_argument(
    "--results",
    type = "character",
    required = TRUE
)

args <- parser$parse_args()

source("/usr/local/bin/validation_scoring_subset.R")
nullmodel1 <- "/models/null_model_sc1_leaderboard_subset.rds"
nullmodel2 <- "/models/null_model_sc2_leaderboard_subset.rds"
if (args$round == "final") {
    nullmodel1 <- "/models/null_model_sc1_subset.rds"
    nullmodel2 <- "/models/null_model_sc2_subset.rds"
}

scores <- score(args$inputfile, args$goldstandard, 
    nullmodel1, nullmodel2, args$round)

result_list <- list(
    "prediction_file_status" = "SCORED",
    "sc1_score" = scores[[1]],
    "sc2_score" = scores[[2]]
)

result_list %>%  
    rjson::toJSON() %>% 
    write(args$results)