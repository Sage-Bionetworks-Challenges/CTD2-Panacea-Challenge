library(argparse)
library(rjson)

source("R/validation_scoring.R")

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
    "--nullmodel1",
    type = "character",
    required = TRUE
)
parser$add_argument(
    "--nullmodel2",
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

scores <- score(args$inputfile, args$goldstandard, 
    args$nullmodel1, args$nullmodel2, args$round)

result_list <- list(
    "prediction_file_status" = "SCORED",
    "sc1_score" = scores$sc1,
    "sc2_score" = scores$sc2
)

result_list %>%  
    rjson::toJSON() %>% 
    write(args$results)