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
parser$add_argument(
    "--metric",
    type = "character",
    required = TRUE
)

args <- parser$parse_args()

if (args$metric == "metric1") {
    source("/usr/local/bin/validation_scoring.R")
} else {
    source("/usr/local/bin/validation_scoring_subset.R")
}

scores <- score(args$inputfile, args$goldstandard, 
    args$nullmodel1, args$nullmodel2, args$round)

result_list <- list(
    "prediction_file_status" = "SCORED",
    "sc1_score" = scores[[1]],
    "sc2_score" = scores[[2]]
)

result_list %>%  
    rjson::toJSON() %>% 
    write(args$results)