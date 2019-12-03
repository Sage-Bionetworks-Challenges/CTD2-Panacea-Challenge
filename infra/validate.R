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
    "--template",
    type = "character",
    required = TRUE
)
parser$add_argument(
    "--results",
    type = "character",
    required = TRUE
)

args <- parser$parse_args()

result_list <- list(
    "prediction_file_status" = "VALID",
    "prediction_file_errors" = "",
    "round" = 1
)

errors <- validate(args$inputfile, args$template)

if (length(errors) > 0) {
    result_list$prediction_file_status = "INVALID"

    errors <- paste(unlist(errors, use.names=F), collapse="\n")
    result_list$prediction_file_errors = errors
}

result_list %>%  
    rjson::toJSON() %>% 
    write(args$results)