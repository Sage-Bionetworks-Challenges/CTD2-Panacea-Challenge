FROM rocker/tidyverse:3.6.1

# Install R packages.
RUN Rscript -e "install.packages('argparse')"
RUN Rscript -e "install.packages('rjson')"

COPY R/validation_scoring_subset.R /usr/local/bin/validation_scoring_subset.R
COPY validate.R /usr/local/bin/validate.R
COPY score.R /usr/local/bin/score.R
COPY models/*.rds /models/