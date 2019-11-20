FROM rocker/tidyverse:3.6.1

# Install R packages.
RUN Rscript -e "install.packages('argparse')"
RUN Rscript -e "install.packages('rjson')"

COPY ../R/validation_scoring.R /usr/local/bin/validation_scoring.R
