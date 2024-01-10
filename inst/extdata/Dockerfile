FROM rocker/tidyverse:4.0.2
LABEL maintainer="cansav09@gmail.com"
WORKDIR /rocker-build/

# Install a package from CRAN
RUN Rscript -e  "options(warn = 2);install.packages( \
    c('styler'), \
    repos = 'https://cloud.r-project.org/')"

# Set final workdir for commands
WORKDIR /home/rstudio
