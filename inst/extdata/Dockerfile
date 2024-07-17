FROM rocker/tidyverse:4.0.2
LABEL maintainer="cansav09@gmail.com"
WORKDIR /rocker-build/

# Install a package from CRAN
RUN Rscript -e  "options(warn = 2);install.packages( \
    c('remotes', \
      'styler'), \
    repos = 'https://cloud.r-project.org/')"

# Install the thing
RUN installGithub.r FredHutch/gimap
  
# Set final workdir for commands
WORKDIR /home/rstudio
