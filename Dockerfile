FROM rocker/rstudio:4.1.2

RUN git clone https://github.com/schalkdaniel/simulations-distr-auc /home/rstudio/simulations-distr-auc
RUN Rscript /home/rstudio/simulations-distr-auc/R/install-pkgs-versions.R
RUN echo "setwd('/home/rstudio/simulations-distr-auc')" >> /home/rstudio/.Rprofile
