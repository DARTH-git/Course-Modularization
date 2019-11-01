################################################################################ 
# This Appendix 0 installs required packages if not located in R's local       #
# library                                                                      #
#                                                                              # 
# Author: Fernando Alarid-Escudero                                             # 
# E-mail: fernando.alarid@cide.edu                                             # 
################################################################################ 

### Function to check if packages are installed and if not, install them
f.install_and_load <- function(packages) {
  # Modified from https://www.listendata.com/2018/12/install-load-multiple-r-packages.html
  # The function below performs the following operations -
  #  - First it finds all the already installed R packages
  #  - Check packages which we want to install are already installed or not.
  #  - If package is already installed, it does not install it again.
  #  - If package is missing (not installed), it installs the package.
  #  - Loop through steps 2, 3 and 4 for multiple packages we want to install
  #  - Load all the packages (both already available and new ones).
  k <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(k)){
    install.packages(k, 
                    repos="https://cran.rstudio.com/", 
                    dependencies = TRUE)
    }
  
  for(package_name in packages){
    library(package_name,
            character.only = TRUE, 
            quietly = TRUE)
  }
}

### Install packages from CRAN
v.packages.to.install <- c("here", "dplyr","devtools","scales","ellipse","ggplot2", "lazyeval","igraph",
                           "truncnorm", "ggraph", "reshape2","knitr")

f.install_and_load(v.packages.to.install)

### Install dampack from GitHub
if (!require(dampack)) {devtools::install_github(repo = "DARTH-git/dampack")}; library(dampack)

### Install dectree from local directory (will install from GitHub once the package is ready)
if (!require(dectree)) {devtools::install_github(repo = "DARTH-git/dectree")}; library(dectree)




