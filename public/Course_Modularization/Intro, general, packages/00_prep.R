### Run all code in this script to download the desired R packages
### from either CRAN or GitHub
# The entire process will take several minutes
# Recommend running the code line by line, in the specified order below

# first download and use this package to conveniently install other packages
install.packages('pacman')
library(pacman) 

# load (install if required) packages from CRAN
p_load("devtools", "dplyr", "scales", "ellipse",
       "lazyeval", "igraph", "truncnorm", 
       "reshape2", "knitr", "markdown", "stringr", 
       "matrixStats", "lhs", "IMIS", "plotrix",
       "psych", "scatterplot3d", "grid", "mgcv", 
       "gridExtra", "gdata", "triangle", "flexsurv", "msm", 
       "survminer", "gems", "jsonlite", "rstudioapi", "tidyverse",
       "shiny", "tidyr", "purrr", "survHE", "mstate", "abind", "diagram", "dampack") 

# Enter an empty line to skip updates when prompted
install_version("ggplot2", version = "3.3.0", repos = "http://cran.us.r-project.org") 

# Install additional packages
p_load("ggrepel")

# load (install if required) packages from GitHub

# install_github("DARTH-git/darthtools", force = TRUE) # (Un)comment if there is a newer version#
p_load_gh("DARTH-git/darthtools")

# install_github("DARTH-git/DARTH-git/OpenTree", force = TRUE) # (Un)comment if there is a newer version
p_load_gh("DARTH-git/OpenTree")


# Install additional packages
p_load("ggraph") 

