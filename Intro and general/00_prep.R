# Welcome to the CE16 course

### Run all code in this script to download the desired R packages 
### from either CRAN or GitHub

# first download and use this package to conveniently install other packages
if (!require('pacman')) install.packages('pacman'); library(pacman) 

# load (install if required) packages from CRAN
p_load("here", "devtools", "dplyr", "scales", "ellipse", "ggplot2", "lazyeval", "igraph", "truncnorm", "ggraph", "reshape2", "knitr", "markdown","stringr")

# load (install if required) packages from GitHub
install_github("DARTH-git/dampack", force = TRUE) # (Un)comment if there is a newer version
p_load_gh("DARTH-git/dampack")


# DECISION TREE
install_github("DARTH-git/dectree", force = TRUE) # (Un)comment if there is a newer version
p_load_gh("DARTH-git/dectree")   


### Below packages are module-specific are ar not relavent at this point. 

# SURVIVAL ANALYSIS
# p_load("gems", "flexsurv", "survHE", "msm", "mstate")

# CALIBRATION
# p_load("lhs", "IMIS", "matrixStats", "plotrix", "psych", "scatterplot3d")

# VOI
#p_load("grid", "mgcv", "gridExtra", "gdata")
          
  
          