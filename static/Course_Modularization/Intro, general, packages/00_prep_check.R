### First run the 00_prep.R script to download the desired R packages
### and to load the package you need
### from either CRAN or GitHub
## Next, run the code from below to check if you succesfully installed the packages 
library(pacman) 

# load (install if required) packages from CRAN
p_load("abind",  "dampack", "data.table", "DES", "devtools", "diagram", "dplyr", 
       "ellipse", "flexsurv", "flexsurvcure", "gdata", "gems", "grid", "gridExtra", 
       "igraph", "jsonlite", "knitr", "lazyeval", "lhs", 
       "markdown", "matrixStats", "mgcv", "msm", "mstate",
       "plotrix", "purrr", "psych", "reshape2", "rstudioapi",   
       "scales", "scatterplot3d", "stringr", "survHE", "survminer", "shiny",
       "tidyverse", "tidyr", "tm", "triangle", "truncnorm")
p_load("ggrepel")
p_load_gh("DARTH-git/OpenTree")
p_load_gh("DARTH-git/darthtools")



p_D <- darthtools::prob_to_rate(p = 0.5, t = 1)
v_p_D <- rtriangle(n = 1000, a = 0, b = 1, c = p_D)
hist(v_p_D)


# ## Base Case
# if you have a base case analysis, can use calculate_icers on that
data(hund_strat)
hund_icers <- calculate_icers(hund_strat$Cost,
                              hund_strat$QALYs,
                              hund_strat$Strategy)

plot(hund_icers) 

# Export this file
getwd() # where is your working directory
saveRDS(hund_icers, "data_hund_icers.rds")  # export the file
# you find the file at the location from the working directory

# Check which version of R you have, do you have 4.1.1?
getRversion()

 