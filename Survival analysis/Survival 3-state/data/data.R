n_pat     <- 550                      # cohort size
n_years   <- 60                       # number of years 

generate  <- gen_data(n_pat,n_years)  # generates true, censored and OS/PFS data 
true_data <- generate$true_data       # stores the true data
sim_data  <- generate$sim_data        # stores the censored data
status    <- generate$status          # stores the censoring status
OS_PFS_data <- generate$OS_PFS_data   # store the OS / PFS structured data
