gen_data <- function(n_pat, n_years)
  # Input:
  # n_pat  : number of patients
  # n_years: follow-up period
  # Output: generated survival data
{
  # specification of hazard functions to generate data from
  hazardf <- gems::generateHazardMatrix(n_states)
  colnames(hazardf@list.matrix) <- 
    rownames(hazardf@list.matrix) <- v_n
  
  # specifying the transition hazard from Healthy -> Sick
  hazardf[["Healthy","Sick"]] <- function (t, r1, r2){
    hweibull(t,r1, r2)
  }
  
  # specifying the transition hazard from Healthy -> Dead 
  hazardf[["Healthy","Dead"]] <- function (t, r1, r2){
    flexsurv::hgompertz(t,r1, r2)
  }
  
  # specifying the transition hazard from Sick -> Dead 
  hazardf[["Sick","Dead"]] <- function (t, r1, r2){
    hweibull(t,r1, r2)
  }
  
  # list of parameters for the hazard functions defined above
  mu        <- gems::generateParameterMatrix(hazardf) 
  rownames(mu@list.matrix) <- 
    colnames(mu@list.matrix) <- v_n
  
  mu[["Healthy", "Sick"]] <- list(1.5, 6)      #  the Weibull parameters for H -> S
  mu[["Healthy", "Dead"]] <- list(0.25, 0.08)  # the Gompertz params for H -> D
  mu[["Sick",    "Dead"]] <- list(0.5,4)       #  the Weibull parameters for S -> D
  
  
  
  # simulate the cohort
  cohort <- gems::simulateCohort(
    transitionFunctions = hazardf,
    parameters = mu,
    cohortSize = n_pat,
    to = n_years)
  
  # extract the simulated true data 
  true_data <- cohort@time.to.state
  colnames(true_data) <- v_n
  
  true_data$Dead[is.na(true_data$Dead)] <- n_years
  true_data$Sick[is.na(true_data$Sick)] <- true_data$Dead[is.na(true_data$Sick)]
  
  
  # create a status variable that will capture the transition events
  true_status         <- matrix(NA, nrow = n_pat, ncol = n_states, dimnames = list(1:n_pat,v_n))
  true_status         <- as.data.frame(true_status)
  true_status$Healthy <- ifelse(is.na(true_data$Healthy),0,1)
  true_status$Dead    <- ifelse(true_data$Dead == n_years, 0, 1)
  true_status$Sick    <- ifelse(true_data$Dead == true_data$Sick, 0, 1)
  
  
  censtime <- runif(n = n_pat, 0, n_years)
  
  censored_Sick <- ifelse(censtime      <= true_data$Sick |
                            true_data$Sick >  5, 1, 0)
  censored_Dead <- ifelse(censtime <= true_data$Dead|
                            true_data$Dead >5, 1, 0)
  
  sim_data <- true_data
  
  sim_data$Sick[censored_Sick == 1] <-  censtime[censored_Sick == 1]
  sim_data$Sick[sim_data$Sick >5 ]  <-  5
  
  sim_data$Dead[censored_Dead == 1] <-  censtime[censored_Dead == 1]
  sim_data$Dead[sim_data$Dead >5] <-  5
  
  status <- true_status
  
  status$Sick[censored_Sick == 1] = 0
  status$Dead[censored_Dead == 1] = 0
  
  # Usually trials report OS/PFS outcomes so we will recreate those
  
  OS_PFS_data <- data.frame(row.names = 1:n_pat)
  
  OS_PFS_data$PFS_time        <- apply(sim_data[, c("Sick","Dead")], 1, min) 
  OS_PFS_data$PFS_status      <- ifelse(status$Dead == 1 | status$Sick == 1, 1, 0 )
  
  OS_PFS_data$OS_time         <- sim_data$Dead
  OS_PFS_data$OS_status       <- status$Dead 
  list(cohort = cohort, true_data = true_data, true_status = true_status, 
       sim_data =  sim_data, status = status, OS_PFS_data = OS_PFS_data)
}

n_pat     <- 550                      # cohort size
n_years   <- 60                       # number of years 

generate  <- gen_data(n_pat,n_years)  # generates true, censored and OS/PFS data 
true_data <- generate$true_data       # stores the true data
sim_data  <- generate$sim_data        # stores the censored data
status    <- generate$status          # stores the censoring status
OS_PFS_data <- generate$OS_PFS_data   # store the OS / PFS structured data
