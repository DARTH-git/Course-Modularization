gen_data <- function(n_pat, n_years)
  # Input:
  # n_pat  : number of patients
  # n_years: follow-up period
  # Output: generated survival data
{
  # specification of hazard functions to generate data from
  hazardf <- gems::generateHazardMatrix(n_s)
  colnames(hazardf@list.matrix) <- 
    rownames(hazardf@list.matrix) <- v_n
  
  # specifying the transition hazard from healthy -> sick
  hazardf[["healthy","sick"]] <- function (t, r1, r2){
    hweibull(t,r1, r2)
  }
  
  # specifying the transition hazard from healthy -> dead 
  hazardf[["healthy","dead"]] <- function (t, r1, r2){
    flexsurv::hgompertz(t,r1, r2)
  }
  
  # specifying the transition hazard from sick -> dead 
  hazardf[["sick","dead"]] <- function (t, r1, r2){
    hweibull(t,r1, r2)
  }
  
  # list of parameters for the hazard functions defined above
  mu        <- gems::generateParameterMatrix(hazardf) 
  rownames(mu@list.matrix) <- 
    colnames(mu@list.matrix) <- v_n
  
  mu[["healthy", "sick"]] <- list(1.5, 6)      #  the Weibull parameters for H -> S
  mu[["healthy", "dead"]] <- list(0.25, 0.08)  # the Gompertz params for H -> D
  mu[["sick",    "dead"]] <- list(0.5,4)       #  the Weibull parameters for S -> D
  
  
  
  # simulate the cohort
  cohort <- gems::simulateCohort(
    transitionFunctions = hazardf,
    parameters = mu,
    cohortSize = n_pat,
    to = n_years)
  
  # extract the simulated true data 
  true_data <- cohort@time.to.state
  colnames(true_data) <- v_n
  
  true_data$dead[is.na(true_data$dead)] <- n_years
  true_data$sick[is.na(true_data$sick)] <- true_data$dead[is.na(true_data$sick)]
  
  
  # create a status variable that will capture the transition events
  true_status         <- matrix(NA, nrow = n_pat, ncol = n_s, dimnames = list(1:n_pat,v_n))
  true_status         <- as.data.frame(true_status)
  true_status$healthy <- ifelse(is.na(true_data$healthy),0,1)
  true_status$dead    <- ifelse(true_data$dead == n_years, 0, 1)
  true_status$sick    <- ifelse(true_data$dead == true_data$sick, 0, 1)
  
  
  censtime <- runif(n = n_pat, 0, n_years)
  
  censored_sick <- ifelse(censtime      <= true_data$sick |
                            true_data$sick >  5, 1, 0)
  censored_dead <- ifelse(censtime <= true_data$dead|
                            true_data$dead >5, 1, 0)
  
  sim_data <- true_data
  
  sim_data$sick[censored_sick == 1] <-  censtime[censored_sick == 1]
  sim_data$sick[sim_data$sick >5 ]  <-  5
  
  sim_data$dead[censored_dead == 1] <-  censtime[censored_dead == 1]
  sim_data$dead[sim_data$dead >5] <-  5
  
  status <- true_status
  
  status$sick[censored_sick == 1] = 0
  status$dead[censored_dead == 1] = 0
  
  # Usually trials report OS/PFS outcomes so we will recreate those
  
  OS_PFS_data <- data.frame(row.names = 1:n_pat)
  
  OS_PFS_data$PFS_time        <- apply(sim_data[, c("sick","dead")], 1, min) 
  OS_PFS_data$PFS_status      <- ifelse(status$dead == 1 | status$sick == 1, 1, 0 )
  
  OS_PFS_data$OS_time         <- sim_data$dead
  OS_PFS_data$OS_status       <- status$dead 
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
