## Function to fit multiple functional forms to survival data using survHE
fit.fun <- function(time, status, data = data, extrapolate = FALSE, times = NULL) {
  # Input
    # time: follow-up time in survival data
    # status: survival status
    # data: survival data
    # extrapolate: whether to extrapolate the survival data or not
    # times: time horizon of extrapolation
  # Output
    # res: a list w/ model object, model results and model AIC, BIC
  require(survHE)
  # Extract the right data columns 
  data$time   <- data[,   time]  
  data$status <- data[, status]  
  
  # Specify time horizon based on whether to extrapolate
  if (extrapolate == TRUE)  {
    plot.times <- max(times)
  } else if  (extrapolate == FALSE) {
    plot.times <- max(data$time)
  }
  
  # Fit parametric survival models
  # Define the vector of models to be used
  mods <- c("exp", "weibull", "gamma", "lnorm", "llogis", "gengamma") 
  # Run the models using MLE via flexsurv
  fit.survHE <- fit.models(formula = Surv(time, status) ~ 1, data = data, distr = mods)
  
  # Extrapolate all models beyond the KM curve and plot
  print(plot(fit.survHE, add.km = T, t = times))
  
  # Compare AIC values
  AIC <- fit.survHE$model.fitting$aic
  AIC <- round(AIC,3)
  
  # Compare BIC values
  BIC <- fit.survHE$model.fitting$bic
  BIC <- round(BIC,3)
  names(AIC) <- names(BIC) <- names(fit.survHE$models)
  
  # Store and return results
  res <- list(fit.survHE  = fit.survHE,
              models      = fit.survHE$models,
              AIC         = AIC,
              BIC         = BIC)
  res
}

## Function to build partitioned survival model with PSA using survHE
partsurv <- function(pfs_survHE, os_survHE, choose_PFS, choose_OS, time = times, n_sim = 100, seed = 421){
  # Input
    # pfs_survHE: survHE obj fitting pfs
    # os_survHE : survHE obj fitting os
    # choose_PFS: preferred PFS distribution
    # choose_OS: preferred OS distribution
    # time = numeric vector of time to estimate probabilities
    # n_sim = number of PSA simulations
    # seed = random number generator seed
  # Output
    # res: a list w/ the cohort trace, the confidence intervals and the mean of the survival probabilities
  set.seed(seed)

  deter <- ifelse(n_sim == 1, 1, 0)
  n_sim <- ifelse(deter == 1, 1000, n_sim)
  
  mod.pfs <- names(pfs_survHE$fit.survHE$models)
  mod.os  <- names(os_survHE$fit.survHE$models)

  # Conduct PSA on survival models
  # fFor PFS 
  fit_PFS <- make.surv(pfs_survHE$fit.survHE,
                       mod = which(mod.pfs == choose_PFS),
                       nsim = n_sim, 
                       t = times)
  # For OS
  fit_OS  <- make.surv(os_survHE$fit.survHE,
                       mod = which(mod.os == choose_OS),
                       nsim = n_sim, 
                       t = times)
        
  pfs.surv <- fit_PFS$mat[[1]][,-1]
  os.surv  <- fit_OS$mat[[1]][,-1]
  
  sick           <- os.surv - pfs.surv  # estimate the probability of remaining in the progressed state
  sick[sick < 0] <- 0                   # in cases where the probability is negative replace with zero
  healthy        <- pfs.surv            # probability of remaining stable
  dead           <- 1 - os.surv         # probability of being dead
  trace <- abind(healthy,
                 sick,
                 dead,rev.along = 0)
  
  trace <- aperm(trace,perm = c(1,3,2))
  
  if(deter ==1 ){  # if running deterministic analysis
    trace = apply(trace,1:2,mean)
    CI <- NULL
    mean.trace = NULL
    } else{  # if running PSA
   # Construct 95% confidence intervals for the survival curves
   CI <- apply(trace,1:2, quantile, probs = c(0.025, 0.975)) 
   mean.trace <- apply(trace,1:2,mean)
   CI = aperm(CI,perm = c(2,3,1))
   dimnames(CI)[[3]] = c("low", "high")
   dimnames(CI)[[2]] <- c("healthy", "sick", "dead")
   dimnames(mean.trace)[[2]] <- c("healthy", "sick", "dead")
  }
  dimnames(trace)[[2]] <- c("healthy", "sick", "dead")
  res   <- list(trace = trace, CI = CI, Mean = mean.trace)

  return(res)
}

## Function to generate survival data
gen_data <- function(n_pat, n_years) {
  # Input
    # n_pat: number of patients to simulate
    # n_years: numer of years of follow-up to simulate
  # Output
    # a list w/ simulated and true survival data
  
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


