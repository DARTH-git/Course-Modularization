## Function to fit multiple functional forms to survival data using survHE
fit.fun <- function(time, status, data = data , add = FALSE, extrapolate = FALSE, times)  
{
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
              models  = fit.survHE$models,
              AIC = AIC,
              BIC = BIC)
  res
}


fit.mstate <- function(time, status, trans,  data = data , add = FALSE, extrapolate = FALSE, times)  
{
  data$time  <- data[, time  ]
  data$tatus <- data[, status]

    if (extrapolate == TRUE)  {
    plot.times <- max(times)
  } else if  (extrapolate == FALSE) {
    plot.times <- max(data$time)
  }
  
  # Progression free survival  
  KM.fit     <-     survfit(Surv(time, status) ~ trans , data = data)                                 # fit Kaplan-Meier curve 
  fit.llogis <- flexsurvreg(Surv(time, status) ~ trans + shape(trans), data = data, dist = "llogis" ) # fit model with loglogistic distribution
  fit.weib   <- flexsurvreg(Surv(time, status) ~ trans + shape(trans), data = data, dist = "weibull") # fit model with Weibull distribution
  fit.lnorm  <- flexsurvreg(Surv(time, status) ~ trans + sdlog(trans), data = data, dist = "lnorm"  ) # fit model with lognormal distribution
  fit.gamma  <- flexsurvreg(Surv(time, status) ~ trans + shape(trans), data = data, dist = "gamma"  ) # fit model with gamma distribution 
  fit.gengamma  <- flexsurvreg(Surv(time, status) ~ trans + Q(trans) + sigma(trans), data = data, dist = "gengamma"  ) # fit model with gamma distribution 
  
  
  # extarapolate all models beyond the KM curve
  if(add){ lines(KM.fit, ylab = "Survival Probability", xlab = "Time", ylim = c(0,1), xlim = c(0, plot.times), conf.int= F)}
  if(!add){ plot(KM.fit, ylab = "Survival Probability", xlab = "Time", ylim = c(0,1), xlim = c(0, plot.times), conf.int= F)}
  lines(fit.llogis,   t = times, col = 2, ci = F)
  lines(fit.weib,     t = times, col = 3, ci = F)
  lines(fit.lnorm,    t = times, col = 4, ci = F)
  lines(fit.gamma,    t = times, col = 5, ci = F)
  lines(fit.gengamma,    t = times, col = 6, ci = F)
  legend("topright", cex = 0.7, c("Kaplan-Meier", "Loglogistic", "Weibull", "Lognormal", "Gen.Gamma"), col = 1:5, lty = rep(1, 5), bty="n")
  
  # compare AIC values
  AIC <- c(    Loglogistic = AIC(fit.llogis),                                         
               Weibull     = AIC(fit.weib), 
               Lognormal   = AIC(fit.lnorm), 
               Gamma       = AIC(fit.gamma),
               GenGamma    = AIC(fit.gengamma))
  AIC= round(AIC,3)
  
  # compare BIC values
  BIC <- c(    Loglogistic = BIC(fit.llogis),                                         
               Weibull     = BIC(fit.weib), 
               Lognormal   = BIC(fit.lnorm), 
               Gamma       = BIC(fit.gamma),
               GenGamma    = BIC(fit.gengamma))
  
  BIC <- round(BIC,3)
  
  res <- list(Loglogistic = fit.llogis,
              Weibull     = fit.weib,
              Lognormal   = fit.lnorm, 
              Gamma       = fit.gamma,
              GenGamma    = fit.gengamma,
              AIC         = AIC,
              BIC         = BIC)
  res
}


trace.DES = function(msm_sim = des_sim, tmat, n_i, times )
{
  # Restructure the data to extract markov trace
  data.mstate.sim <- data.frame(cbind(matrix(t(msm_sim$st), ncol=1),
                                      matrix(t(msm_sim$t) , ncol=1)))
  colnames(data.mstate.sim) <- c("state","time")
  data.mstate.sim$subject <- rep(1:n_i, each = ncol(msm_sim$st))
  
  data.mstate.sim = na.omit(data.mstate.sim)
  data.mstate.sim = data.mstate.sim[!duplicated(data.mstate.sim), ] # remove duplicate entries in the dataset
  
  # create transition intensitiy matrix with initial values based on the structure of tmat
  twoway7.q               <- tmat
  twoway7.q[!is.na(tmat)] <- 0.5
  twoway7.q[is.na(tmat)]  <- 0
  # fit msm model only so that we can extract the prevalence (i.e. trace) thrrough the prevalence.msm function
  
  fit.msm.sim <- msm(state ~ time,subject = subject, data = data.mstate.sim, qmatrix = twoway7.q, 
                     exacttimes = T, use.deriv = TRUE, analyticp = FALSE, fixedpars = TRUE, hessian = F)
 
  M.tr.des <- prevalence.msm(fit.msm.sim, times = times) # Markov trace when DES model is used
  
  return(M.tr.des[[3]]/100)
}


partsurv <- function(pfs_survHE, os_survHE,  choose_PFS, choose_OS, time = times, n_sim = 1, seed = 421){
  # Input
  # pfs_survHE: survHE obj fitting pfs
  # os_survHE : survHE obj fitting os
  # choose_PFS: preferred PFS distribution
  # choose_OS: preferred OS distribution
  # time = numeric vector of time to estimate probabilities
  # n_sim = number of PSA simulations
    # output:
  #  res a list w/ one entry of a data frame w/ probabilities associated w/ stable ,prog and dead.
  set.seed(seed)

  deter <- ifelse(n_sim ==1,1,0)
  n_sim <- ifelse(deter == 1, 1000,n_sim)
    mod.pfs <- names(pfs_survHE$fit.survHE$models)
    mod.os  <- names(os_survHE$fit.survHE$models)

        fit_PFS   <- make.surv(pfs_survHE$fit.survHE,
                          mod = which(mod.pfs == choose_PFS),
                          nsim = n_sim, 
                          t = times)
        fit_OS    <- make.surv(os_survHE$fit.survHE,
                               mod = which(mod.os == choose_OS),
                               nsim = n_sim, 
                               t = times)
        
    pfs.surv <- fit_PFS$mat[[1]][,-1]
    os.surv  <- fit_OS$mat[[1]][,-1]


  
  sick                 <- os.surv - pfs.surv      # estimate the probability of remaining in the progressed state
  sick[sick < 0]       <- 0                       # in cases where the probability is negative replace with zero
  healthy              <- pfs.surv                # probability of remaining stable
  dead                 <- 1 - os.surv             # probability of being dead
  trace <- abind( healthy,
                  sick,
                  dead,rev.along = 0)
  
  
      trace = aperm(trace,perm = c(1,3,2))
  #browser()
  if(deter ==1 ){
    trace = apply(trace,1:2,mean)
    CI <- NULL
    mean.trace = NULL
    }else{  
   CI <- apply(trace,1:2, quantile, probs = c(0.025,0.975)) 
   mean.trace <- apply(trace,1:2,mean)
   CI = aperm(CI,perm = c(2,3,1))
   dimnames(CI)[[3]] = c("low", "high")
   dimnames(CI)[[2]] <- c("healthy",   "sick",  "dead")
   dimnames(mean.trace)[[2]] <- c("healthy",   "sick",  "dead")
   
    }
      
  dimnames(trace)[[2]] <- c("healthy",   "sick",  "dead")
  res   <- list(trace = trace, CI = CI, Mean = mean.trace)

  return(res)
}


flexsurvreg_prob <- function(object, newparams = NULL, times){

  if(is.null(newparams) == T ){
  params <- object$res[,1]  
  params <- as.matrix(t(params))
  }else {
    params <- newparams 
    params <- as.matrix(params)
    }

  if (ncol(params)== 1){
  surv <- object$dfns$p(times, params[,1], lower.tail = F)
  }else if (ncol(params)== 2){
    surv <- object$dfns$p(times,params[,1],params[,2], lower.tail = F)
   }else if (ncol(params)== 2){
     surv <- object$dfns$p(times,params[,1],params[,2],params[,3], lower.tail = F)
   } else{
     surv <- object$dfns$p(times,params[,1],params[,2],params[,3], lower.tail = F)
   }
   
  t.p <- 1- surv[-1]/(surv[-length(surv)])
return(t.p = t.p)
  }


gen_data <- function(n_pat, n_years)
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
        sim_data =  sim_data,      status = status, OS_PFS_data = OS_PFS_data)
}


