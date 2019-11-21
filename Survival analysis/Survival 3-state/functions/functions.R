## Function to fit multiple functional forms to survival data
fit.fun <- function(time, status, data = data , add = FALSE, extrapolate = FALSE, times)  
{
  #Extact the right data columns 
  data$time   <-   data[,   time]  
  data$status <-   data[, status]  

    if (extrapolate == TRUE)  {
    plot.times <- max(times)
  } else if  (extrapolate == FALSE) {
    plot.times <- max(data$time)
  }
  
  # Progression free survival  
  KM.fit     <-     survfit(Surv(time, status) ~ 1, data = data)                   # fit Kaplan-Meier curve 
  fit.llogis <- flexsurvreg(Surv(time, status) ~ 1, data = data, dist = "llogis" ) # fit model with loglogistic distribution
  fit.weib   <- flexsurvreg(Surv(time, status) ~ 1, data = data, dist = "weibull") # fit model with Weibull distribution
  fit.lnorm  <- flexsurvreg(Surv(time, status) ~ 1, data = data, dist = "lnorm"  ) # fit model with lognormal distribution
  fit.gamma  <- flexsurvreg(Surv(time, status) ~ 1, data = data, dist = "gamma"  ) # fit model with gamma distribution 
  fit.exp    <- flexsurvreg(Surv(time, status) ~ 1, data = data, dist = "exp"    ) # fit model with exponential distribution
  fit.gengamma  <- flexsurvreg(Surv(time, status) ~ 1, data = data, dist = "gengamma"  ) # fit model with gamma distribution  
  
  
  # extarapolate all models beyond the KM curve
  if(add){ lines(KM.fit, ylab = "Survival Probability", xlab = "Time", ylim = c(0,1), xlim = c(0, plot.times), conf.int= F)}
  if(!add){ plot(KM.fit, ylab = "Survival Probability", xlab = "Time", ylim = c(0,1), xlim = c(0, plot.times), conf.int= F, mark.time= T)}
  lines(fit.llogis,   t = times, col = 2, ci = F)
  lines(fit.weib,     t = times, col = 3, ci = F)
  lines(fit.lnorm,    t = times, col = 4, ci = F)
  lines(fit.gamma,    t = times, col = 5, ci = F)
  lines(fit.gengamma,    t = times, col = 6, ci = F)
  lines(fit.exp,      t = times, col = 7, ci = F)
  legend("topright", cex = 0.7, c("Kaplan-Meier", "Loglogistic", "Weibull", "Lognormal", "Gamma","GenGamma", "Exponential"), col = 1:7, lty = rep(1, 7), bty="n")
  
  # compare AIC values
  AIC <- c(    Loglogistic = AIC(fit.llogis),                                         
               Weibull     = AIC(fit.weib), 
               Lognormal   = AIC(fit.lnorm), 
               Gamma       = AIC(fit.gamma),
               GenGamma       = AIC(fit.gengamma),
               Exponentail = AIC(fit.exp))
  AIC= round(AIC,3)
  
  # compare BIC values
  BIC <- c(    Loglogistic = BIC(fit.llogis),                                         
               Weibull     = BIC(fit.weib), 
               Lognormal   = BIC(fit.lnorm), 
               Gamma       = BIC(fit.gamma),
               GenGamma    = BIC(fit.gengamma),
               Exponential = BIC(fit.exp))
  
  BIC <- round(BIC,3)
  
  res <- list(Loglogistic = fit.llogis,
              Weibull     = fit.weib,
              Lognormal   = fit.lnorm, 
              Gamma       = fit.gamma,
              GenGamma       = fit.gengamma,
              Exponential = fit.exp, 
              AIC         = AIC,
              BIC         = BIC)
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
  KM.fit     <-     survfit(Surv(time, status) ~ trans , data = data)                   # fit Kaplan-Meier curve 
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


partsurv <- function(fit.pfs, fit.os, time = times){
  # Input
  # fit.pfs: flexsurv obj fitting pfs
  # fit.os: flexsurv obj fitting os
  # title:
  # time = numeric vector of time to estimate probabilities
  # output:
  #  res a list w/ one entry of a data frame w/ probabilities associated w/ stable ,prog and dead.
  
  pfs.surv <- summary(fit.pfs, t = time, ci = F)[[1]]$est
  os.surv  <- summary(fit.os,  t = time, ci = F)[[1]]$est
  sick                 <- os.surv - pfs.surv      # estimate the probability of remaining in the progressed state
  sick[sick < 0]       <- 0                       # in cases where the probability is negative replace with zero
  healthy               <- pfs.surv                # probability of remaining stable
  dead                 <- 1 - os.surv             # probability of being dead
  trace <- cbind(healthy, sick, dead)
  res   <- list(trace = trace)

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
  
  mu[["healthy", "sick"]] <- list(1.5, 6)   #  the Weibull parameters for H -> S
  mu[["healthy", "dead"]] <- list(0.25, 0.08)  # the Gompertz params for H -> D
  mu[["sick",    "dead"]] <- list(0.5,4)     #  the Weibull parameters for S -> D
  
  
  
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

samplev <- function(m.Probs, m) {
  # Arguments
  # m.Probs: matrix with probabilities (n.i * n.s)
  # m:       number of states than need to be sampled per individual  
  # Return
  # ran:    n.i x m matrix filled with sampled health state(s) per individual
  
  d <- dim(m.Probs)  # dimensions of the matrix filled with the multinomical probabilities for the health states 
  n <- d[1]          # first dimension - number of rows (number of individuals to sample for)
  k <- d[2]          # second dimension - number of columns (number of health states considered)
  lev <- dimnames(m.Probs)[[2]]  # extract the names of the health states considered for sampling
  if (!length(lev))  # in case names for the health states are missing, use numbers to specify the health states
    lev <- 1:k       # create a sequence from 1:k (number of health states considered)
  # create a matrix 
  ran <- matrix(lev[1], ncol = m, nrow = n) # create the matrix ran, filled with the first health state of the levels 
  U <- t(m.Probs)    # transposed m.Probs matrix n.i x n.s --> n.s x n.i 
  
  for(i in 2:k) {    # start loop, from the 2nd health states
    U[i, ] <- U[i, ] + U[i - 1, ] # start summing the probabilities of the different health states per individual 
  }
  if (any((U[k, ] - 1) > 1e-05))  # sum of all probs per individual - 1 should be 0 (use 1e-05 for rounding issues), else print the error statement
    stop("error in multinom: probabilities do not sum to 1")
  
  for (j in 1:m) {   # start loop of the state that needs to be sampled (m)
    un <- rep(runif(n), rep(k, n))       # sample from a uniform distribution of length n*k
    ran[, j] <- lev[1 + colSums(un > U)] # store the health state at the jth column of the U matrix
  }
  ran # return the new health state per individual n.i x m
} # close the function 


#plot health state trace
plot_m_TR <- function(m_M) {
  # plot the distribution of the population across health states over time (trace)
  # count the number of individuals in each health state at each cycle
  m_TR <- t(apply(m_M, 2, function(x) table(factor(x, levels = v_n, ordered = TRUE)))) 
  m_TR <- m_TR / n_i                                       # calculate the proportion of individuals 
  colnames(m_TR) <- v_n                                    # name the rows of the matrix
  rownames(m_TR) <- paste("Cycle", 0:n_t, sep = " ")       # name the columns of the matrix
  # Plot trace of first health state
  matplot(m_TR, type = "l", main = "Health state trace", col= 1:n_s,
          ylim = c(0, 1), ylab = "Proportion of cohort", xlab = "Cycle")
  legend("topright", v_n, col = 1:n_s,    # add a legend to current plot
         lty = rep(1, 3), bty = "n", cex = 0.65)
  
}
