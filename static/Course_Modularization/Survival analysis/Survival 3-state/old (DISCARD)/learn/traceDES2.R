trace.DES2 <- function (cohorts, times, M = 100, stateNames = paste("State", 
                                                    as.list(1:dim(cohorts)[1]))) 
{
  require(gems)
  statesNumber <- dim(cohorts)[1]
  if (dim(cohorts)[2] >= 1000) {
    dd <- data.frame(t(cohorts))
    dd$cc <- c(rep(1:M, times = floor(dim(cohorts)[2])/M), 
               1:(dim(cohorts)[2]%%M + 1))[1:dim(cohorts)[2]]
    dd <- split(dd, dd$cc)
    prev <- list()
    ppp <- NULL
    for (ciind in 1:M) {
      cohorts <- t(dd[[ciind]][, 1:statesNumber])
      if (ncol(cohorts) == 0) {
        prev[[ciind]] <- matrix(NA, ncol = nrow(cohorts), 
                                nrow = length(times))
      }
      if (ncol(cohorts) == 1) {
        prev[[ciind]] <- matrix(0, ncol = nrow(cohorts), 
                                nrow = length(times))
        for (tInd in 1:(nrow(cohorts) - 1)) {
          prev[[ciind]][times < min(c(max(times) + 1, 
                                      cohorts[(tInd + 1):nrow(cohorts)]), na.rm = TRUE) & 
                          times >= cohorts[tInd], tInd] <- 1
        }
        prev[[ciind]][times >= cohorts[nrow(cohorts)], 
                      nrow(cohorts)] <- 1
      }
      if (ncol(cohorts) >= 2) {
        prep1 <- gems:::data_prep(cohorts, max(times))
        prep1 <- prep1[!(prep1$Tstart == prep1$Tstop), 
        ]
        prepData <- gems:::msmDataPrep(prep1)
        prev[[ciind]] <- prevalence(prepData, times, 
                                    dim(cohorts)[1])
      }
      ppp <- rbind(ppp, prev[[ciind]])
    }
    ddd <- data.frame(ppp)
    ddd$times <- times
    prev <- as.matrix(ddply(ddd, .(times), function(x) colMeans(x)))[, 
                                                                     1:statesNumber]
    lower <- as.matrix(ddply(ddd, .(times), function(x) apply(x, 
                                                              2, function(y) {
                                                                stats::quantile(y, 0.025)
                                                              })))[, 1:statesNumber]
    upper <- as.matrix(ddply(ddd, .(times), function(x) apply(x, 
                                                              2, function(y) {
                                                                stats::quantile(y, 0.975)
                                                              })))[, 1:statesNumber]
  }
  else {
    if (ncol(cohorts) == 0) {
      prev <- matrix(NA, ncol = nrow(cohorts), nrow = length(times))
    }
    if (ncol(cohorts) == 1) {
      prev <- matrix(0, ncol = nrow(cohorts), nrow = length(times))
      for (tInd in 1:(nrow(cohorts) - 1)) {
        prev[times < min(c(max(times) + 1, cohorts[(tInd + 
                                                      1):nrow(cohorts)]), na.rm = TRUE) & times >= 
               cohorts[tInd], tInd] <- 1
      }
      prev[times >= cohorts[nrow(cohorts)], nrow(cohorts)] <- 1
    }
    if (ncol(cohorts) >= 2) {
      prep1 <- gems:::data_prep(cohorts, max(times))
      prep1 <- prep1[!(prep1$Tstart == prep1$Tstop), ]
      prepData <- gems:::msmDataPrep(prep1)
      prev <- prevalence(prepData, times, dim(cohorts)[1])
    }
    lower <- matrix(NA, ncol = nrow(cohorts), nrow = length(times))
    upper <- matrix(NA, ncol = nrow(cohorts), nrow = length(times))
  }
  dimnames(prev) <- list(paste("Time", times), stateNames)
  dimnames(lower) <- list(paste("Time", times), stateNames)
  dimnames(upper) <- list(paste("Time", times), stateNames)
  pp = methods::new("PosteriorProbabilities")
  pp@states = stateNames
  pp@times = times
  pp@probabilities = prev
  pp@lower = lower
  pp@upper = upper
  pp@type = "Transition probabilities"
  return(pp)
}
