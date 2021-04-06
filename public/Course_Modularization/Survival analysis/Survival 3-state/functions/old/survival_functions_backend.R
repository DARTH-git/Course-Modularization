### These functions are part of the msm package
### Will need these to speed up Alan's `prevalence.msm1` function

observed.msm <- function(x, times=NULL, interp=c("start","midpoint"), censtime=Inf, subset=NULL)
{
  if (!inherits(x, "msm")) stop("expected x to be a msm model")
  ## For general HMMs use the Viterbi estimate of the observed state.
  if (!is.null(x$pci)) {
    state <- x$data$mf$"(state)"[!x$data$mf$"(pci.imp)"]
    time <- x$data$mf$"(time)"[!x$data$mf$"(pci.imp)"]
    subject <- x$data$mf$"(subject)"
    subject <- subject[!x$data$mf$"(pci.imp)"]
  } else {
    if ((x$hmodel$hidden && !x$emodel$misc) ||
        (!x$emodel$misc && x$cmodel$ncens>0) )
      state <- viterbi.msm(x)$fitted
    else if (x$emodel$misc && x$cmodel$ncens>0) {
      vit <- viterbi.msm(x)$fitted
      state <- x$data$mf$"(state)"
      state[state %in% x$cmodel$censor] <- vit[state %in% x$cmodel$censor]
      ## TODO for misc models with censoring, impute only censored obs states from viterbi
    }  else
      state <- x$data$mf$"(state)"
    
    time <- x$data$mf$"(time)"; subject <- x$data$mf$"(subject)"
  }
  if (is.null(subset)) subset <- unique(subject)
  time <- time[subject %in% subset]
  state <- state[subject %in% subset]
  subject <- subject[subject %in% subset]
  if (is.null(times))
    times <- seq(min(time), max(time), (max(time) - min(time))/10)
  states.expand <- matrix(nrow=length(unique(subject)), ncol=length(times))
  pts <- unique(subject)
  absorb <- absorbing.msm(x)
  interp <- match.arg(interp)
  if (!is.numeric(censtime)) stop("censtime should be numeric")
  if (length(censtime)==1) censtime <- rep(censtime, length(pts))
  else if (length(censtime)!=length(pts)) stop("censtime of length ", length(censtime), ", should be 1 or ", length(pts))
  for (i in seq_along(pts)){
    state.i <- state[(subject==pts[i])]
    time.i <- time[(subject==pts[i])]
    j <- 1
    while(j <= length(times)) {
      if (times[j] < time.i[1]) {
        mtime <- max(which(times-time.i[1] < 0))
        states.expand[i, j:mtime] <- NA
        j <- mtime + 1
        next;
      } else if (times[j] > time.i[length(time.i)]) {
        if (state.i[length(time.i)] %in% absorb && (times[j] <= censtime[i])) {
          states.expand[i, j:(length(times))] <-  state.i[length(time.i)]
        } else states.expand[i, j:(length(times))] <-  NA
        break;
      } else {
        prevtime.ind <- max(which(time.i <= times[j]))
        prevtime <- time.i[prevtime.ind]
        if (interp=="midpoint") {
          nexttime.ind <- min(which(time.i >= times[j]))
          nexttime <- time.i[nexttime.ind]
          midpoint <- (prevtime + nexttime) / 2
          states.expand[i,j] <- state.i[if (times[j] <= midpoint) prevtime.ind else nexttime.ind]
        } else
          states.expand[i,j] <- state.i[prevtime.ind]
      }
      j <- j+1
    }
  }
  obstab <- t(apply(states.expand, 2, function(y) table(factor(y, levels=seq(length=x$qmodel$nstates)))))
  obsperc <- 100*obstab / rep(rowSums(obstab), ncol(obstab))
  dimnames(obstab) <- dimnames(obsperc) <- list(times, paste("State", 1:x$qmodel$nstates))
  obstab <- cbind(obstab, Total=rowSums(obstab))
  
  covhist <- get.covhist(x, subset)
  covcat <- ## distinct covariate history group each subject falls into (ordinal)
    if (is.null(covhist)) rep(1, length(unique(subject)))
  else match(covhist$hist, unique(covhist$hist))
  risk <- matrix(nrow=length(times), ncol=length(unique(covcat)), dimnames = list(times, unique(covhist$hist)))
  for (i in seq_along(unique(covcat))) {
    obst <- t(apply(states.expand[covcat==unique(covcat)[i],,drop=FALSE], 2,
                    function(y) table(factor(y, levels=seq(length=x$qmodel$nstates)))))
    risk[,i] <- rowSums(obst)
  }
  
  list(obstab=obstab, obsperc=obsperc, risk=risk)
}


get.covhist <- function(x, subset=NULL) {
  ## Keep only times where the covariate changes, or first or last obs
  mf <- x$data$mf
  if (x$qcmodel$ncovs > 0) {
    if (!is.null(subset)) {
      subs <- mf$"(subject)" %in% subset
      mf <- mf[subs,,drop=FALSE]
    }
    subj <- match(mf$"(subject)", unique(mf$"(subject)"))
    n <- length(subj)
    apaste <- do.call("paste", mf[,attr(mf,"covnames"),drop=FALSE])
    first <- !duplicated(subj); last <- rev(!duplicated(rev(subj)))
    keep <- (c(0, apaste[1:(n-1)]) != apaste) | first | last
    ## Keep and tabulate unique covariate series
    covseries <- split(apaste[keep], subj[keep]) # as a list of char vectors
    covseries <- sapply(covseries, paste, collapse=" , ") # as one char vector, one series per pt.
    ## also need p matrices for different times as well as different covs.
    ## but only interested in cov change times if there's more than one
    ## transition (at least one times change point)
    change.times <- mf$"(time)"; change.times[first] <- change.times[last] <- 0
    change.times <- split(change.times[keep & (!(first|last))], subj[keep & (!(first|last))])
    change.times <- sapply(change.times, paste, collapse= " , ")
    covseries.t <- paste(covseries, change.times, sep="; ")
    ids <- unique(subj)[!duplicated(covseries.t)] # subj ids, one with each distinct series
    ncombs <- table(covseries.t)[unique(covseries.t)]# how many per series
    covmat <- cbind(subject=subj, time=mf$"(time)", mf[,attr(mf,"covnames"),drop=FALSE])
    covmat <- covmat[(subj %in% ids) & keep,]
    list(example=covmat, # rows of the original data sufficient to define the distinct histories
         hist=covseries.t) # one per subject listing their covariate history as a string
  }
  else NULL
}


### Alan will try to speed up `prevalence.msm` function in the `msm` package
prevalence.msm1 <- function(x,
                            times=NULL,
                            timezero=NULL,
                            initstates=NULL,
                            covariates="population",
                            misccovariates="mean",
                            piecewise.times=NULL,
                            piecewise.covariates=NULL,
                            ci=c("none","normal","bootstrap"),
                            cl = 0.95,
                            B = 1000,
                            cores = NULL,
                            interp=c("start","midpoint"),
                            censtime=Inf,
                            subset=NULL,
                            plot = FALSE, ...
)
{
  if (!inherits(x, "msm")) stop("expected x to be a msm model")
  ## Estimate observed state occupancies in the data at a series of times
  time <- model.extract(x$data$mf, "time")
  if (is.null(times))
    times <- seq(min(time), max(time), (max(time) - min(time))/10)
  obs <- observed.msm(x, times, interp, censtime, subset)
  res <- list(obsperc=obs$obsperc)
  names(res) <- c("Observed percentages")
  res
}


trace.DES1 = function(msm_sim = des_sim, tmat, n_i, times )
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

  M.tr.des1 <- prevalence.msm1(fit.msm.sim, times = times) # Markov trace when DES model is used
  
  return(M.tr.des1[[1]]/100)
}

