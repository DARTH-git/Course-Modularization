#### Formatting functions ####
# Run them all before continuing!
# Function for number of axis ticks in ggplot2 graphs
number_ticks <- function(n) {function(limits) pretty(limits, n)} 
# Total population affected by the decision
TotPop <- function(time, prev, incid, disc = 0){
  # Computes total population afected by technology
  #
  # Args:
  #   time:  vector with time points defining technology lifetime
  #   prev:  present prevalence
  #   incid: incidence
  #   disc:  discount factor; deafult = 0.
  #
  # Returns:
  #   tot.pop: total population afected by technology over technology lifetime
  #  
  # Technology Life Time, the last entry of vector `time`
  LT            <- time[length(time)]
  # Vector with population afected by the technolgy at each time point
  pop.time      <- c(prev, rep(incid, (length(time)-1))) 
  # Vector with present value of population afected by the technolgy at each time point
  disc.pop.time <- pop.time/(1+disc)^time
  # Total population afected by the technology
  tot.pop <-sum(disc.pop.time)
}
# Cost of Research
CostRes <- function(fixed.cost = 0, 
                    samp.size, 
                    cost.per.patient, 
                    INMB, 
                    clin.trial = TRUE, n.arms = 2){
  # Computes the cost of collecting information (i.e., through a research study)
  #
  # Args:
  #   fixed.cost:       fixed cost of collecting information
  #                     (e.g., fixed cost of a clinical trial); default = 0
  #   samp.size:               vector with sample sizes
  #   cost.per.patient: cost per patient in research study
  #   INMB:             Incremental Net Monetary Benefit
  #   clin.trial:       indicator whether calculation is for a clinical trial;
  #                     default = TRUE
  #   n.arms:           Number of arms in research study design; default = 2
  #
  # Returns:
  #   cost.res: vector with the total cost of collecting information for each simple size
  #
  if (clin.trial){
    Cost.Res <- fixed.cost + n.arms*samp.size*cost.per.patient + samp.size*INMB
  } else { # E.g., cohort study
    Cost.Res <- fixed.cost + samp.size*cost.per.patient
  }
  return(Cost.Res)
}