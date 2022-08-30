#------------------------------------------------------------------------------#
####                         Decision Model                                 ####
#------------------------------------------------------------------------------#
#' Decision Model
#'
#' \code{decision_model} implements the decision model used.
#'
#' @param l_params_all List with all parameters of decision model
#' @param verbose Logical variable to indicate print out of messages
#' @return The transition probability array and the cohort trace matrix.
#' 
decision_model <- function(l_params_all, verbose = FALSE) {
  with(as.list(l_params_all), {
    # compute internal paramters as a function of external parameters
    r_HD    = - log(1 - p_HD)  # rate of death in healthy
    r_S1D   = hr_S1 * r_HD 	   # rate of death in sick
    r_S2D   = hr_S2 * r_HD  	 # rate of death in sicker
    p_S1D   = 1 - exp(-r_S1D)  # probability to die in sick
    p_S2D   = 1 - exp(-r_S2D)  # probability to die in sicker
    
    ####### INITIALIZATION ##########################################
    
    ############# PROCESS ###########################################
    
    ####### RETURN OUTPUT  ###########################################
    out <- 
    
    return(out)
  }
  )
}

#------------------------------------------------------------------------------#
####              Calculate cost-effectiveness outcomes                     ####
#------------------------------------------------------------------------------#
#' Calculate cost-effectiveness outcomes
#'
#' \code{calculate_ce_out} calculates costs and effects for a given vector of parameters using a simulation model.
#' @param l_params_all List with all parameters of decision model
#' @param n_wtp Willingness-to-pay threshold to compute net benefits
#' @return A data frame with discounted costs, effectiveness and NMB.
#' 
calculate_ce_out <- function(l_params_all, n_wtp = 100000){ # User defined
  with(as.list(l_params_all), {
    ## Create discounting vectors
    
    ## Run STM model at a parameter set for each intervention
    
    ## Cohort trace by treatment
    
    ## Vectors with costs and utilities by treatment
    
    ## Mean Costs and QALYs for Treatment and NO Treatment
    
    ## Total discounted mean Costs and QALYs
    
    ## Vector with total discounted mean Costs and QALYs
    
    ## Vector with discounted net monetary benefits (NMB)
    
    ## Dataframe with discounted costs, effectiveness and NMB
    df_ce <- data.frame(Strategy = v_names_str,
                        Cost     = v_tc_d,
                        Effect   = v_tu_d,
                        NMB      = v_nmb_d)
    
    return(df_ce)
  }
  )
}
