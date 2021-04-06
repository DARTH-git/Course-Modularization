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
    # create the Markov trace matrix M capturing the proportion of the cohort 
    # in each state at each cycle
    m_M_notrt <- m_M_trt <- matrix(NA, 
                                   nrow     = n_t + 1, ncol = n_states,
                                   dimnames = list(paste("cycle", 0:n_t, sep = " "), v_names_states))
    
    # The cohort starts as healthy
    m_M_notrt[1, ] <- m_M_trt[1, ] <- c(1, 0, 0, 0) # initiate first cycle of cohort trace 
    
    # create transition probability matrix for NO treatment
    m_P_notrt  <- matrix(0,
                         nrow = n_states,
                         ncol = n_states,
                         dimnames = list(v_names_states, v_names_states)) # name the columns and rows of the matrix

    # fill in the transition probability matrix
    # from Healthy
    m_P_notrt["H", "H"  ] <- (1 - p_HD) * (1 - p_HS1)
    m_P_notrt["H", "S1" ] <- (1 - p_HD) * p_HS1
    m_P_notrt["H", "D"  ] <- p_HD
    # from Sick
    m_P_notrt["S1", "H" ] <- (1 - p_S1D) * p_S1H
    m_P_notrt["S1", "S1"] <- (1 - p_S1D) * (1 - (p_S1H + p_S1S2))
    m_P_notrt["S1", "S2"] <- (1 - p_S1D) * p_S1S2
    m_P_notrt["S1", "D" ] <- p_S1D
    # from Sicker
    m_P_notrt["S2", "S2"] <- 1 - p_S2D
    m_P_notrt["S2", "D" ] <- p_S2D
    # from Dead
    m_P_notrt["D", "D"  ] <- 1
    
    # Check that transition probabilities are in [0, 1]
    check_transition_probability(m_P_notrt, verbose = TRUE)
    # Check that all rows sum to 1
    check_sum_of_transition_array(m_P_notrt, n_states = n_states, n_cycles = n_t, verbose = TRUE)
    
    # create transition probability matrix for treatment same as no treatment
    m_P_trt <- m_P_notrt
    
    ############# PROCESS ###########################################
    
    for (t in 1:n_t){     # loop through the number of cycles
      m_M_notrt[t + 1, ] <- t(m_M_notrt[t, ]) %*% m_P_notrt  # estimate the Markov trace 
      # for the next cycle (t + 1)
      m_M_trt[t + 1, ]   <- t(m_M_trt[t, ])   %*% m_P_trt    # estimate the Markov trace 
      # for the next cycle (t + 1)
    } # close the loop
    
    ####### RETURN OUTPUT  ###########################################
    out <- list(m_M = m_M_notrt,
                m_P = m_P_notrt)
    
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
    v_dwc <- 1 / ((1 + d_e) ^ (0:(n_t))) # vector with discount weights for costs
    v_dwe <- 1 / ((1 + d_c) ^ (0:(n_t))) # vector with discount weights for QALYs
    
    ## Run STM model at a parameter set for each intervention
    l_model_out_no_trt <- decision_model(l_params_all = l_params_all)
    l_model_out_trt    <- decision_model(l_params_all = l_params_all)
    
    ## Cohort trace by treatment
    m_M_no_trt  <- l_model_out_no_trt$m_M # No treatment
    m_M_trt     <- l_model_out_trt$m_M    # Treatment
    
    # Vectors with costs and utilities by treatment
    v_u_notrt   <- c(u_H, u_S1,  u_S2, u_D)
    v_u_trt     <- c(u_H, u_trt, u_S2, u_D)
    
    v_c_notrt   <- c(c_H, c_S1, c_S2, c_D)
    v_c_trt     <- c(c_H, c_S1 + c_trt, c_S2 + c_trt, c_D)
    
    ## Mean Costs and QALYs for Treatment and NO Treatment
    v_tu_notrt  <- m_M_notrt   %*%  v_u_notrt
    v_tu_trt    <- m_M_trt     %*%  v_u_trt
    
    v_tc_notrt  <- m_M_notrt   %*%  v_c_notrt
    v_tc_trt    <- m_M_trt     %*%  v_c_trt 
    
    ## Total discounted mean Costs and QALYs
    tu_d_notrt  <- t(v_tu_notrt)   %*%  v_dwe   
    tu_d_trt    <- t(v_tu_trt)     %*%  v_dwe
    
    tc_d_notrt  <- t(v_tc_notrt)   %*%  v_dwc
    tc_d_trt    <- t(v_tc_trt)     %*%  v_dwc
    
    # store them into a vector
    v_tc_d      <- c(tc_d_notrt, tc_d_trt)
    v_tu_d      <- c(tu_d_notrt, tu_d_trt)
    
    ## Vector with discounted net monetary benefits (NMB)
    v_nmb_d     <- v_tu_d * n_wtp - v_tc_d
    
    ## Dataframe with discounted costs, effectiveness and NMB
    df_ce <- data.frame(Strategy = v_names_str,
                        Cost     = v_tc_d,
                        Effect   = v_tu_d,
                        NMB      = v_nmb_d)
    
    return(df_ce)
  }
  )
}
