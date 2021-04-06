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
    r_HD    = - log(1 - p_HD) # rate of death in healthy
    r_S1D   = hr_S1 * r_HD 	  # rate of death in sick
    r_S2D   = hr_S2 * r_HD  	# rate of death in sicker
    p_S1D   = 1 - exp(-r_S1D) # probability to die in sick
    p_S2D   = 1 - exp(-r_S2D) # probability to die in sicker
    
    ####### INITIALIZATION ##########################################
    # create the cohort trace
    m_M <- matrix(NA, nrow = n_t + 1 , 
                  ncol = n_s,
                  dimnames = list(0:n_t, v_names_states))  # create Markov trace (n_t + 1 because R doesn't understand  Cycle 0)
    
    m_M[1, ] <- c(1, 0, 0, 0)  # initialize Markov trace
    
    # create transition probability matrix for NO treatment
    m_P <- matrix(0,
                  nrow = n_s, 
                  ncol = n_s,
                  dimnames = list(v_names_states, v_names_states))
    # fill in the transition probability array
    ### From Healthy
    m_P["H", "H"]   <- (1 - p_HD) * (1 - p_HS1)
    m_P["H", "S1"]  <- (1 - p_HD) * p_HS1
    m_P["H", "D"]   <- p_HD
    ### From Sick
    m_P["S1", "H"]  <- (1- p_S1D)  * p_S1H
    m_P["S1", "S1"] <- (1 - p_S1D) * (1 - (p_S1H + p_S1S2))
    m_P["S1", "S2"] <- (1 - p_S1D) * p_S1S2
    m_P["S1", "D"]  <- p_S1D
    ### From Sicker
    m_P["S2", "S2"] <- 1 - p_S2D
    m_P["S2", "D"]  <- p_S2D
    ### From Dead
    m_P["D", "D"]   <- 1
  
    ############# PROCESS ###########################################
    
    for (t in 1:n_t){                     # throughout the number of cycles
      m_M[t + 1, ] <- m_M[t, ] %*% m_P    # estimate the Markov trace for cycle the next cycle (t + 1)
    }
    
    ####### EPIDEMIOLOGICAL OUTPUT  ###########################################
    #### Overall Survival (OS) ####
    v_os      <- 1 - m_M[, "D"]           # calculate the overall survival (OS) probability for no treatment
    
    #### Disease prevalence #####
    v_prev    <- rowSums(m_M[, c("S1", "S2")])/v_os
    
    #### Proportion of sick in S1 state #####
    v_prop_S1 <- m_M[, "S1"] / v_prev
    
    ####### RETURN OUTPUT  ###########################################
    out <- list(m_M      = m_M,
                m_P      = m_P,
                Surv     = v_os[-1],
                Prev     = v_prev[-1],
                PropSick = v_prop_S1[c(11, 21, 31)])
    
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
    
    ## Vectors with costs and utilities by treatment
    v_u_no_trt  <- c(u_H, u_S1,  u_S2, u_D)
    v_u_trt     <- c(u_H, u_trt, u_S2, u_D)
    
    v_c_no_trt  <- c(c_H, c_S1,         c_S2,         c_D)
    v_c_trt     <- c(c_H, c_S1 + c_trt, c_S2 + c_trt, c_D)
    
    ## Mean Costs and QALYs for Treatment and NO Treatment
    v_tu_no_trt <- m_M_no_trt %*% v_u_no_trt
    v_tu_trt    <- m_M_trt %*% v_u_trt
    
    v_tc_no_trt <- m_M_no_trt %*% v_c_no_trt
    v_tc_trt    <- m_M_trt %*% v_c_trt
    
    ## Total discounted mean Costs and QALYs
    tu_d_no_trt <- t(v_tu_no_trt) %*% v_dwe 
    tu_d_trt    <- t(v_tu_trt) %*% v_dwe
    
    tc_d_no_trt <- t(v_tc_no_trt) %*% v_dwc
    tc_d_trt    <- t(v_tc_trt)    %*% v_dwc
    
    ## Vector with total discounted mean Costs and QALYs
    v_tc_d      <- c(tc_d_no_trt, tc_d_trt)
    v_tu_d      <- c(tu_d_no_trt, tu_d_trt)
    
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
