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
    
    ####### INITIALIZATION ##########################################
    # create the cohort trace
    m_M <- m_M_trt <-  matrix(NA, 
                              nrow = n_t + 1 ,  # create Markov trace (n.t + 1 because R doesn't 
                              # understand Cycle 0)
                              ncol = n_states, 
                              dimnames = list(0:n_t, v_n))
    
    m_M[1, ] <- m_M_trt[1, ] <- v_init          # initialize first cycle of Markov trace
    
    # create the transition probability matrix
    m_P  <- matrix(0,
                   nrow = n_states, ncol = n_states,
                   dimnames = list(v_n, v_n))   # name the columns and rows of the transition 
                                                # probability matrix

    # fill in the transition probability matrix
    # from Healthy
    m_P["Healthy", "Healthy"] <- (1 - p_HD) * (1 - p_HS)
    m_P["Healthy", "Sick"]    <- (1 - p_HD) * p_HS
    m_P["Healthy", "Dead"]    <- p_HD
    
    # from Sick
    m_P["Sick", "Sick"] <- 1 - p_SD
    m_P["Sick", "Dead"] <- p_SD
    
    # from Dead
    m_P["Dead", "Dead"] <- 1
    
    # Under treatment
    m_P_trt <- m_P
    m_P_trt["Healthy", "Healthy"] <- (1 - p_HD) * (1 - p_HS_trt)
    m_P_trt["Healthy", "Sick"]    <- (1 - p_HD) * p_HS_trt
    
    # Check that transition probabilities are in [0, 1]
    check_transition_probability(m_P, verbose = TRUE)
    check_transition_probability(m_P_trt, verbose = TRUE)
    # Check that all rows sum to 1
    check_sum_of_transition_array(m_P, n_states = n_states, n_cycles = n_t, verbose = TRUE)
    check_transition_probability(m_P_trt, verbose = TRUE)
    
    ############# PROCESS ###########################################
    
    for (t in 1:n_t){                               # loop through the number of cycles
      m_M[t + 1, ]     <- m_M[t, ]     %*% m_P      # estimate the state vector for the next cycle (t + 1)
      m_M_trt[t + 1, ] <- m_M_trt[t, ] %*% m_P_trt  # for treatment
    }
    
    ####### RETURN OUTPUT  ###########################################
    out <- list(m_M      = m_M,
                m_M_trt  = m_M_trt,
                m_P      = m_P,
                m_P_trt  = m_P_trt)
    
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
calculate_ce_out <- function(l_params_all, n_wtp = 10000){ # User defined
  with(as.list(l_params_all), {
    ## Create discounting vectors
    v_dwc <- 1 / ((1 + d_e) ^ (0:(n_t))) # vector with discount weights for costs
    v_dwe <- 1 / ((1 + d_c) ^ (0:(n_t))) # vector with discount weights for QALYs
    
    ## Run STM model at a parameter set 
    l_model_out <- decision_model(l_params_all = l_params_all)
    
    ## Cohort trace 
    m_M     <- l_model_out$m_M 
    m_M_trt <- l_model_out$m_M_trt
    
    # per cycle
    # calculate expected costs by multiplying m_M with the cost vector for the different 
    # health states   
    v_tc     <- m_M     %*% c(c_H, c_S, c_D)          # Standard of Care
    v_tc_trt <- m_M_trt %*% c(c_H, c_S + c_trt, c_D)  # Treatment
    # calculate expected QALYs  by multiplying m_M with the utilities for the different 
    # health states   
    v_tu     <- m_M     %*% c(u_H, u_S, u_D)          # Standard of Care
    v_tu_trt <- m_M_trt %*% c(u_H, u_S, u_D)          # Treatment
    
    # Discount costs by multiplying the cost vector with discount weights  
    tc_d     <-  t(v_tc)     %*% v_dwc      # Standard of Care
    tc_d_trt <-  t(v_tc_trt) %*% v_dwc      # Treatment
    # Discount QALYS by multiplying the QALYs vector with discount weights 
    tu_d     <-  t(v_tu)     %*% v_dwe      # Standard of Care
    tu_d_trt <-  t(v_tu_trt) %*% v_dwe      # Treatment
    
    # store them into a vector
    v_tc_d   <- c(tc_d, tc_d_trt)
    v_tu_d   <- c(tu_d, tu_d_trt)
    
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
