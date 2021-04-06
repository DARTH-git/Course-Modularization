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
    ## Initial state vector
    # All starting healthy
    v_s_init <- c("Healthy" = 1, "Sick" = 0, "Dead" = 0)  
    v_s_init
    
    ## Initialize cohort trace for cSTM for all strategies
    m_M_SoC <- matrix(0, 
                      nrow = (n_cycles + 1), ncol = n_states, 
                      dimnames = list(v_names_cycles, v_names_states))
    # Store the initial state vector in the first row of the cohort trace
    m_M_SoC[1, ] <- v_s_init
    ## Initialize cohort traces
    m_M_trtA <- m_M_trtB <- m_M_SoC # structure and initial states remain the same
    
    ## Initialize transition probability matrix 
    # all transitions to a non-death state are assumed to be conditional on survival 
    m_P_SoC  <- matrix(0,
                       nrow = n_states, ncol = n_states,
                       dimnames = list(v_names_states, v_names_states)) # define row and column names
    
    ## Standard of Care
    # from Healthy
    m_P_SoC["Healthy", "Healthy"] <- (1 - p_HD) * (1 - p_HS_SoC)
    m_P_SoC["Healthy", "Sick"]    <- (1 - p_HD) *      p_HS_SoC
    m_P_SoC["Healthy", "Dead"]    <-      p_HD
    
    # from Sick
    m_P_SoC["Sick", "Sick"] <- 1 - p_SD
    m_P_SoC["Sick", "Dead"] <-     p_SD
    
    # from Dead
    m_P_SoC["Dead", "Dead"] <- 1
    
    ## Treatment A
    m_P_trtA <- m_P_SoC
    m_P_trtA["Healthy", "Healthy"] <- (1 - p_HD) * (1 - p_HS_trtA)
    m_P_trtA["Healthy", "Sick"]    <- (1 - p_HD) *      p_HS_trtA
    
    ## Treatment B
    m_P_trtB <- m_P_SoC
    m_P_trtB["Healthy", "Healthy"] <- (1 - p_HD) * (1 - p_HS_trtB)
    m_P_trtB["Healthy", "Sick"]    <- (1 - p_HD) *      p_HS_trtB
    
    # Check that transition probabilities are in [0, 1]
    check_transition_probability(m_P_SoC,  verbose = TRUE)
    check_transition_probability(m_P_trtA, verbose = TRUE)
    check_transition_probability(m_P_trtB, verbose = TRUE)
    # Check that all rows sum to 1
    check_sum_of_transition_array(m_P_SoC,  n_states = n_states, verbose = TRUE)
    check_sum_of_transition_array(m_P_trtA, n_states = n_states, verbose = TRUE)
    check_sum_of_transition_array(m_P_trtB, n_states = n_states, verbose = TRUE)

    ############# PROCESS ###########################################
    for (t in 1:n_cycles){  # loop through the number of cycles
      m_M_SoC [t + 1, ] <- m_M_SoC [t, ] %*% m_P_SoC   # estimate the state vector for the next cycle (t + 1)
      m_M_trtA[t + 1, ] <- m_M_trtA[t, ] %*% m_P_trtA  # estimate the state vector for the next cycle (t + 1)
      m_M_trtB[t + 1, ] <- m_M_trtB[t, ] %*% m_P_trtB  # estimate the state vector for the next cycle (t + 1)
    }
    
    ####### RETURN OUTPUT  ###########################################
    out <- list(m_M_SoC   = m_M_SoC,
                m_M_trtA  = m_M_trtA,
                m_M_trtB  = m_M_trtB,
                m_P_SoC   = m_P_SoC,
                m_P_trtA  = m_P_trtA,
                m_P_trtB  = m_P_trtB)
    
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
#' @param n_wtp Willingness-to-pay threshold to compute net monetary benefits (
#' NMB)
#' @return A dataframe with discounted costs, effectiveness and NMB.
#' 
calculate_ce_out <- function(l_params_all, n_wtp = 10000){ # User defined
  with(as.list(l_params_all), {
    # discount weights for costs and effects
    v_dwc <- 1 / (1 + d_c) ^ (0:n_cycles) 
    v_dwe <- 1 / (1 + d_e) ^ (0:n_cycles) 
    
    ## Run STM model at a parameter set 
    l_model_out <- decision_model(l_params_all = l_params_all)
    
    ## Cohort traces 
    m_M_SoC  <- l_model_out$m_M_SoC 
    m_M_trtA <- l_model_out$m_M_trtA
    m_M_trtB <- l_model_out$m_M_trtB
    
    # per cycle
    # calculate expected costs by multiplying cohort trace with the cost vector for the different health states   
    v_tc_SoC  <- m_M_SoC  %*% c(c_H, c_S, c_D)  
    v_tc_trtA <- m_M_trtA %*% c(c_H + c_trtA, c_S, c_D)  
    v_tc_trtB <- m_M_trtB %*% c(c_H + c_trtB, c_S, c_D)  
    
    # calculate expected QALYs by multiplying cohort trace with the utilities for the different health states   
    v_tu_SoC  <- m_M_SoC  %*% c(u_H, u_S, u_D)  
    v_tu_trtA <- m_M_trtA %*% c(u_H, u_S, u_D) 
    v_tu_trtB <- m_M_trtB %*% c(u_H, u_S, u_D) 
    
    # Discount costs by multiplying the cost vector with discount weights (v_dw) 
    tc_d_SoC  <-  t(v_tc_SoC)  %*% v_dwc
    tc_d_trtA <-  t(v_tc_trtA) %*% v_dwc
    tc_d_trtB <-  t(v_tc_trtB) %*% v_dwc
    
    # Discount QALYS by multiplying the QALYs vector with discount weights (v_dw)
    tu_d_SoC  <-  t(v_tu_SoC)  %*% v_dwe
    tu_d_trtA <-  t(v_tu_trtA) %*% v_dwe
    tu_d_trtB <-  t(v_tu_trtB) %*% v_dwe
    
    # Store them into a vector
    v_tc_d <- c(tc_d_SoC, tc_d_trtA, tc_d_trtB)
    v_tu_d <- c(tu_d_SoC, tu_d_trtA, tu_d_trtB)
    
    # Vector with discounted net monetary benefits (NMB)
    v_nmb_d   <- v_tu_d * n_wtp - v_tc_d
    
    ## Dataframe with discounted costs, effectiveness and NMB
    df_ce <- data.frame(Strategy = v_names_str,
                        Cost     = v_tc_d,
                        Effect   = v_tu_d,
                        NMB      = v_nmb_d)
    
    return(df_ce)
  }
  )
}
