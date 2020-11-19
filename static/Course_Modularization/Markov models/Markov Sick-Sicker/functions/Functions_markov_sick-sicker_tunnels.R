#-----------------------------------------#
####          Decision Model           ####
#-----------------------------------------#
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
    r_HD    = - log(1 - p_HD)          # rate of death in healthy
    r_S1D   = hr_S1 * r_HD 	           # rate of death in sick
    r_S2D   = hr_S2 * r_HD  	         # rate of death in sicker
    p_S1D   = 1 - exp(-r_S1D)          # probability to die in sick
    p_S2D   = 1 - exp(-r_S2D)          # probability to die in sicker
    
    ########################################## INITIALIZATION ##########################################
    ### create the markov trace matrix M capturing the proportion of the cohort in each state at each cycle
    m_M_no_trt <- m_M_trt <- matrix(NA, 
                                    nrow = n_t + 1, ncol = n_s_td, # Update!
                                    dimnames = list(paste("cycle", 0:n_t, sep = " "), v_n_td)) # Update!
    
    head(m_M_no_trt) # show first 6 rows of the matrix 
    
    # The cohort starts as healthy
    m_M_no_trt[1, ] <- m_M_trt[1, ] <- c(1, rep(0, tunnel_size), 0, 0) # initialize first cycle of Markov trace accounting for the tunnels
    
    ### create transition probability array for NO treatment
    a_P_no_trt <- array(0,                                          # Create 3-D array
                       dim = c(n_s_td, n_s_td, n_t),
                       dimnames = list(v_n_td, v_n_td, 0:(n_t-1)))  # name dimensions of the transition probability array
    
    ### fill in the transition probability array
    # From Healthy
    a_P_no_trt["H", "H", ]            <- 1 - (p_HS1 + p_HD)
    a_P_no_trt["H", Sick_tunnel[1], ] <- p_HS1
    a_P_no_trt["H", "D", ]            <- p_HD
    
    # From Sick
    for(i in 1:(tunnel_size - 1)){
      a_P_no_trt[Sick_tunnel[i], "H", ]  <- p_S1H
      a_P_no_trt[Sick_tunnel[i], Sick_tunnel[i + 1], ] <- 1 - (p_S1H + p_S1S2[i] + p_S1D)
      a_P_no_trt[Sick_tunnel[i], "S2", ] <- p_S1S2[i]
      a_P_no_trt[Sick_tunnel[i], "D", ]  <- p_S1D
    }
    a_P_no_trt[Sick_tunnel[tunnel_size], "H", ]  <- p_S1H
    a_P_no_trt[Sick_tunnel[tunnel_size], Sick_tunnel[tunnel_size], ] <- 1 - (p_S1H + p_S1S2[tunnel_size] + p_S1D)
    a_P_no_trt[Sick_tunnel[tunnel_size], "S2", ] <- p_S1S2[tunnel_size]
    a_P_no_trt[Sick_tunnel[tunnel_size], "D", ]  <- p_S1D
    
    # From Sicker
    a_P_no_trt["S2", "S2", ] <- 1 - p_S2D
    a_P_no_trt["S2", "D", ]  <- p_S2D
    
    # From Dead
    a_P_no_trt["D", "D", ]   <- 1
    
    ### Check if transition matrix is valid (i_e_, each row should add up to 1)
    valid <- apply(a_P_no_trt, 3, function(x) sum(rowSums(x))==n_s_td)
    if (!isTRUE(all.equal(as.numeric(sum(valid)), as.numeric(n_t)))) {
      stop("This is not a valid transition Matrix")
    }
    
    ### create transition probability matrix for treatment same as NO treatment
    a_P_trt <- a_P_no_trt
    
    ########################################## PROCESS ##########################################
    
    for (t in 1:n_t){                                                   # loop through the number of cycles
      m_M_no_trt[t + 1, ] <- t(m_M_no_trt[t, ]) %*% a_P_no_trt[ , , t]  # estimate the Markov trace for cycle the next cycle (t + 1)
      m_M_trt[t + 1, ]    <- t(m_M_trt[t, ])    %*% a_P_trt[, , t]      # estimate the Markov trace for cycle the next cycle (t + 1)
    } # lose the loop
    
    ### create aggregated traces
    m_M_td_no_trt <- cbind(H  = m_M_no_trt[, "H"], 
                           S1 = rowSums(m_M_no_trt[, 2:(tunnel_size +1)]), 
                           S2 = m_M_no_trt[, "S2"],
                           D  = m_M_no_trt[, "D"])
    head(m_M_td_no_trt)
    m_M_td_trt    <- cbind(H  = m_M_trt[, "H"], 
                           S1 = rowSums(m_M_trt[, 2:(tunnel_size +1)]), 
                           S2 = m_M_trt[, "S2"],
                           D  = m_M_trt[, "D"])
    
    ########################################## EPIDEMIOLOGICAL OUTPUT  ##########################################
    #### Overall Survival (OS) ####
    v_os_no_trt_tunnels <- 1 - m_M_no_trt[, "D"]       # calculate the overall survival (OS) probability for no treatment
    v_os_no_trt_tunnels <- rowSums(m_M_no_trt[, 1:3])  # alternative way of calculating the OS probability   
    
    #### Life Expectancy (LE) ####
    v_le_tunnels <- sum(v_os_no_trt_tunnels)           # summing probablity of OS over time  (i_e_ life expectancy)
    
    #### Disease prevalence ####
    v_prev_tunnels <- rowSums(m_M_td_no_trt[, c("S1", "S2")]) / v_os_no_trt_tunnels
    
    #### ratio of sick(S1) vs sicker(S2) ####
    v_ratio_S1S2_tunnels <- m_M_td_no_trt[, "S1"] / m_M_td_no_trt[, "S2"]
    
    ########################################## RETURN OUTPUT  ##########################################
    out <- list(m_M_td_no_trt = m_M_td_no_trt,
                m_M_td_trt    = m_M_td_trt,
                a_P_no_trt    = a_P_no_trt,
                a_P_trt       = a_P_trt,
                Surv          = v_os_no_trt_tunnels[-1],
                Prev          = v_prev_tunnels[-1],
                Ratio_S1S2    = v_ratio_S1S2_tunnels)
    
    return(out)
  }
  )
}

#---------------------------------------------#
#### Calculate cost-effectiveness outcomes ####
#---------------------------------------------#
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
    m_M_no_trt  <- l_model_out_no_trt$m_M_td_no_trt # No treatment
    m_M_trt     <- l_model_out_trt$m_M_td_trt    # Treatment
    
    ## Vectors with costs and utilities by treatment
    v_u_no_trt  <- c(u_H, u_S1, u_S2, u_D)
    v_u_trt     <- c(u_H, u_trt, u_S2, u_D)
    
    v_c_no_trt  <- c(c_H, c_S1, c_S2, c_D)
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
