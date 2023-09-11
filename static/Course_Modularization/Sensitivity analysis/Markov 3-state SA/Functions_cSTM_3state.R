
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
#' @export
calculate_ce_out <- function(l_params_all, n_wtp = 10000 ,verbose= verbose){ # User defined
  with(as.list(l_params_all), {
    ########################### Process model inputs ###########################
    ## Model states
    n_cycles        <- n_time_horizon_yr / cycle_length   # number of cycles
    v_names_cycles  <- paste("cycle", 0:n_cycles)         # cycle names
    v_names_states  <- c("Healthy", "Sick", "Dead")       # state names
    n_states        <- length(v_names_states)             # number of health states 
    
    ### Cycle-specific discount weight for costs and effects 
    v_dwc   <- 1 / ((1 + (d_e * cycle_length)) ^ (0:n_cycles))
    v_dwe   <- 1 / ((1 + (d_c * cycle_length)) ^ (0:n_cycles))
    
    ### Calculate cycle-specific transition probabilities
    p_HD <- 1 - exp(-(-log(1 - p_HD_yr) * cycle_length)) # probability of dying from Healthy
    p_SD <- 1 - exp(-(-log(1 - p_SD_yr) * cycle_length)) # probability of dying from Sick
    
    # probability of becoming sick from Healthy
    p_HS_SoC  <- 1 - exp(-(-log(1 - p_HS_yr_SoC) * cycle_length))  # Standard of Care
    p_HS_trtA <- 1 - exp(-(-log(1 - p_HS_yr_trtA) * cycle_length)) # Treatment A
    p_HS_trtB <- 1 - exp(-(-log(1 - p_HS_yr_trtB) * cycle_length)) # Treatment B
    
    ### Calculate cycle-specific state rewards
    #### Costs
    v_c_SoC  <- c(c_H_yr, c_S_yr, c_D_yr) * cycle_length             # Standard of Care
    v_c_trtA <- c(c_H_yr + c_trtA_yr, c_S_yr, c_D_yr) * cycle_length # Treatment A
    v_c_trtB <- c(c_H_yr + c_trtB_yr, c_S_yr, c_D_yr) * cycle_length # Treatment B
    #### QALYs
    v_q_SoC  <- c(u_H, u_S, u_D) * cycle_length # Standord of Care
    v_q_trtA <- v_q_trtB <- v_q_SoC             # Treatments A and B have same utilities as SoC
    
    
    # All starting healthy
    v_m_init <- c("Healthy" = 1, "Sick" = 0, "Dead" = 0)  
    
    ###################### Construct state-transition models ###################
    ### Initialize cohort trace for SoC 
    m_M_SoC <- matrix(0, 
                      nrow = (n_cycles + 1), ncol = n_states, 
                      dimnames = list(v_names_cycles, v_names_states))
    # Store the initial state vector in the first row of the cohort trace
    m_M_SoC[1, ] <- v_m_init
    
    ## Initialize cohort traces for treatments A and B
    # Structure and initial states are the same as for SoC
    m_M_trtA <- m_M_trtB <- m_M_SoC
    
    ## Create transition probability arrays for strategy SoC 
    ### Initialize transition probability array for strategy SoC 
    # All transitions to a non-death state are assumed to be conditional on survival
    m_P_SoC <- matrix(0,  # Create transition probability matrix 
                     nrow = n_states,
                     ncol = n_states,
                     dimnames = list(v_names_states, v_names_states)) # name the dimensions of the array 
    
    ### Fill in array
    ## Standard of Care
    # from Healthy
    ### Fill in matrix 
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
    # Start with same matrix as SoC, but replace parameters that differ for trtA
    m_P_trtA <- m_P_SoC 
    m_P_trtA["Healthy", "Healthy"] <- (1 - p_HD) * (1 - p_HS_trtA)
    m_P_trtA["Healthy", "Sick"]    <- (1 - p_HD) *      p_HS_trtA
    
    ## Treatment B
    # Start with same matrix as SoC, but replace parameters that differ for trtB
    m_P_trtB <- m_P_SoC
    m_P_trtB["Healthy", "Healthy"] <- (1 - p_HD) * (1 - p_HS_trtB)
    m_P_trtB["Healthy", "Sick"]    <- (1 - p_HD) *      p_HS_trtB
    
    ## Check if transition array and probabilities are valid
    # Check that transition probabilities are in [0, 1]
    check_transition_probability(m_P_SoC,  verbose = verbose)
    check_transition_probability(m_P_trtA, verbose = verbose)
    check_transition_probability(m_P_trtB, verbose = verbose)
    # Check that all rows sum to 1
    check_sum_of_transition_array(m_P_SoC,  n_states = n_states, n_cycles = n_cycles, verbose = verbose)
    check_sum_of_transition_array(m_P_trtA, n_states = n_states, n_cycles = n_cycles, verbose = verbose)
    check_sum_of_transition_array(m_P_trtB, n_states = n_states, n_cycles = n_cycles, verbose = verbose)
    
    # Iterative solution of age-dependent cSTM
    for(t in 1:n_cycles){
      ## Fill in cohort trace
      # For SoC
      m_M_SoC[t + 1, ]  <- m_M_SoC[t, ]  %*% m_P_SoC
      # For strategy A
      m_M_trtA[t + 1, ] <- m_M_trtA[t, ] %*% m_P_trtA
      # For strategy B
      m_M_trtB[t + 1, ] <- m_M_trtB[t, ] %*% m_P_trtB
    }
    
    ## Store the cohort traces in a list 
    l_m_M <- list(SoC =  m_M_SoC,
                  A   =  m_M_trtA,
                  B   =  m_M_trtB)
    names(l_m_M) <- v_names_str
    
  
    ### State rewards
    ## Scale by the cycle length 
    # Vector of state utilities under strategy SoC
    v_u_SoC    <- c(H  = u_H, 
                    S  = u_S,
                    D  = u_D) * cycle_length
    # Vector of state costs under strategy SoC
    v_c_SoC    <- c(H  = c_H_yr, 
                    S  = c_S_yr,
                    D  = c_D_yr) * cycle_length
    # Vector of state utilities under treatment A
    v_u_trtA   <- c(H  = u_H, 
                    S  = u_S, 
                    D  = u_D) * cycle_length
    # Vector of state costs under treatment A
    v_c_trtA   <- c(H  = c_H_yr + c_trtA_yr, 
                    S  = c_S_yr, 
                    D  = c_D_yr) * cycle_length
    # Vector of state utilities under treatment B
    v_u_trtB   <- c(H  = u_H, 
                    S  = u_S, 
                    D  = u_D) * cycle_length
    # Vector of state costs under treatment B
    v_c_trtB   <- c(H  = c_H_yr + c_trtB_yr, 
                    S  = c_S_yr, 
                    D  = c_D_yr) * cycle_length
    
    ## Store state rewards 
    # Store the vectors of state utilities for each strategy in a list 
    l_u   <- list(SQ = v_u_SoC,
                  A  = v_u_trtA,
                  B  = v_u_trtB)
    # Store the vectors of state cost for each strategy in a list 
    l_c   <- list(SQ = v_c_SoC,
                  A  = v_c_trtA,
                  B  = v_c_trtB)
    
    # assign strategy names to matching items in the lists
    names(l_u) <- names(l_c) <- v_names_str
    
   
    ## Loop through each strategy and calculate total utilities and costs 
    v_tot_cost <- v_tot_qaly <- vector(mode = "numeric", length = n_str)
    names(v_tot_qaly) <- names(v_tot_cost) <- v_names_str
    
    
    ### Expected discounted cost for each strategy
    v_tot_cost["Standard of Care"] <- t(m_M_SoC  %*% v_c_SoC)  %*% (v_dwe * v_wcc)
    v_tot_cost["Treatment A"]      <- t(m_M_trtA %*% v_c_trtA) %*% (v_dwe * v_wcc)
    v_tot_cost["Treatment B"]      <- t(m_M_trtB %*% v_c_trtB) %*% (v_dwe * v_wcc)
    
    ### Expected discounted QALYs for each strategy
    v_tot_qaly["Standard of Care"] <- t(m_M_SoC  %*% v_q_SoC)  %*% (v_dwe * v_wcc)
    v_tot_qaly["Treatment A"]      <- t(m_M_trtA %*% v_q_trtA) %*% (v_dwe * v_wcc)
    v_tot_qaly["Treatment B"]      <- t(m_M_trtB %*% v_q_trtB) %*% (v_dwe * v_wcc)
    
    ## Vector with discounted net monetary benefits (NMB)
    v_nmb <- v_tot_qaly * n_wtp - v_tot_cost
    
    ## data.frame with discounted costs, effectiveness and NMB
    df_ce <- data.frame(Strategy = v_names_str,
                        Cost     = v_tot_cost,
                        Effect   = v_tot_qaly,
                        NMB      = v_nmb)
    
    return(df_ce)
  }
  )
}


#------------------------------------------------------------------------------#
####             Generate a PSA input parameter dataset                     ####
#------------------------------------------------------------------------------#
#' Generate parameter sets for the probabilistic sensitivity analysis (PSA)
#'
#' \code{generate_psa_params} generates a PSA dataset of the parameters of the 
#' cost-effectiveness analysis.
#' @param n_sim Number of parameter sets for the PSA dataset
#' @param seed Seed for the random number generation
#' @return A data.frame with a PSA dataset of he parameters of the 
#' cost-effectiveness analysis
#' @export
generate_psa_params <- function(n_sim = 1000, seed = 071818){
  set.seed(seed) # set a seed to be able to reproduce the same results
  df_psa <- data.frame(
    # Transition probabilities (per cycle)
    # probability of dying
    p_HD_yr      = rbeta(n_sim, shape1 = 1, shape2 = 99),  # from healthy          
    p_SD_yr      = rbeta(n_sim, shape1 = 22.4, shape2 = 201.6),  # from sick          
    # probability of becoming sick when healthy, conditional on surviving
    p_HS_yr_SoC  = rbeta(n_sim, shape1 = 24, shape2 = 450),      # standard of care
    p_HS_yr_trtA = rbeta(n_sim, shape1 = 15, shape2 = 368),      # treatment A
    p_HS_yr_trtB = rbeta(n_sim, shape1 = 16, shape2 = 767),      # treatment B    

    ## State rewards
    # Costs
    c_H_yr       = rgamma(n_sim, shape = 16, scale = 25),        # cost of one cycle in healthy state
    c_S_yr       = rgamma(n_sim, shape = 100, scale = 10),       # cost of one cycle in sick state
    c_D       = 0,                                            # cost of one cycle in dead state
    c_trtA_yr    = 800,                                          # cost of treatment A (per cycle) in healthy state
    c_trtB_yr    = 1500,                                         # cost of treatment B (per cycle) in healthy state
    
    # Utilities
    u_H       = rbeta(n_sim, shape1 =  1.5, shape2 = 0.0015), # utility when healthy 
    u_S       = rbeta(n_sim, shape1 = 49.5, shape2 = 49.5),   # utility when sick
    u_D       = 0                                             # utility when dead
  )
  return(df_psa)
}
