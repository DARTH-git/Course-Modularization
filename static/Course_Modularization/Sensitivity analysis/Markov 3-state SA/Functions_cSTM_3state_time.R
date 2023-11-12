
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
calculate_ce_out <- function(l_params_all, n_wtp = 10000, verbose = FALSE){ # User defined
  with(as.list(l_params_all), {
    
    v_names_states  <- c("H", "S", "D")       # state names, Healthy (H), Sick (S), Dead(D)
    n_states        <- length(v_names_states) # number of health states 
    
    ## Cycle names
    v_names_cycles  <- paste("cycle", 0:n_cycles)
    
    ## Rate of dying when healthy (age-dependent) - this is now a sequence of numbers
    v_r_HD    <- r_base * rr_annual ^ (cycle_length *  (0:(n_cycles - 1)))
    
    ###  Transition Probabilities
    
    ### Converting rates to probabilities
    # p = 1 - exp( -r * cycle_length)
    p_HS_SoC  <- rate_to_prob(r = r_HS_SoC, t = cycle_length)  # probability  of becoming sick when healthy, under SoC
    p_HS_trtA <- rate_to_prob(r = r_HS_trtA, t = cycle_length) # probability of becoming sick when healthy, under treatment A
    p_HS_trtB <- rate_to_prob(r = r_HS_trtB, t = cycle_length) # probability of becoming sick when healthy, under treatment B
    p_SD      <- rate_to_prob(r = r_SD, t = cycle_length)      # probability of dying when sick
    v_p_HD    <- rate_to_prob(r = v_r_HD, t = cycle_length)    # probability of dying when healthy (vector)
    
    
    # All starting healthy
    v_m_init <- c("H" = 1, "S" = 0, "D" = 0)  
    
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
    a_P_SoC <- array(0,  # Create 3-D array
                     dim = c(n_states, n_states, n_cycles),
                     dimnames = list(v_names_states, v_names_states, 
                                     v_names_cycles[-length(v_names_cycles)])) # name the dimensions of the array 
    
    ### Fill in array
    ## Standard of Care
    # from Healthy
    a_P_SoC["H", "H", ]    <- (1 - v_p_HD) * (1 - p_HS_SoC)
    a_P_SoC["H", "S",    ] <- (1 - v_p_HD) *      p_HS_SoC
    a_P_SoC["H", "D",    ] <-      v_p_HD
    
    # from Sick
    a_P_SoC["S", "S", ] <- 1 - p_SD
    a_P_SoC["S", "D", ] <-     p_SD
    
    # from Dead
    a_P_SoC["D", "D", ] <- 1
    
    ## Treatment A
    a_P_trtA <- a_P_SoC
    a_P_trtA["H", "H", ]    <- (1 - v_p_HD) * (1 - p_HS_trtA)
    a_P_trtA["H", "S",    ] <- (1 - v_p_HD) *      p_HS_trtA
    
    ## Treatment B
    a_P_trtB <- a_P_SoC
    a_P_trtB["H", "H", ]    <- (1 - v_p_HD) * (1 - p_HS_trtB)
    a_P_trtB["H", "S",    ] <- (1 - v_p_HD) *      p_HS_trtB
    
    ## Check if transition array and probabilities are valid
    # Check that transition probabilities are in [0, 1]
    check_transition_probability(a_P_SoC,  verbose = verbose)
    check_transition_probability(a_P_trtA, verbose = verbose)
    check_transition_probability(a_P_trtB, verbose = verbose)
    # Check that all rows sum to 1
    check_sum_of_transition_array(a_P_SoC,  n_states = n_states, n_cycles = n_cycles, verbose = verbose)
    check_sum_of_transition_array(a_P_trtA, n_states = n_states, n_cycles = n_cycles, verbose = verbose)
    check_sum_of_transition_array(a_P_trtB, n_states = n_states, n_cycles = n_cycles, verbose = verbose)
    
    # Iterative solution of age-dependent cSTM
    for(t in 1:n_cycles){
      ## Fill in cohort trace
      # For SoC
      m_M_SoC[t + 1, ]  <- m_M_SoC[t, ]  %*% a_P_SoC[, , t]
      # For strategy A
      m_M_trtA[t + 1, ] <- m_M_trtA[t, ] %*% a_P_trtA[, , t]
      # For strategy B
      m_M_trtB[t + 1, ] <- m_M_trtB[t, ] %*% a_P_trtB[, , t]
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
    v_c_SoC    <- c(H  = c_H, 
                    S  = c_S,
                    D  = c_D) * cycle_length
    # Vector of state utilities under treatment A
    v_u_trtA   <- c(H  = u_H, 
                    S  = u_S, 
                    D  = u_D) * cycle_length
    # Vector of state costs under treatment A
    v_c_trtA   <- c(H  = c_H + c_trtA, 
                    S  = c_S, 
                    D  = c_D) * cycle_length
    # Vector of state utilities under treatment B
    v_u_trtB   <- c(H  = u_H, 
                    S  = u_S, 
                    D  = u_D) * cycle_length
    # Vector of state costs under treatment B
    v_c_trtB   <- c(H  = c_H + c_trtB, 
                    S  = c_S, 
                    D  = c_D) * cycle_length
    
    ## Store state rewards 
    # Store the vectors of state utilities for each strategy in a list 
    l_u   <- list(SoQ = v_u_SoC,
                  A   = v_u_trtA,
                  B   = v_u_trtB)
    # Store the vectors of state cost for each strategy in a list 
    l_c   <- list(SoQ = v_c_SoC,
                  A   = v_c_trtA,
                  B   = v_c_trtB)
    
    # assign strategy names to matching items in the lists
    names(l_u) <- names(l_c) <- v_names_str
    
    # Create empty vectors to store total utilities and costs 
    v_tot_qaly <- v_tot_cost <- vector(mode = "numeric", length = n_str)
    names(v_tot_qaly) <- names(v_tot_cost) <- v_names_str
    
    ## Loop through each strategy and calculate total utilities and costs 
    for (i in 1:n_str) {
      v_u_str <- l_u[[i]]   # select the vector of state utilities for the i-th strategy
      v_c_str <- l_c[[i]]   # select the vector of state costs for the i-th strategy
      
      ### Expected QALYs and costs per cycle 
      ## Vector of QALYs and Costs
      # Apply state rewards 
      v_qaly_str <- l_m_M[[i]] %*% v_u_str # sum the utilities of all states for each cycle
      v_cost_str <- l_m_M[[i]] %*% v_c_str # sum the costs of all states for each cycle
      
      ### Discounted total expected QALYs and Costs per strategy and apply within-cycle correction if applicable
      # QALYs
      v_tot_qaly[i] <- t(v_qaly_str) %*% (v_dwe * v_wcc)
      # Costs
      v_tot_cost[i] <- t(v_cost_str) %*% (v_dwc * v_wcc)
    }
    
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
    # Transition rates
    # rate of dying
    r_SD      = rlnorm(n_sim,meanlog =   log(0.1), sdlog =  (0.1)),  # from sick
    
    r_base     = 0.003, # rate of dying when healthy at t = 0
    rr_annual  = 1.1,  # annual increase of mortality rate
    
    # probability of becoming sick when healthy, conditional on surviving
    r_HS_SoC  = rlnorm(n_sim,meanlog =   log(0.05), sdlog = 0.05) ,      # standard of care
    r_HS_trtA =  rlnorm(n_sim,meanlog =   log(0.04), sdlog = 0.04) ,
    r_HS_trtB =  rlnorm(n_sim,meanlog =   log(0.02), sdlog = 0.02) ,
    
    ## State rewards
    # Costs
    c_H       = rgamma(n_sim, shape = 16, scale = 25),        # cost of one cycle in healthy state
    c_S       = rgamma(n_sim, shape = 100, scale = 10),       # cost of one cycle in sick state
    c_D       = 0,                                            # cost of one cycle in dead state
    c_trtA    = 800,                                          # cost of treatment A (per cycle) in healthy state
    c_trtB    = 1500,                                         # cost of treatment B (per cycle) in healthy state
    
    # Utilities
    u_H       = rbeta(n_sim, shape1 =  1.5, shape2 = 0.0015), # utility when healthy 
    u_S       = rbeta(n_sim, shape1 = 49.5, shape2 = 49.5),   # utility when sick
    u_D       = 0                                             # utility when dead
  )
  return(df_psa)
}
