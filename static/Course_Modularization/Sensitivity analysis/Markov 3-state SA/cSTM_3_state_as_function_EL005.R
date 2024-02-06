## 10.1 Model input for SA
library("darthtools")
library("dampack")




l_params_all <- list(
  # Transition probabilities
  # probability of dying
  r_HS_SoC   = 0.05,  # rate of becoming sick when healthy, under standard of care
  r_HS_trtA  = 0.04,  # rate of becoming sick when healthy, under treatment A
  r_HS_trtB  = 0.02,  # rate of becoming sick when healthy, under treatment B
  r_SD       = 0.1,   # rate of dying when sick
  r_base     = 0.003, # rate of dying when healthy at t = 0
  rr_annual  = 1.1,  # annual increase of mortality rate
  # Costs 
  c_H       = 400,   # cost of one cycle in healthy state
  c_S       = 1000,  # cost of one cycle in sick state
  c_D       = 0,     # cost of one cycle in dead state
  c_trtA    = 800,   # cost of treatment A (per cycle) in healthy state
  c_trtB    = 1500,  # cost of treatment B (per cycle) in healthy state
  # Utilities 
  u_H       = 1,     # utility when healthy 
  u_S       = 0.5,   # utility when sick
  u_D       = 0,     # utility when dead
  # Discount rates
  d_e       = 0.03,  # discount rate per cycle equal discount of costs and QALYs by 3%
  d_c       = 0.03,  # discount rate per cycle equal discount of costs and QALYs by 3%
  # Time horizon
  n_cycles  = 60,     # simulation time horizon (number of cycles)
  cycle_length = 1
)




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
calculate_ce_out <- function(l_params_all, n_wtp = 10000 , verbose = verbose){ # User defined
  with(as.list(l_params_all), {
    ########################### Process model inputs ###########################
    ## Model states
    v_names_cycles  <- paste("cycle", 0:n_cycles)         # cycle names
    v_names_states  <- c("H", "S", "D")       # state names
    n_states        <- length(v_names_states)             # number of health states 
    
    ### Strategies 
    v_names_str     <- c("Standard of Care",         # store the strategy names
                         "Treatment A", 
                         "Treatment B")  
    n_str           <- length(v_names_str)           # number of strategies
   
     ## Within-cycle correction (WCC) using Simpson's 1/3 rule, alternative method = "half-cycle"
    v_wcc <- gen_wcc(n_cycles = n_cycles,  method = "Simpson1/3")
    
    ### Cycle-specific discount weight for costs and effects 
    v_dwc   <- 1 / ((1 + (d_e * cycle_length)) ^ (0:n_cycles))
    v_dwe   <- 1 / ((1 + (d_c * cycle_length)) ^ (0:n_cycles))
    
    # rate of dying when healthy (time-dependent) - this is now a sequence of numbers
    v_r_HD    <- r_base * rr_annual^(cycle_length * (0:(n_cycles - 1)))
    
    ###  Transition Probabilities
  
    ### Converting rates to probabilities
    # function in darthtoolds, rate_to_prob based on p = 1 - exp( -r * cycle_length)
    p_HS_SoC  <- rate_to_prob(r = r_HS_SoC,  t = cycle_length) # probability  of becoming sick when healthy, under SoC
    p_HS_trtA <- rate_to_prob(r = r_HS_trtA, t = cycle_length) # probability of becoming sick when healthy, under treatment A
    p_HS_trtB <- rate_to_prob(r = r_HS_trtB, t = cycle_length) # probability of becoming sick when healthy, under treatment B
    p_SD      <- rate_to_prob(r = r_SD,      t = cycle_length) # probability of dying when sick
    v_p_HD    <- rate_to_prob(r = v_r_HD,    t = cycle_length) # probability of dying when healthy (vector)
    

    
    
    # All starting healthy
    v_m_init <- c("H" = 1, "S" = 0, "D" = 0)  
    
    ###################### Construct state-transition models ###################
    ### Initialize cohort trace for SoC. Set all values to zero
    # dimensions are number of cycles + 1 (for initial state) and number of states
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
    a_P_SoC["H", "H", ] <- (1 - v_p_HD) * (1 - p_HS_SoC)
    a_P_SoC["H", "S", ] <- (1 - v_p_HD) *      p_HS_SoC
    a_P_SoC["H", "D", ] <-      v_p_HD
    
    # from Sick
    a_P_SoC["S", "S", ] <- 1 - p_SD
    a_P_SoC["S", "D", ] <-     p_SD
    
    # from Dead
    a_P_SoC["D", "D", ] <- 1
    
    ## Treatment A
    a_P_trtA <- a_P_SoC # Store array, only replace values that are different in the treatment
    a_P_trtA["H", "H", ] <- (1 - v_p_HD) * (1 - p_HS_trtA)
    a_P_trtA["H", "S", ] <- (1 - v_p_HD) *      p_HS_trtA
    
    ## Treatment B
    a_P_trtB <- a_P_SoC
    a_P_trtB["H", "H", ] <- (1 - v_p_HD) * (1 - p_HS_trtB)
    a_P_trtB["H", "S", ] <- (1 - v_p_HD) *      p_HS_trtB
    
    ## Check if transition array and probabilities are valid
    # Check that transition probabilities are in [0, 1]
    check_transition_probability(a_P_SoC,  verbose = TRUE)
    check_transition_probability(a_P_trtA, verbose = TRUE)
    check_transition_probability(a_P_trtB, verbose = TRUE)
    # Check that all rows sum to 1
    check_sum_of_transition_array(a_P_SoC,  n_states = n_states, n_cycles = n_cycles, verbose = TRUE)
    check_sum_of_transition_array(a_P_trtA, n_states = n_states, n_cycles = n_cycles, verbose = TRUE)
    check_sum_of_transition_array(a_P_trtB, n_states = n_states, n_cycles = n_cycles, verbose = TRUE)
    
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
    

    ## Scale by the cycle length
    
    # Standard of Care
    # vector of state QALYs accrued each cycle
    v_u_SoC    <- c(H  = u_H, 
                    S  = u_S,
                    D  = u_D) * cycle_length
    # vector of state costs accrued each cycle
    v_c_SoC    <- c(H  = c_H, 
                    S  = c_S,
                    D  = c_D) * cycle_length
    
    # Treatment A
    # vector of state QALYs accrued each cycle
    v_u_trtA   <- c(H  = u_H, 
                    S  = u_S, 
                    D  = u_D) * cycle_length
    # vector of state costs accrued each cycle
    v_c_trtA   <- c(H  = c_H + c_trtA, 
                    S  = c_S, 
                    D  = c_D) * cycle_length
    
    # Treatment B
    # vector of state QALYs accrued each cycle
    v_u_trtB   <- c(H  = u_H, 
                    S  = u_S, 
                    D  = u_D) * cycle_length
    # vector of state costs accrued each cycle
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
    
    
    ## data.frame with discounted costs, effectiveness and NMB
    df_ce <- data.frame(Strategy = v_names_str,
                        Cost     = v_tot_cost,
                        Effect   = v_tot_qaly)
    
    return(df_ce)
  }
  )
}


df_ce_out <- calculate_ce_out(l_params_all)


## Incremental cost-effectiveness ratios (ICERs) 
df_cea <- calculate_icers(cost       = df_ce_out$Cost, 
                          effect     = df_ce_out$Effect,
                          strategies = df_ce_out$Strategy)
df_cea

