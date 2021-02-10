
# Developed by the Decision Analysis in R for Technologies in Health (DARTH) group
# Fernando Alarid-Escudero, PhD (1) 
# Eva A. Enns, MS, PhD (2)	
# M.G. Myriam Hunink, MD, PhD (3,4)
# Hawre J. Jalal, MD, PhD (5) 
# Eline M. Krijkamp, MSc (3)
# Petros Pechlivanoglou, PhD (6) 

# In collaboration of: 		
# 1 Center for Research and Teaching in Economics (CIDE), Drug Policy Program, Mexico
# 2 University of Minnesota School of Public Health, Minneapolis, MN, USA
# 3 Erasmus MC, Rotterdam, The Netherlands
# 4 Harvard T.H. Chan School of Public Health, Boston, USA
# 5 University of Pittsburgh Graduate School of Public Health, Pittsburgh, PA, USA
# 6 The Hospital for Sick Children, Toronto and University of Toronto, Toronto ON, Canada

# used for sensitivity analysis

calculate_ce_out <- function (l_params_all, n_wtp = 100000) {
  with(as.list(l_params_all), {
    
    ## 04.1 Static characteristics
    v_x <- runif(n_i, min = 0.95, max = 1.05)  # treatment effect modifier at baseline               
    # Dynamic characteristics
    # sample from age distribution an initial age for every individual
    v_age0 <- sample(x = dist_Age$age, prob = dist_Age$prop, size = n_i, replace = TRUE) 
    
    ## 04.2 Dynamic characteristics 
    # sample from age distribution an initial age for every individual
    v_age0 <- sample(x = dist_Age$age, prob = dist_Age$prop, size = n_i, replace = TRUE) 
    # a vector with the time of being sick at the start of the model  
    
    # Specify the initial health state of the individuals 
    # everyone begins in the healthy state (in this example)
    # a vector with the initial health state for all individuals
    v_M_init  <- rep("H", n_i)
    v_Ts_init <- rep(0, n_i)  # since all individuals start healthy this value is zero for everyone
    
    ## 04.3 Create a dataframe with the individual characteristics 
    df_X <- data.frame(ID = 1:n_i, x = v_x, Age = v_age0, n_ts = v_Ts_init) # create a dataframe with an ID number for every individual, the individual treatment effect modifier and the age of the individuals 
    
    ## 05.1 Probability function
    Probs <- function(M_t, df_X, t) { 
      # Arguments:
      # M_t:  health state occupied by individual i at cycle t (character variable)
      # df_X: data frame with individual characteristics data 
      # t:    current cycle 
      # Returns: 
      # transition probabilities for that cycle
      
      # create matrix of state transition probabilities  
      m_p_t           <- matrix(0, nrow = n_states, ncol = n_i) 
      rownames(m_p_t) <-  v_names_states  # give the state names to the rows
      
      # look up baseline probability and rate of dying based on individual characteristics
      p_HD_all <- inner_join(df_X, p_mort, by = c("Age"))
      p_HD     <- p_HD_all[M_t == "H","p_HD"]
      
      # update the m_p with the appropriate probabilities   
      # transition probabilities when healthy
      m_p_t[, M_t == "H"]  <- rbind((1 - p_HD) * (1 - p_HS1),
                                    (1 - p_HD) *      p_HS1 , 
                                    0 ,
                                    p_HD               )                              
      # transition probabilities when sick
      m_p_t[, M_t == "S1"] <- rbind((1 - p_S1D[df_X$n_ts]) * p_S1H,
                                    (1 - p_S1D[df_X$n_ts]) * (1 - p_S1H - p_S1S2),
                                    (1 - p_S1D[df_X$n_ts]) *              p_S1S2 , 
                                         p_S1D[df_X$n_ts]                        )  
      # transition probabilities when sicker
      m_p_t[, M_t == "S2"] <- rbind(0, 
                                    0, 
                                    1 - p_S2D, 
                                        p_S2D)                                            
      # transition probabilities when dead   
      m_p_t[, M_t == "D"]  <- rbind(0, 0, 0, 1)                                                        
      
      return(t(m_p_t))
    }           
    
    ## 05.2 Cost function
    Costs <- function (M_t, Trt = FALSE) {
      # Arguments:
      # M_t: health state occupied by individual i at cycle t (character variable)
      # Trt: is the individual being treated? (default is FALSE) 
      # Returns:
      # costs accrued in this cycle
      
      c_t <- 0                                # by default the cost for everyone is zero 
      c_t[M_t == "H"]  <- c_H                 # update the cost if healthy
      c_t[M_t == "S1"] <- c_S1 + c_Trt * Trt  # update the cost if sick conditional on treatment
      c_t[M_t == "S2"] <- c_S2 + c_Trt * Trt  # update the cost if sicker conditional on treatment
      c_t[M_t == "D"]  <- c_D                 # update the cost if dead
      
      return(c_t)  # return the costs
    }
    
    ## 05.3 Health outcome function
    Effs <- function (M_t, df_X, Trt = FALSE, cl = 1) {
      # Arguments:  
      # M_t:  health state occupied by individual i at cycle t (character variable)
      # df_X: data frame with individual characteristics data 
      # Trt:  is the individual treated? (default is FALSE) 
      # cl:   cycle length (default is 1)
      # Returns:
      # QALYs accrued this cycle
      
      u_t <- 0                                 # by default the utility for everyone is zero
      u_t[M_t == "H"]  <- u_H                  # update the utility if healthy
      u_t[M_t == "S1" & Trt == FALSE] <- u_S1  # update the utility if sick
      # update the utility if sick but on treatment (adjust for individual effect modifier) 
      u_t[M_t == "S1" & Trt == TRUE]  <- u_Trt * df_X$x[M_t == "S1"]  
      u_t[M_t == "S2"] <- u_S2                 # update the utility if sicker
      u_t[M_t == "D"]  <- u_D                  # update the utility if dead
      
      QALYs <-  u_t * cl  # calculate the QALYs during cycle t
      return(QALYs)       # return the QALYs
    }
    
    # 06 Run Microsimulation
    MicroSim <- function(n_i, df_X, Trt = FALSE, seed = 1) {
      # Arguments:  
      # n_i:  number of individuals
      # df_X: data frame with individual characteristics data 
      # Trt:  is this the individual receiving treatment? (default is FALSE)
      # seed: default is 1
      # Returns:
      # results: data frame with total cost and QALYs  
      
      set.seed(seed) # set the seed
      
      n_states <- length(v_names_states) # the number of health states
      
      # create three matrices called m_M, m_C and m_E
      # number of rows is equal to the n_i, the number of columns is equal to n_t  
      # (the initial state and all the n_t cycles)
      # m_M is used to store the health state information over time for every individual
      # m_C is used to store the costs information over time for every individual
      # m_E is used to store the effects information over time for every individual
      
      m_M <- m_C <- m_E <- m_Ts <-  matrix(nrow = n_i, ncol = n_t + 1, 
                                           dimnames = list(paste("ind"  , 1:n_i, sep = " "), 
                                                           paste("cycle", 0:n_t, sep = " ")))  
      
      m_M [, 1] <- v_M_init    # initial health state at cycle 0 for individual i
      
      # calculate costs per individual during cycle 0
      m_C[, 1]  <- Costs(m_M[, 1], Trt)     
      # calculate QALYs per individual during cycle 0
      m_E[, 1]  <- Effs (m_M[, 1], df_X, Trt)   
      
      # open a loop for time running cycles 1 to n_t 
      for (t in 1:n_t) {
        # calculate the transition probabilities for the cycle based on  health state t
        m_P <- Probs(m_M[, t], df_X, t)             
        # check if transition probabilities are between 0 and 1
        check_transition_probability(m_P, verbose = TRUE)
        # check if checks if each of the rows of the transition probabilities matrix sum to one
        check_sum_of_transition_array(m_P, n_rows = n_i, n_cycles = n_t, verbose = TRUE)
        # sample the current health state and store that state in matrix m_M 
        m_M[, t + 1]  <- samplev(m_P)                  
        # calculate costs per individual during cycle t + 1
        m_C[, t + 1]  <- Costs(m_M[, t + 1], Trt)         
        # calculate QALYs per individual during cycle t + 1
        m_E[, t + 1]  <- Effs(m_M[, t + 1], df_X, Trt)    
        
        # update time since illness onset for t + 1 
        df_X$n_ts <- if_else(m_M[, t + 1] == "S1", df_X$n_ts + 1, 0) 
        # update the age of individuals that are alive
        df_X$Age[m_M[, t + 1] != "D"]  <- df_X$Age[m_M[, t + 1] != "D"] + 1
        
        # Display simulation progress
        if(t/(n_t/10) == round(t/(n_t/10), 0)) { # display progress every 10%
          cat('\r', paste(t/n_t * 100, "% done", sep = " "))
        }
        
      } # close the loop for the time points 
      
      # calculate  
      tc      <- m_C %*% v_dwc  # total (discounted) cost per individual
      te      <- m_E %*% v_dwe  # total (discounted) QALYs per individual 
      tc_hat  <- mean(tc)       # average (discounted) cost 
      te_hat  <- mean(te)       # average (discounted) QALYs
      
      # store the results from the simulation in a list
      results <- list(m_M = m_M, m_C = m_C, m_E = m_E, tc = tc , te = te, 
                      tc_hat = tc_hat, te_hat = te_hat)   
      
      return(results)  # return the results
      
    } # end of the MicroSim function  
    
    # By specifying all the arguments in the `MicroSim()` the simulation can be started
    # In this example the outcomes are of the simulation are stored in the variables `outcomes_no_tr` and `outcomes_trt`.
    
    # Run the simulation for both no treatment and treatment options
    outcomes_no_trt <- MicroSim(n_i, df_X, Trt = FALSE, seed = 1)
    outcomes_trt    <- MicroSim(n_i, df_X, Trt = TRUE,  seed = 1)
    
    # 08 Cost-Effectiveness Analysis
    # store the mean costs of each strategy in a new variable C (vector of costs)
    v_C <- c(outcomes_no_trt$tc_hat, outcomes_trt$tc_hat)
    # store the mean QALYs of each strategy in a new variable E (vector of effects)
    v_E <- c(outcomes_no_trt$te_hat, outcomes_trt$te_hat)
    
    ## Vector with discounted net monetary benefits (NMB)
    v_nmb_d <- v_E * n_wtp - v_C
    
    ## Dataframe with discounted costs, effectiveness and NMB
    df_ce <- data.frame(Strategy = v_names_str,
                        Cost     = v_C,
                        Effect   = v_E,
                        NMB      = v_nmb_d)
    return(df_ce)
  }
  )
}