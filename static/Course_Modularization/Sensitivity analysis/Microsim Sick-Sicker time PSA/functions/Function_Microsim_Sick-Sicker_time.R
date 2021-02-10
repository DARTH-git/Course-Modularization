
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
    
    # Static characteristics
    v_x      <- runif(n_i, min = 0.95, max = 1.05) # treatment effect modifier at baseline                
    # Dynami characteristics
    v_Ts_init <- rep(0, n_i)         # since all individuals start healthy this value is zero for everyone
    v_age0   <- sample(x = dist_Age$age, prob = dist_Age$prop, size = n_i, replace = TRUE) # sample from age distribution an initial age for every individual
    ## Create the data.frame
    df_X    <- data.frame(ID = 1:n_i, x = v_x, Age = v_age0, n_ts = v_Ts_init) # create a dataframe with an ID number for every individual, the individual treatment effect modifier and the age of the individuals 
    
    v_M_init  <- rep("H", n_i)
    
    # The Probs function that updates the transition probabilities of every cycle is shown below.
    Probs <- function(M_t, df_X, t) { 
      # Arguments:
        # M_t:  health state occupied at cycle t (character vector)
        # df_X: dataframe with individual characteristics
        # t:    current cycle 
      # Returns: 
        # transition probabilities for that cycle
      
      m_p_t           <- matrix(0, nrow = n_states, ncol = n_i)  # create matrix of state transition probabilities
      rownames(m_p_t) <-  v_names_states                         # give the state names to the rows
      
      # look up baseline probability and rate of dying based on individual characteristics
      p_HD_all <- inner_join(df_X, p_mort, by = c("Age"))
      p_HD     <- p_HD_all[M_t == "H","p_HD"]
      
      # update the v_p with the appropriate probabilities   
      # transition probabilities when healthy
      m_p_t[, M_t == "H"]  <- rbind((1 - p_HD) * (1 - p_HS1), 
                                    (1 - p_HD) *      p_HS1 , 
                                                          0 ,
                                         p_HD               )                             
      m_p_t[, M_t == "S1"] <- rbind((1 - p_S1D[df_X$n_ts]) *      p_S1H,
                                    (1 - p_S1D[df_X$n_ts]) * (1 - p_S1H - p_S1S2),
                                    (1 - p_S1D[df_X$n_ts]) *              p_S1S2 , 
                                         p_S1D[df_X$n_ts]                        ) 
      m_p_t[, M_t == "S2"] <- rbind(0, 
                                    0,
                                    1 - p_S2D,
                                        p_S2D)           
      m_p_t[, M_t == "D"]  <- rbind(0, 0, 0, 1)         
      return(t(m_p_t))
    }       
    
    # The Costs function estimates the costs at every cycle.
    Costs <- function (M_t, Trt = FALSE) {
      # Arguments:
        # M_t: health state occupied by individual i at cycle t (character variable)
        # Trt: is the individual being treated? (default is FALSE) 
      # Returns:
        # costs accrued in this cycle
      
      c_t <- 0                                 # by default the cost for everyone is zero 
      c_t[M_t == "H"]  <- c_H                  # update the cost if healthy
      c_t[M_t == "S1"] <- c_S1 + c_Trt * Trt   # update the cost if sick conditional on treatment
      c_t[M_t == "S2"] <- c_S2 + c_Trt * Trt   # update the cost if sicker conditional on treatment
      c_t[M_t == "D"]  <- c_D                  # update the cost if dead
      
      return(c_t)        		                   # return the costs
    }
    
    # The Effs function to update the utilities at every cycle.
    Effs <- function (M_t, df_X, Trt = FALSE, cl = 1) {
      # Arguments:
        # M_t:  health state occupied by individual i at cycle t (character variable)
        # df_X: data frame with individual characteristics data 
        # t:    current cycle 
      # Returns: 
        # transition probabilities for that cycle
      
      u_t <- 0                                       # by default the utility for everyone is zero
      u_t[M_t == "H"]  <- u_H                        # update the utility if healthy
      u_t[M_t == "S1" & Trt == FALSE] <- u_S1        # update the utility if sick
      u_t[M_t == "S1" & Trt == TRUE]  <- u_Trt * df_X$x[M_t == "S1"]  # update the utility if sick but on treatment (adjust for individual effect modifier) 
      u_t[M_t == "S2"] <- u_S2                       # update the utility if sicker
      u_t[M_t == "D"]  <- u_D                        # update the utility if dead
      
      QALYs <-  u_t * cl            # calculate the QALYs during cycle t
      return(QALYs)                 # return the QALYs
    }
    
    #### 06 Run Microsimulation ####
    MicroSim <- function(n_i, df_X, Trt = FALSE, seed = 1) {
      # Arguments:  
        # n_i:     number of individuals
        # df_X     data frame with individual characteristics data 
        # Trt:     is this the individual receiving treatment? (default is FALSE)
        # seed:    default is 1
      # Returns:
        # results: data frame with total cost and QALYs  
      
      set.seed(seed) # set the seed
      
      # create three matrices called m_M, m_C and m_E
      # number of rows is equal to the n_i, the number of columns is equal to n_t  (the initial state and all the n_t cycles)
      # m_M is used to store the health state information over time for every individual
      # m_C is used to store the costs information over time for evey individual
      # m_E is used to store the effects information over time for every individual
      
      m_M <- m_C <- m_E <- m_Ts <-  matrix(nrow = n_i, ncol = n_t + 1, 
                                           dimnames = list(paste("ind"  , 1:n_i, sep = " "), 
                                                           paste("cycle", 0:n_t, sep = " ")))  
      
      m_M [, 1] <- v_M_init    # initial health state at cycle 0 for individual i
      
      m_C[, 1]  <- Costs(m_M[, 1], Trt)        # calculate costs per individual during cycle 0
      m_E[, 1]  <- Effs (m_M[, 1], df_X, Trt)  # calculate QALYs per individual during cycle 0
      
      # open a loop for time running cycles 1 to n_t 
      for (t in 1:n_t) {
        v_p <- Probs(m_M[, t], df_X, t)                 # calculate the transition probabilities for the cycle based on  health state t
        m_M[, t + 1]  <- samplev(v_p)                   # sample the current health state and store that state in matrix m_M 
        m_C[, t + 1]  <- Costs(m_M[, t + 1], Trt)       # calculate costs per individual during cycle t + 1
        m_E[, t + 1]  <- Effs(m_M[, t + 1], df_X, Trt)  # calculate QALYs per individual during cycle t + 1
     
        df_X$n_ts <- if_else(m_M[, t + 1] == "S1", df_X$n_ts + 1, 0) # update time since illness onset for t + 1 
        df_X$Age[m_M[,t + 1] != "D"]  <- df_X$Age[m_M[, t + 1] != "D"] + 1
        
        
      } # close the loop for the time points 
      
      # calculate  
      tc     <- m_C %*% v_dwc  # total (discounted) cost per individual
      te     <- m_E %*% v_dwe  # total (discounted) QALYs per individual 
      tc_hat <- mean(tc)       # average (discounted) cost 
      te_hat <- mean(te)       # average (discounted) QALYs
      
      # store the results from the simulation in a list
      results <- list(m_M = m_M, m_C = m_C, m_E = m_E, tc = tc , te = te, tc_hat = tc_hat, te_hat = te_hat)   
      return(results)  # return the results
    } # end of the MicroSim function  
    
    ### Run the simulation for both no treatment and treatment options
    outcomes_no_trt  <- MicroSim(n_i, df_X, Trt = FALSE, seed = 1)
    outcomes_trt     <- MicroSim(n_i, df_X, Trt = TRUE, seed = 1)
    
    ## Vector with total discounted mean Costs and QALYs
    v_tc_d  <- c(outcomes_no_trt$tc_hat, outcomes_trt$tc_hat)
    v_tu_d  <- c(outcomes_no_trt$te_hat, outcomes_trt$te_hat)
    
    ## Vector with discounted net monetary benefits (NMB)
    v_nmb_d <- v_tu_d * n_wtp - v_tc_d
    
    ## Dataframe with discounted costs, effectiveness and NMB
    df_ce <- data.frame(Strategy = v_names_str,
                        Cost     = v_tc_d,
                        Effect   = v_tu_d,
                        NMB      = v_nmb_d)
    return(df_ce)
  }
  )
}