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
# 3 Erasmus University, Rotterdam, The Netherlands
# 4 Harvard T.H. Chan School of Public Health, Boston, USA
# 5 University of Pittsburgh Graduate School of Public Health, Pittsburgh, PA, USA
# 6 The Hospital for Sick Children, Toronto and University of Toronto, Toronto ON, Canada

# used for sensitivity analysis

calculate_ce_out <- function (l_params_all, n_wtp = 100000) {
  with(as.list(l_params_all), {
    
    ## 03.2 Calculate internal model parameters
    
    ### Discount weight for costs and effects 
    v_dwc   <- 1 / ((1 + (d_e * cycle_length)) ^ (0:n_cycles))
    v_dwe   <- 1 / ((1 + (d_c * cycle_length)) ^ (0:n_cycles))
    
    # Within-cycle correction (WCC) - method  options Simpson's 1/3 rule, "half-cycle" or "none" 
    v_wcc    <- darthtools::gen_wcc(n_cycles = n_cycles, 
                                    method = "Simpson1/3") # vector of wcc
    
    ## 04.1 Static characteristics
    # randomly sample the sex of an individual (50% female)
    v_sex <- sample(x = c("Female", "Male"), prob = c(0.5, 0.5), size = n_i, replace = TRUE) 
    
    ## 04.2 Dynamic characteristics 
    # Specify the initial health state of the individuals 
    # everyone begins in the healthy state (in this example)
    v_M_init  <- rep("Healthy", times = n_i)   
    v_Ts_init <- rep(0, n_i)  # a vector with the time of being sick at the start of the model 
    
    
    ## 04.3 Create a dataframe with the individual characteristics 
    # create a data frame with each individual's 
    # ID number, treatment effect modifier, age and initial time in sick state 
    df_X  <- data.frame(ID = 1:n_i, Sex = v_sex, n_cycles_s = v_Ts_init, M_init = v_M_init)
    # NOTE: we use n_cycles_s for the number of times being sick, we start the data frame with the initial "history" of time being sick by saving v_Ts_init. However, during the simulation this value is updated and therefore called number of times being sick.
    head(df_X) # print the first rows of the dataframe
    
    # 05 Define Simulation Functions
    
    ## 05.1 Probability function
    #The `Probs` function updates the transition probabilities of every cycle is shown below 
    Probs <- function(M_t, df_X, Trt = "SoC") { 
      # Arguments:
      # M_t:  health state occupied at cycle t (character variable)
      # df_X: data frame with individual characteristics data 
      # Trt:  treatment
      # Returns: 
      # transition probabilities for that cycle
      
      # Treatment specific transition probabilities
      if (Trt == "SoC") {
        p_HS <- p_HS_SoC
      } else if (Trt == "A") {
        p_HS <- p_HS_trtA 
      } else if (Trt == "B") {
        p_HS <- p_HS_trtB
      }
      
      # create matrix of state transition probabilities
      m_p_t           <- matrix(0, nrow = n_states, ncol = n_i)  
      # give the state names to the rows
      rownames(m_p_t) <-  v_names_states                               
      
      # lookup baseline probability and rate of dying based on individual characteristics
      p_HD_all <- inner_join(df_X, df_p_HD, by = c("Sex"))
      p_HD     <- p_HD_all[M_t == "Healthy", "p_HD"]
      
      # update m_p_t with the appropriate probabilities 
      # (all non-death probabilities are conditional on survival)
      # transition probabilities when Healthy 
      m_p_t["Healthy", M_t == "Healthy"] <- (1 - p_HD) * (1 - p_HS)
      m_p_t["Sick",    M_t == "Healthy"] <- (1 - p_HD) *      p_HS 
      m_p_t["Dead",    M_t == "Healthy"] <-      p_HD              
      
      # transition probabilities when Sick 
      m_p_t["Healthy", M_t == "Sick"]    <-  0
      m_p_t["Sick",    M_t == "Sick"]    <-  1 - p_SD[df_X$n_cycles_s]
      m_p_t["Dead",    M_t == "Sick"]    <-      p_SD[df_X$n_cycles_s]     
      
      # transition probabilities when Dead
      m_p_t["Healthy", M_t == "Dead"]    <- 0
      m_p_t["Sick",    M_t == "Dead"]    <- 0
      m_p_t["Dead",    M_t == "Dead"]    <- 1  
      
      return(t(m_p_t))
    }      
    
    ## 05.2 Cost function
    # The `Costs` function estimates the costs at every cycle.
    Costs <- function (M_t, Trt = "SoC") {
      # Arguments:
      # M_t: health state occupied at cycle t (character variable)
      # Returns: 
      # costs accrued in this cycle
      # Trt:  treatment
      
      # Treatment specific transition costs
      if (Trt == "SoC") {
        c_trt <- 0
      } else if (Trt == "A") {
        c_trt <- c_trtA
      } else if (Trt == "B") {
        c_trt <- c_trtB
      }
      
      c_t <- c()
      c_t[M_t == "Healthy"] <- c_H + c_trt  # costs accrued by being healthy this cycle
      c_t[M_t == "Sick"]    <- c_S          # costs accrued by being sick this cycle
      c_t[M_t == "Dead"]    <- c_D          # costs at dead state
      
      return(c_t)  # return costs accrued this cycle
    }
    
    ## 05.3 Health outcome function
    #The `Effs` function to update the utilities at every cycle.
    Effs <- function (M_t, cl = 1) {
      # Arguments:
      # M_t: health state occupied at cycle t (character variable)
      # cl:   cycle length (default is 1)
      # Returns: 
      # QALYs accrued this cycle 
      
      q_t <- c() 
      q_t[M_t == "Healthy"] <- u_H  # utility for being healthy this cycle
      q_t[M_t == "Sick"]    <- u_S  # utility for being sick this cycle
      q_t[M_t == "Dead"]    <- u_D  # utility for dead state
      
      QALYs <- q_t * cl  # calculate the QALYs during cycle t
      return(QALYs)      # return the QALYs accrued this cycle
    }
    
    # Microsimulation
    MicroSim <- function(n_i, df_X, seed = 1, Trt = "SoC") {
      # Arguments:  
      # n_i: number of individuals
      # df_X: data frame with individual data 
      # seed: seed for the random number generator, default is 1
      # Trt: treatment 
      # Returns:
      # results: data frame with total cost and QALYs
      
      set.seed(seed) # set a seed to be able to reproduce the same results
      
      # create three matrices called m_M, m_C and m_E
      # number of rows is equal to the n_i, the number of columns is equal to n_cycles 
      # (the initial state and all the n_cycles cycles)
      # m_M is used to store the health state information over time for every individual
      # m_C is used to store the costs information over time for every individual
      # m_E is used to store the effects information over time for every individual
      
      m_M <- m_C <- m_E <-  matrix(nrow = n_i, ncol = n_cycles + 1, 
                                   dimnames = list(paste("ind"  , 1:n_i, sep = " "), 
                                                   paste("cycle", 0:n_cycles, sep = " ")))  
      
      m_M[, 1] <- as.character(df_X$M_init) # initial health state
      m_C[, 1] <- Costs(m_M[, 1])  # costs accrued during cycle 0
      m_E[, 1] <- Effs(m_M[, 1], cl = 1)   # QALYs accrued during cycle 0
      
      # open a loop for time running cycles 1 to n_cycles 
      for (t in 1:n_cycles) {
        # calculate the transition probabilities for the cycle based on health state t
        m_P <- Probs(m_M[, t], df_X, Trt = Trt)
        # check if transition probabilities are between 0 and 1
        check_transition_probability(m_P, verbose = TRUE)
        # check if each of the rows of the transition probabilities matrix sum to one
        check_sum_of_transition_array(m_P, n_rows = n_i, n_cycles = n_cycles, verbose = TRUE)
        
        # sample the next health state and store that state in matrix m_M
        m_M[, t + 1]  <- samplev(m_P, 1)    
        # calculate costs per individual during cycle t + 1
        m_C[, t + 1]  <- Costs(m_M[, t + 1], Trt = Trt)  
        # calculate QALYs per individual during cycle t + 1
        m_E[, t + 1]  <- Effs (m_M[, t + 1], cl = 1)  
        
        # update time since illness onset for t + 1 
        # NOTE: this code has a "reset of history" for time being sick
        # once someone is not "Sick" anymore, we reset n_cycles_s (set back to zero)
        # when you don't want a "reset" replace the last zero with teh current value
        df_X$n_cycles_s <- if_else(m_M[, t + 1] == "Sick", df_X$n_cycles_s + 1, 0) 
        
        # Display simulation progress
        if(t/(n_cycles/10) == round(t/(n_cycles/10), 0)) { # display progress every 10%
          cat('\r', paste(t/n_cycles * 100, "% done", sep = " "))
        }
        
      } # close the loop for the time points 
      
      # calculate  
      tc      <- m_C %*% (v_dwc * v_wcc)   # total (discounted and cycle corrected) cost per individual
      te      <- m_E %*% (v_dwe * v_wcc)   # total (discounted and cycle corrected) QALYs per individual 
      tc_hat  <- mean(tc)       # average (discounted) cost 
      te_hat  <- mean(te)       # average (discounted) QALY  
      # store the results from the simulation in a list
      results <- list(m_M = m_M, m_C = m_C, m_E = m_E, tc = tc , te = te, 
                      tc_hat = tc_hat, te_hat = te_hat)   
      
      return(results)  # return the results
      
    } # end of the `MicroSim` function  
    
    # Run the simulation 
    outcomes_SoC  <- MicroSim(n_i = n_i, df_X = df_X, seed = 1, Trt = "SoC")
    outcomes_trtA <- MicroSim(n_i = n_i, df_X = df_X, seed = 1, Trt = "A")
    outcomes_trtB <- MicroSim(n_i = n_i, df_X = df_X, seed = 1, Trt = "B")
    
    # Cost-Effectiveness Analysis
    # store the mean costs of each strategy in a new variable C (vector of costs)
    v_C <- c(outcomes_SoC$tc_hat, outcomes_trtA$tc_hat, outcomes_trtB$tc_hat)
    # store the mean QALYs of each strategy in a new variable E (vector of effects)
    v_E <- c(outcomes_SoC$te_hat, outcomes_trtA$te_hat, outcomes_trtB$te_hat)
    
    ## Vector with discounted net monetary benefits (NMB)
    v_nmb_d <- v_E * n_wtp - v_C
    
    # Dataframe with discounted costs, effectiveness and NMB
    df_ce <- data.frame(Strategy = v_names_str,
                        Cost     = v_C,
                        Effect   = v_E,
                        NMB      = v_nmb_d)
    return(df_ce)
  }
  )
}
    