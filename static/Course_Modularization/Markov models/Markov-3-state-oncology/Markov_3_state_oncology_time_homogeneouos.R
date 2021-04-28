################################ Initial setup ############################### 
rm(list = ls())    # remove any variables in R's memory 

### Load packages
library(darthtools)
library(dampack) 
library(scales)
library(doParallel)

### Load supplementary functions
source("funtions/Functions.R")

################################ Model input ################################# 
## General setup
cycle_length <- 1  # in months
n_cycles     <- 60 # time horizon, number of cycles
v_names_states <- c("SD", "PD", "D") # the 3 health states of the model:
# Stable Disease (SD), Progressed Disease (PD). Dead (D)
n_states    <- length(v_names_states)     # number of health states 

# Discounting factors
d_c         <- 0.035  # discount rate for costs 
d_e         <- 0.035  # discount rate for QALYs

# Strategies
v_names_str <- c("Best Supportive Care (BSC)",  # store the strategy names
                 "Comparator 1 (C1)", 
                 "Intervention (Int)") 
n_str       <- length(v_names_str)        # number of strategies

# Within-cycle correction (WCC) using Simpson's 1/3 rule
v_wcc <- darthtools::gen_wcc(n_cycles = n_cycles, 
                             method = "Simpson1/3") # vector of wcc

## Transition rates (per cycle) under BSC 
r_SD_PD_BSC <- 1 / (200/28) # constant rate of progressing
r_SD_D_BSC  <- 1 / (300/28) # constant rate of dying from SD
r_PD_D_BSC  <- 1 / (100/28) # constant rate of dying from PD

## Effectiveness of interventions as hazard ratios (HR)
## Comparator 1 vs. BSC
hr_SD_PD_C1  <- 0.80
hr_SD_D_C1   <- 0.90
hr_PD_D_C1   <- 0.95
## Intervention vs. BSC
hr_SD_PD_Int <- 0.75
hr_SD_D_Int  <- 0.85
hr_PD_D_Int  <- 0.90

### State rewards
## Costs
# BSC
c_SD_MedRes_BSC  <- 600 # Medical Resource Utilization Costs in SD (monitoring, tests, etc)
c_PD_MedRes_BSC  <- 700 # Medical Resource Utilization Costs in PD (monitoring, tests, etc)
c_SD_DrugAcq_BSC <- 500 # Drug acquisition 1L (treatment duration= time spent in SD)
c_SD_DrugAdm_BSC <- 10  # Drug administration 1L (treatment duration= time spent in SD)
c_AE_BSC         <- 1200 # AE costs (one off at the start)
# Comparator 1
c_SD_MedRes_C1  <- 600 # Medical Resource Utilization Costs in SD (monitoring, tests, etc)
c_PD_MedRes_C1  <- 700 # Medical Resource Utilization Costs in PD (monitoring, tests, etc)
c_SD_DrugAcq_C1 <- 1500 # Drug acquisition 1L (treatment duration= time spent in SD)
c_SD_DrugAdm_C1 <- 300  # Drug administration 1L (treatment duration= time spent in SD)
c_AE_C1         <- 1200 # AE costs (one off at the start)
# Intervention
c_SD_MedRes_Int  <- 600 # Medical Resource Utilization Costs in SD (monitoring, tests, etc)
c_PD_MedRes_Int  <- 700 # Medical Resource Utilization Costs in PD (monitoring, tests, etc)
c_SD_DrugAcq_Int <- 2500 # Drug acquisition 1L (treatment duration= time spent in SD)
c_SD_DrugAdm_Int <- 300  # Drug administration 1L (treatment duration= time spent in SD)
c_AE_Int         <- 1500 # AE costs (one off at the start)

# Vector with intervention-specific AE costs
v_c_AE <- c(BSC = c_AE_BSC, 
            C1  = c_AE_C1,
            Int = c_AE_Int)
# cost of being dead (per cycle)
c_D    <- 0     

## Utilities
# SD state
u_SD_BSC <- 0.890 # under BSC
u_SD_C1  <- 0.890 # under C1
u_SD_Int <- 0.898 # under Int
# PD state
u_PD_BSC <- 0.797 # under BSC
u_PD_C1  <- 0.797 # under C1
u_PD_Int <- 0.797 # under Int
# D state
u_D    <- 0     # utility when Dead 
# Disutility from AEs applied in the 1st cycle
du_BSC <- -0.100 # under BSC
du_C1  <- -0.025 # under C1
du_Int <- -0.050 # under Int

# Vector with intervention-specific AE costs
v_du_AE <- c(BSC = du_BSC, 
             C1  = du_C1,
             Int = du_Int)

# Discount weight for costs and effects
v_dwc  <- 1 / ((1 + d_e) ^ (0:n_cycles))
v_dwe  <- 1 / ((1 + d_c) ^ (0:n_cycles))

### Process model inputs
## Compute transition probabilities from rates under BSC
p_SD_PD_BSC <- rate_to_prob(r_SD_PD_BSC, t = cycle_length)
p_SD_D_BSC  <- rate_to_prob(r_SD_D_BSC , t = cycle_length)
p_PD_D_BSC  <- rate_to_prob(r_PD_D_BSC , t = cycle_length)

## Compute transition probabilities from rates under Comparator 1
p_SD_PD_C1 <- rate_to_prob(r_SD_PD_BSC * hr_SD_PD_C1, t = cycle_length)
p_SD_D_C1  <- rate_to_prob(r_SD_D_BSC  * hr_SD_D_C1 , t = cycle_length)
p_PD_D_C1  <- rate_to_prob(r_PD_D_BSC  * hr_PD_D_C1 , t = cycle_length)
## Compute transition probabilities from rates under Intervention
p_SD_PD_Int <- rate_to_prob(r_SD_PD_BSC * hr_SD_PD_Int, t = cycle_length)
p_SD_D_Int  <- rate_to_prob(r_SD_D_BSC  * hr_SD_D_Int , t = cycle_length)
p_PD_D_Int  <- rate_to_prob(r_PD_D_BSC  * hr_PD_D_Int , t = cycle_length)

####################### Construct state-transition models ######################
## Initial state vector
# All starting healthy
v_m_init <- c(SD = 1, PD = 0, D = 0) # initial state vector
v_m_init

## Initialize cohort trace for BSC
m_M <- matrix(NA, 
              nrow     = (n_cycles + 1), ncol = n_states, 
              dimnames = list(0:n_cycles, v_names_states))
# Store the initial state vector in the first row of the cohort trace
m_M[1, ] <- v_m_init
## Initialize cohort trace for Comparator 1 and Intervention
# Structure and initial states are the same as for BSC
m_M_C1   <- m_M # Comparator 1
m_M_Int  <- m_M # Intervention

## Initialize transition probability matrix for strategy BSC
# all transitions to a non-death state are assumed to be conditional on survival 
m_P <- matrix(NA, 
              nrow     = n_states, 
              ncol     = n_states, 
              dimnames = list(v_names_states, 
                              v_names_states)) # define row and column names
## Initialize transition probability matrix for Comparator 1 and Intervention 
## as a copy of BSC's
m_P_C1  <- m_P
m_P_Int <- m_P

### Fill in matrices
## Under BSC
# From H
m_P["SD", "SD"] <- (1 - p_SD_PD_BSC - p_SD_D_BSC)
m_P["SD", "PD"] <- p_SD_PD_BSC
m_P["SD", "D"]  <- p_SD_D_BSC
# From S1
m_P["PD", "SD"] <- 0
m_P["PD", "PD"] <- (1 - p_PD_D_BSC)
m_P["PD", "D"]  <- p_PD_D_BSC
# From D
m_P["D", "SD"] <- 0
m_P["D", "PD"] <- 0
m_P["D", "D"]  <- 1

## Under Comparator 1
# From H
m_P_C1["SD", "SD"] <- (1 - p_SD_PD_C1 - p_SD_D_C1)
m_P_C1["SD", "PD"] <- p_SD_PD_C1
m_P_C1["SD", "D"]  <- p_SD_D_C1
# From S1
m_P_C1["PD", "SD"] <- 0
m_P_C1["PD", "PD"] <- (1 - p_PD_D_C1)
m_P_C1["PD", "D"]  <- p_PD_D_C1
# From D
m_P_C1["D", "SD"] <- 0
m_P_C1["D", "PD"] <- 0
m_P_C1["D", "D"]  <- 1

## Under Intervention
# From H
m_P_Int["SD", "SD"] <- (1 - p_SD_PD_Int - p_SD_D_Int)
m_P_Int["SD", "PD"] <- p_SD_PD_Int
m_P_Int["SD", "D"]  <- p_SD_D_Int
# From S1
m_P_Int["PD", "SD"] <- 0
m_P_Int["PD", "PD"] <- (1 - p_PD_D_Int)
m_P_Int["PD", "D"]  <-      p_PD_D_Int
# From D
m_P_Int["D", "SD"] <- 0
m_P_Int["D", "PD"] <- 0
m_P_Int["D", "D"]  <- 1

### Check if transition probability matrices are valid
## Check that transition probabilities are [0, 1]
check_transition_probability(m_P,     verbose = TRUE)
check_transition_probability(m_P_C1,  verbose = TRUE)
check_transition_probability(m_P_Int, verbose = TRUE)
## Check that all rows sum to 1
check_sum_of_transition_array(m_P,     n_states = n_states, n_cycles = n_cycles, verbose = TRUE)
check_sum_of_transition_array(m_P_C1,  n_states = n_states, n_cycles = n_cycles, verbose = TRUE)
check_sum_of_transition_array(m_P_Int, n_states = n_states, n_cycles = n_cycles, verbose = TRUE)

#### Run Markov model ####
# Iterative solution of time-independent cSTM
for(t in 1:n_cycles){
  # For SoC
  m_M[t + 1, ]     <- m_M[t, ]     %*% m_P
  # For Comparator 1
  m_M_C1[t + 1, ]  <- m_M_C1[t, ]  %*% m_P_C1
  # For Intervention
  m_M_Int[t + 1, ] <- m_M_Int[t, ] %*% m_P_Int
}

## Store the cohort traces in a list
l_m_M <- list(m_M,
              m_M_C1,
              m_M_Int)
names(l_m_M) <- v_names_str

#### Plot Outputs ####
### Plot the cohort trace for strategy BSC
plot_trace(m_M)

#### State Rewards ####
## Vector of state utilities under strategy BSC
v_u_BSC    <- c(SD = u_SD_BSC, 
                PD = u_PD_BSC,
                D  = u_D)
## Vector of state costs under strategy BSC
v_c_BSC    <- c(SD = c_SD_MedRes_BSC + c_SD_DrugAcq_BSC + c_SD_DrugAdm_BSC, 
                PD = c_PD_MedRes_BSC,
                D  = c_D)
## Vector of state utilities under Comparator 1
v_u_C1    <- c(SD = u_SD_C1, 
               PD = u_PD_C1,
               D  = u_D)
## Vector of state costs under strategy A
v_c_C1   <- c( SD = c_SD_MedRes_C1 + c_SD_DrugAcq_C1 + c_SD_DrugAdm_C1, 
               PD = c_PD_MedRes_C1,
               D  = c_D)
## Vector of state utilities under Intervention
v_u_Int   <- c(SD = u_SD_Int, 
               PD = u_PD_Int, 
               D  = u_D)
## Vector of state costs under strategy B
v_c_Int   <- c(SD = c_SD_MedRes_Int + c_SD_DrugAcq_Int + c_SD_DrugAdm_Int, 
               PD = c_PD_MedRes_Int,
               D  = c_D)

## Store the vectors of state utilities for each strategy in a list 
l_u   <- list(v_u_BSC,
              v_u_C1,
              v_u_Int)
## Store the vectors of state cost for each strategy in a list 
l_c   <- list(v_c_BSC,
              v_c_C1,
              v_c_Int)

# assign strategy names to matching items in the lists
names(l_u) <- names(l_c) <- v_names_str
  
## create empty vectors to store total utilities and costs 
v_tot_qaly <- v_tot_cost <- vector(mode = "numeric", length = n_str)
names(v_tot_qaly) <- names(v_tot_cost) <- v_names_str

#### Loop through each strategy and calculate total utilities and costs ####
for (i in 1:n_str) {
  v_u_str <- l_u[[i]]   # select the vector of state utilities for the i-th strategy
  v_c_str <- l_c[[i]]   # select the vector of state costs for the i-th strategy
  
  #### Expected QALYs and costs per cycle ####
  ### Vector of QALYs and Costs
  ## Apply state rewards ###
  v_qaly_str <- l_m_M[[i]] %*% v_u_str # sum the utilities of all states for each cycle
  v_cost_str <- l_m_M[[i]] %*% v_c_str # sum the costs of all states for each cycle
  
  #### Discounted total expected QALYs and Costs per strategy and apply half-cycle correction if applicable ####
  ## QALYs
  v_tot_qaly[i] <- t(v_qaly_str) %*% (v_dwe * v_wcc) + v_du_AE[i] # Subtract AE disutility (one-off, applied in the 1st cycle)
  ## Costs
  v_tot_cost[i] <- t(v_cost_str) %*% (v_dwc * v_wcc) + v_c_AE[i] # Add AE costs (one-off, applied in the 1st cycle)
}

########################## Cost-effectiveness analysis #######################
### Calculate incremental cost-effectiveness ratios (ICERs)
df_cea <- calculate_icers(cost       = v_tot_cost, 
                          effect     = v_tot_qaly,
                          strategies = v_names_str)
df_cea

### Create CEA table in proper format
table_cea <- format_table_cea(df_cea)
table_cea

### CEA frontier
plot(df_cea, label = "all") +
  expand_limits(x = max(table_cea$QALYs) + 0.5) 
