### Load packages
library(darthtools)
library(dampack) 
library(doParallel)

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
                 "Intervention (Int") 
n_str       <- length(v_names_str)        # number of strategies

# Within-cycle correction (WCC) using Simpson's 1/3 rule
v_wcc <- darthtools::gen_wcc(n_cycles = n_cycles, 
                             method = "Simpson1/3") # vector of wcc

## Transition rates (per cycle), hazard ratios 
r_SD_PD <- 1/200 # constant rate of progressing
r_SD_D  <- 1/300 # constant rate of dying from SD
r_PD_D  <- 1/100 # constant rate of dying from PD

## Effectiveness of interventions as hazard ratios (HR)
## Comparator 1 vs. BSC
hr_SD_PD_C1 <- 0.80
hr_SD_D_C1  <- 0.90
hr_PD_D_C1  <- 0.95
## Intervention vs. BSC
hr_SD_PD_Int <- 0.75
hr_SD_D_Int  <- 0.85
hr_PD_D_Int  <- 0.90

### State rewards
## Costs
c_H    <- 2000  # cost of remaining one cycle in Healthy 
c_S1   <- 4000  # cost of remaining one cycle in Sick 
c_S2   <- 15000 # cost of remaining one cycle in Sicker 
c_D    <- 0     # cost of being dead (per cycle)
c_trtA <- 12000 # cost of treatment A
c_trtB <- 13000 # cost of treatment B
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
du_SD_BSC <- -0.100 # under BSC
du_SD_C1  <- -0.025 # under C1
du_SD_Int <- -0.050 # under Int

# Discount weight for costs and effects
v_dwc  <- 1 / ((1 + d_e) ^ (0:n_cycles))
v_dwe  <- 1 / ((1 + d_c) ^ (0:n_cycles))

### Process model inputs
## Compute transition probabilities from rates
p_SD_PD <- rate_to_prob(r_SD_PD)
p_SD_D  <- rate_to_prob(r_SD_D)
p_PD_D  <- rate_to_prob(r_PD_D)

####################### Construct state-transition models ######################
## Initial state vector
# All starting healthy
v_m_init <- c(SD = 1, PD = 0, D = 0) # initial state vector
v_m_init
