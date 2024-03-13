# imagine 5 years of follow up on a trial
# 30 of 100 people become sick
# 30 of 100 people die
# 40 people are left healthy

p_HD_5 <- 0.3
p_HS_5 <- 0.3



cycle_length    <- 1                            # cycle length equal to one year (use 1/12 for monthly)
n_cycles        <- 5/cycle_length                   # number of cycles
v_names_cycles  <- paste("cycle", 0:n_cycles)    # cycle names
v_names_states  <- c("H", "S", "D")  # state names
n_states        <- length(v_names_states)        # number of health states 

#

p_HD <- prob_to_prob(p_HD_5,1/n_cycles)
p_HS <- prob_to_prob(p_HS_5,1/n_cycles)


v_m_init <- c("H" = 1, "S" = 0, "D" = 0)  

m_M_SoC <- matrix(0, 
                  nrow = (n_cycles + 1), ncol = n_states, 
                  dimnames = list(v_names_cycles, v_names_states))
# Store the initial state vector in the first row of the cohort trace
m_M_SoC[1, ] <- v_m_init


a_P_SoC <- array(0,  # Create 3-D array
                 dim = c(n_states, n_states, n_cycles),
                 dimnames = list(v_names_states, v_names_states, 
                                 v_names_cycles[-length(v_names_cycles)])) # name the dimensions of the array 

### Fill in array
## Standard of Care
# from Healthy
a_P_SoC["H", "H", ] <- (1 - p_HD - p_HS)
a_P_SoC["H", "S", ] <-  p_HS
a_P_SoC["H", "D", ] <-  p_HD

# from Sick
a_P_SoC["S", "S", ] <- 1
a_P_SoC["D", "D", ] <- 1
# Iterative solution of age-dependent cSTM
for(t in 1:n_cycles){
  ## Fill in cohort trace
  # For SoC
  m_M_SoC[t + 1, ]  <- m_M_SoC[t, ]  %*% a_P_SoC[, , t]
  # For strategy A
}

plot_trace(m_M_SoC)


# what if we were to work with rates

t_S <- runif(30,1,n_cycles)
t_D <- runif(30,1,n_cycles)
t_H <- rep(5,40)

personyears <- (sum(t_H)+sum(t_D)+sum(t_S))
r_HD <- 30 / personyears 
#similarly for disease
r_HS <- 30 / personyears



# #or 
#library(flexsurv)
# dur <- c(t_H,t_S, t_D)
# sick  <- c(rep(0,40),rep(1,30),rep(0,30))
# death <- c(rep(0,40),rep(0,30), rep(1,30)) 
# fit_HD<-flexsurvreg(Surv(dur, death )~1, dist = "exp")
# fit_HS<-flexsurvreg(Surv(dur, sick)~1, dist = "exp")
# r_HD <- fit_HD$res[1]
# r_HS <- fit_HS$res[1]

p_HD<- prob_to_rate(r_HD)
p_HS<- prob_to_rate(r_HS)



m_M_SoC2 <- matrix(0, 
                  nrow = (n_cycles + 1), ncol = n_states, 
                  dimnames = list(v_names_cycles, v_names_states))
# Store the initial state vector in the first row of the cohort trace
m_M_SoC2[1, ] <- v_m_init


a_P_SoC2 <- array(0,  # Create 3-D array
                 dim = c(n_states, n_states, n_cycles),
                 dimnames = list(v_names_states, v_names_states, 
                                 v_names_cycles[-length(v_names_cycles)])) # name the dimensions of the array 

### Fill in array
## Standard of Care
# from Healthy
a_P_SoC2["H", "H", ] <- (1 - p_HD - p_HS)
a_P_SoC2["H", "S", ] <-  p_HS
a_P_SoC2["H", "D", ] <-  p_HD

# from Sick
a_P_SoC2["S", "S", ] <- 1
a_P_SoC2["D", "D", ] <- 1
# Iterative solution of age-dependent cSTM
for(t in 1:n_cycles){
  ## Fill in cohort trace
  # For SoC
  m_M_SoC2[t + 1, ]  <- m_M_SoC2[t, ]  %*% a_P_SoC2[, , t]
  
}

plot_trace(m_M_SoC2)

m_M_SoC2


