# this function generates parameters for the microsimulation model
# STEPS:
# 1) define all initial input in a list 
# 2) modify this input based on shiny values
# 3) process the inital parameter input to generate the final model input

generate_params <- function(input_list){
 
# 1) Define all initial parameter input in a list   
init_params <- list(
  v_n         = c("healthy", "sick", "dead"),  # vector with state names
  n_t         = 60,                            # number of cycles
  n_i         = 10000,                         # number of individuals
  d_e         = 0.03,                        # equal discount of costs and QALYs by 3% 
  d_c         = 0.03,
  p_HS        = 0.05,    # probability healthy -> sick
  p_HD_female = 0.0382,  # probability health -> dead when female
  p_HD_male   = 0.0463,  # probability health -> dead when male
  p_SD1_5     = 0.1:0.5, # probability sick   -> dead in cycles 1-5
  p_SD6       = 0.7,     # probability sick   -> dead in cycle 6
  # Costs inputs
  c_H         = 1500,  # cost of one cycle in healthy state
  c_S         = 5000,  # cost of one cycle in sick state
  c_D         = 0,
  
  # utility inputs
  u_H         = 1,    # utility when healthy 
  u_S         = 0.85,  # utility when sick 
  u_D         = 0     # utility when dead
)

# 2) replace initial values with Shiny parameter values
init_params <- modifyList(init_params,input_list)


# 2) process parameter values to generate the final model input
process_params<- with(init_params,{
  list(
n_states    = length(v_n),                   # number of states
# probability to die in sick state by cycle of being sick
p_SD        = c(p_SD1_5, rep(p_SD6, n_t - 5)) ,

# calculate discount weights for costs for each cycle based on discount rate d_c
v_dwc       = 1 / (1 + d_e) ^ (0:n_t) ,
# calculate discount weights for effectiveness for each cycle based on discount rate d_e
v_dwe       = 1 / (1 + d_c) ^ (0:n_t) ,
#### Deterministic analysis ####

# Transition probabilities 
# (all non-dead probabilities are conditional on survival)
m_p_HD      = data.frame(Sex = c("Female", "Male"), p_HD = c(p_HD_female, p_HD_male)),

# Specify the initial health state of the individuals 
# everyone begins in the healthy state (in this example)
v_M_init    = rep("healthy", times = n_i)   ,
v_Ts_init   = rep(0, n_i)  # a vector with the time of being sick at the start of the model  
)})

# return the finla model inputs
return(c(init_params,process_params))
}