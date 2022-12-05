library('OpenTree')
library('remotes')


# Decision Tree Setup
setwd("/Users/hjalal/files/4_Teaching/CDC2022/PartII/Day1")


# The first time create the tree 
#create_tree(file_name="Mosquito_example",dir_name=getwd())

# Then open the tree
open_tree(file_name="Mosquito_example",dir_name=getwd())

#Parameter Names

v_names_str<-c("Do nothing", "Spray", "Test and Spray") #Names of Strategies
n_str


#Evaluate Decision Tree
c_spray <- 1500
c_test <- 4000
c_infection_survived <- 10000
c_infection_died <- 20000
c_toxicity_death <- 5000
df_tree<-evaluate_model("Mosquito_example", n_payoffs = 2)
df_tree


# vector of total cost and QALYs
v_total_qaly <- v_total_cost <- vector(mode = "numeric", length = n_str)

# Calculate total costs and QALYs for each strategy 
for (i in 1:n_str) {
  v_total_qaly[i] <- df_tree[[i]]$prob %*% df_tree[[i]]$payoff1
  v_total_cost[i] <- df_tree[[i]]$prob %*% df_tree[[i]]$payoff2
}

# calculate vector of nmb
wtp <- 100000
v_nmb        <- v_total_qaly * wtp - v_total_cost                      

df_output <- data.frame(Strategy =  v_names_str,
                        Cost     =  v_total_cost,
                        Effect   =  v_total_qaly,
                        NMB      =  v_nmb)

# model output
df_output

## 05 Cost-Effectiveness Analysis

# create the transition probability matrix for NO treatment
library(dampack)
decision_tree_rockclimber_cea  <- calculate_icers(cost       = df_output$Cost,
                                                  effect     = df_output$Effect,
                                                  strategies = df_output$Strategy)
decision_tree_rockclimber_cea


## 05.1 Plot frontier of Decision Tree

plot(decision_tree_rockclimber_cea, effect_units = "QALYs", label="all")
