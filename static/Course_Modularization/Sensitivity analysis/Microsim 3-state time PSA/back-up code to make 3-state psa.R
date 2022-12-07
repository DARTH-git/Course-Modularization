

# 09 Probabilistic Sensitivity Analysis (PSA)

```{r}
# Function that generates random sample for PSA

gen_psa <- function(n_sim = 1000, seed = 071818){
  set.seed(seed) # set a seed to be able to reproduce the same results
  
  df_psa <- data.frame(
    # Transition probabilities (per cycle)
    
    # (all non-dead probabilities are conditional on survival)
    # probability of becoming sick when healthy, conditional on surviving
    p_HS_SoC  = rbeta(n_sim, shape1 = 24, shape2 = 450),      # under standard of care
    p_HS_trtA = rbeta(n_sim, shape1 = 170, shape2 = 4095),      # under treatment A
    p_HS_trtB = rbeta(n_sim, shape1 = 16, shape2 = 767),      # under treatment B    
    
    p_HD_female = 0.0382,  # probability healthy -> dead when female
    p_HD_male   = 0.0463,  # probability healthy -> dead when male
    
    ## State rewards
    # Costs
    c_H       = rgamma(n_sim, shape = 16, scale = 25),        # cost of one cycle in healthy state
    c_S       = rgamma(n_sim, shape = 100, scale = 10),       # cost of one cycle in sick state
    c_D       = 0,                                            # cost of one cycle in dead state
    c_trtA    = rgamma(n_sim, shape = 400, scale = 2),        # cost of treatment A (per cycle) in healthy state
    c_trtB    = rgamma(n_sim, shape = 501, scale = 3),        # cost of treatment B (per cycle) in healthy state
    
    # Utilities
    u_H       = rbeta(n_sim, shape1 =  1.5, shape2 = 0.0015), # utility when healthy 
    u_S       = rbeta(n_sim, shape1 = 49.5, shape2 = 49.5),   # utility when sick
    u_D       = 0                                             # utility when dead
  )
  
  return(df_psa)
}

# Try it
gen_psa(10) 


n_i <- 1000   # Decrease number of individuals since PSA takes a lot of time
n_sim <- 100  # Number of PSA simulations

# Dynamic characteristics 
# Specify the initial health state of the individuals 
# everyone begins in the healthy state (in this example)
v_M_init  <- rep("H", n_i)    # a vector with the initial health state for all individuals 

# Generate PSA input dataset
df_psa_input <- gen_psa(n_sim = n_sim)
# First six observations
head(df_psa_input)

## Histogram of parameters
# Make sure the Plots window is large enough to plot all the histograms
ggplot(melt(df_psa_input, variable.name = "Parameter"), aes(x = value)) +
  facet_wrap(~Parameter, scales = "free") +
  geom_histogram(aes(y = ..density..)) +
  theme_bw(base_size = 16)     

# Initialize dataframes with PSA output 
# Dataframe of costs
df_c <- as.data.frame(matrix(0, 
                             nrow = n_sim,
                             ncol = n_str))
colnames(df_c) <- v_names_str
# Dataframe of effectiveness
df_e <- as.data.frame(matrix(0, 
                             nrow = n_sim,
                             ncol = n_str))
colnames(df_e) <- v_names_str
```



## 09.2 Run microsimulation model on each parameter set of PSA input dataset

```{r}
start.time <- proc.time()

for(i in 1:n_sim){
  df_ce_psa <- calculate_ce_out(df_psa_input[i, ])
  df_c[i, ] <- df_ce_psa$Cost   # take the cost from the PSA run and store in df_c
  df_e[i, ] <- df_ce_psa$Effec  # take the cost from the PSA run and store in df_e
  # Display simulation progress
  if(i/(n_sim/10) == round(i/(n_sim/10),0)) { # display progress every 10%
    cat('\r', paste(' ', 'Overall progress: ', i/n_sim * 100, "% done", 
                    sep = " "))
  }
}
elapsed.time <-proc.time() - start.time

### Creae PSA object for dampack
l_psa <- make_psa_obj(cost          = df_c, 
                      effectiveness = df_e, 
                      parameters    = df_psa_input, 
                      strategies    = v_names_str)
```

## 09.3 Cost Effectiveness Analysis

Vector with willingness-to-pay (WTP) thresholds your considering and would like to have in your plot.

```{r}
v_wtp <- seq(0, 300000, by = 10000)
```

## 09.3.1 ICER

```{r}
# use dampack to calculate the ICER
calculate_icers(cost       = df_ce_psa$Cost,
                effect     = df_ce_psa$Effect,
                strategies = df_ce_psa$Strategy)
```

## 09.3.2 Cost-Effectiveness Acceptability Curves (CEAC) and Frontier (CEAF)

```{r}
out_ceaf <- ceac(v_wtp, l_psa)
plot(out_ceaf)
```

## 09.3.3 Cost-Effectiveness Scatter plot

```{r}
plot(l_psa)
```

## 09.4.4 Expected value of perfect information (EVPI)

```{r}
evpi <- calc_evpi(wtp = v_wtp, psa = l_psa)
# EVPI plot
plot(evpi, effect_units = "QALY")
```
