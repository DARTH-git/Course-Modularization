# Analysis
library(survival) # load the 
library(dplyr)
library(tidyr)
library(utils)
library(flexsurv)

# load data + describe info

# show good coding practice and a summary of the data set and the parameters that are used
# also how they relate to the sick-sicker case 

# X -> Sick
# x -> Dead


df_colon_wide <- colon %>%
  pivot_wider(names_from = etype, values_from = c(time, status))

# make a subset e-type equal to 1 ()
df_colon_wide <- df_colon_wide %>% subset(rx != "Obs")
# if status is 1 and etype is 1 = transition to recurrence 


df_colon_wide <- df_colon_wide %>% mutate(event_HS1  = if_else(status_1 == 1, 1, 0))  %>%
  mutate(event_S1D  = if_else(status_1 == 1 & status_2 == 1, 1, 0)) %>% 
  mutate(event_HD   = if_else(status_1 == 0 & status_2 == 1, 1, 0)) %>% 
  mutate(time_HS1  = time_1) %>% 
  mutate(time_S1D  = if_else(event_S1D == 1, time_2 - time_1, NA)) %>% 
                               mutate(time_HD   = if_else(event_HD == 1, time_2 , time_1)) 
                         


head(df_colon_wide)


# Data analysis to fit parameteric surival models



# sick to sick - based on recurrence 
head(data)

# What do we do for the treatment 
# Treatment effect used 

# sick to dead based on death

colon_wide <- colon %>%
  pivot_wider(names_from = etype, values_from = c(time, status))

# 