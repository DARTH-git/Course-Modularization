### This R script illustrates how you can generate synthetic datasets in R
# The generated synthetic dataset will not contain any real values from the original dataset
# The DARTH workgroup 2020


### Step 1: Read in your source dataset
# you must specify the full file path to your dataset
# specify any other arguments based on the structure of your dataset
# i.e. source_data <- read.csv("full file path to your dataset", header=T)

# From now we will use this sample dataset to illustrate the next steps
# You can skip specifying the file path if you put this script and your source dataset in the same location and
# Click: Session -> Set Working Directory -> To Source File Location in the top bar of your RStudio
source_data <- read.csv('sample_data.csv', header = T)


### Step 2: Extract the identification variable (ID variable, this could be named differently in your dataset)
if (!require("dplyr")) install.packages("dplyr")
library(dplyr)
IDs <- source_data$ID
# remove the ID variable from the source dataset
source_data <- dplyr::select(source_data, -ID)

# Any other data cleaning/manipulation procedure should be done at this stage
# Most importantly, only keep the variables you will need, remove the ones you don't need 
# (no need to synthesize the variables you don't need)


### Step 3: Work on the first variable
# if the 1st variable/column in your dataset (after you remove the ID variable) is not numeric/continuous, 
# re-order the columns of your dataset so that the first column is a continuous variable. 
# If the 1st variable is already continuous, skip this above steps

# Now, because our synthetic data generation method synthesizes the first variable by taking a bootstrap sample
# from its values in the original data, this might draw some concerns with data privacy issues
# To eliminate this concern, do the following (so that we do not have ANY original data values)

# 1) Extract the first variable from your dataset
first_var  <- source_data[, 1]

# 2) Sample from the marginal distribution of the first variable
# Take a random sample of the first variable from a normal distribution using sample mean and sample standard deviation
# The first argument is the sample size, it should be the number of rows in your source dataset
set.seed(123) # set the seed to ensure reproducible results
first_var_fitted <- rnorm(nrow(source_data), mean=mean(first_var), sd=sd(first_var))

# 3) Replace the 1st variable in your source dataset with the fitted first variable from the previous step
source_data[, 1] <- first_var_fitted


### Step 4: Generate the synthetic dataset using functions in `synthpop`
# Install and load the `synthpop` package
if (!require("synthpop")) install.packages("synthpop")
library(synthpop)

# Use the syn() function to synthesize data
# the 1st argument is your pre-processed original dataset before step 3
# the 2nd argument is the method of synthesis - use 'parametric'
# the 3rd argument is the random seed, you should set it for reproducible results
syn_param <- syn(source_data, method = 'parametric', seed = 123)  # synthesize the data (process)
# obtain the actual synthetic dataset 
synthetic_data <- syn_param$syn # obtain the actual synthetic dataset 

# Put the ID variable back in
synthetic_data <- cbind(IDs, synthetic_data)

# Rename your ID variable to what it was originally in your source dataset
colnames(synthetic_data)[1] <- 'ID'


### Step 6: Export the synthetic dataset as a .csv file
# 1st argument is the synthetic data
# 2nd argument is the name of the .csv file (you specify)
write.csv(synthetic_data, file='synthetic_sample_data.csv')
