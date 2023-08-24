# Make function
# A function in R is an object containing multiple interrelated statements that are run together in a predefined order every time the function is called. 

# Packages contain many function
# You can write them yourself
# Main purpose of creating a user-defined function is to optimize our program, avoid the repetition of the same block of code used for a specific task that is frequently performed in a particular project

# A good practice is creating a function whenever we're supposed to run a certain set of commands more than twice.


# Simple examples of function in R
x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
mean(x)
min(x)
max(x)
length(x)



c_S
c_D

calculate_BMI <- function(weight, height){
  # Argument:
  # weight: weight of an individual in kg
  # height: height of an individual in meters
  
 # browser()
  BMI <- weight / (height^2)
  BMI_rounded <- round(BMI, digits = 0 )
  return(BMI_rounded)
  
}

calculate_BMI(height = 1.68, weight = 70)











calculate_BMI <- function(height, weight){
  #Arguments
  # height in meters
  # weight in kg
  # Returns
  # the BMI of a person 
  
  BMI <- weight / (height^2)
    
  return(BMI)
  
}


