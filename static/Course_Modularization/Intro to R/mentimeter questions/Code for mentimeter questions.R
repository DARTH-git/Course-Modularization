
my_vector <- c(1:10)
length(my_vector)

# data type questions
my_matrix <- matrix(1:10, 
                    nrow = 5, 
                    ncol = 2)
length(my_matrix)
dim(my_matrix)



# Data manipulation question
##1
month <- list(
  x = c("Jan", "Feb", "March"),
  y = c("Aug", "Sep", "Oct")
)

month$x[1]


##2 
my_matrix <- matrix(11:26, 
                    nrow = 4, 
                    ncol = 4)

my_matrix
my_matrix[, 3]


# Loops and functions

# loops
x <- 2020
if(x < 2000){
  price <- 1000 + 20
} else {
  price <- 1000 + 3000
}

price


# functions 

v_lenght <- c(152, 180, 199, 176 ) # length in cm
v_weight <- c(60, 80, 95, 75)      # weight in kg

calculate_BMI <- function(lenght, weight){
  
  BMI <- weight / (lenght/100)^2 # Calculate BMI
  BMI <- round(BMI)
  return(BMI) # return the BMI
}

calculate_BMI(lenght = v_lenght, weight = v_weight)

