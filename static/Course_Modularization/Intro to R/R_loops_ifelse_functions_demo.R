## for loop ========
# R favors vector math, but also supports simple for loops
for (x in 1:5){
  print(x)
}

# variable reset within loop?
for (x in 1:5){
  x <- 2
  print(x)
}

# break a loop
for (x in 1:5){
  print(x)
  if (x == 2) break
}

# loop through string
for (x in c("one", "two", "three")){
  print(x)
}

# loop through a dataframe
head(cars)
?cars
View(cars)
for (x in cars){
  print(x)
}

# loop through a matrix
y <- matrix(1:20, nrow = 5)
y
for (x in y){
  print(x)
}

## if else statements ========
x <- 2
if (x == 2) {
  print("Light")
} else {
  print("Dark")
}

# another format
ifelse(FALSE, "Light", "Dark")

# yet another format
if (TRUE) print("Light") else print("DARK")

# what happens if your condition has more than one element?
my_condition <- c(TRUE, FALSE)
if (my_condition) {
  print("Light")
} else {
  print("Dark")
}

ifelse(my_condition, "Light", "Dark")

# switch statement vs. case_when() ==== 
switch(2, "one", "two", "three")



# functions ============
x <- 1:10
mean(x)

# my first function 
my_mean <- function(x){
  y <- sum(x)/length(x)
  return(y)
}
my_mean(x)

# variable scope 
x <- 10
my_mean <- function(){
  print(x)
}
my_mean()
print(x)

# what will happen if the variable changes within a function?
x <- 10
my_mean <- function(){
  x <- 20
  print(x)
}
my_mean()
print(x)

# dealing with functions in libraries 
y <- matrix(1:20, nrow = 10)
# without loading the library, we can use library::function()
matrixStats::colSds(y)

# or we can load the library
library(matrixStats)
colSds(y)
# to examine the content of a function in a library, we can use
matrixStats::colSds
