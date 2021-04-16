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
for (x in cars){
  print(x)
}

# loop through a matrix
y <- matrix(1:20, nrow = 5)
for (x in y){
  print(x)
}

## if else statements ========
if (TRUE) {
  print("Light")
} else {
  print("Dark")
}

# another format
ifelse(TRUE, "Light", "Dark")

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


# functions ============
x <- 1:10
mean(x)

my_mean <- function(x){
  y <- sum(x)/length(x)
  return(y)
}
my_mean(x)

