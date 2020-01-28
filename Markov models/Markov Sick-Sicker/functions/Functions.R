#------------------------------------------------------------------------------#
#### R function to extract the parameters of a beta distribution            ####
####                from mean and st. deviation                             ####
#------------------------------------------------------------------------------#
#' @param m mean 
#' @param s standard deviation
#' 
betaPar <- function(m, s) 
{
  a <- m * ((m * (1 - m) / s ^ 2) - 1)
  b <- (1 - m) * ((m * (1 - m) / s ^ 2) - 1)
  list(a = a, b = b)
}

#------------------------------------------------------------------------------#
#### R function to extract the parameters of a gamma distribution           ####
####                   from mean and st. deviation                          ####
#------------------------------------------------------------------------------#
#' @param m mean 
#' @param s standard deviation
#' 
gammaPar <- function(m, s) {   
  # m: mean  
  # s: standard deviation 
  shape <- m ^ 2 / s ^ 2
  scale <- s ^ 2 / m
  list(shape = shape, scale = scale)
}
