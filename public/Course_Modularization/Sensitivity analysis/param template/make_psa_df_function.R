
# Make a function to make the PSA data set 
make_psa_df <- function(df_param, n_iter, seed = 123){
  # Arguments:
  ## df_param: a dataframe with the parameters
  ## n_iter:   the number of PSA iterations
  ## seed :    seed to be able to reproduce the results, default = 123
  # Return:
  ## param_psa: dataframe with estimated PSA parameters
  
  set.seed(seed)    # set the seed
  l_shapes <- list()

  # check if all parameters have a unique name
  if(length(df_param$parameter) != length(unique(df_param$parameter))){
    print("the data contains with at least two parameters with similar names")
  }
  
  # create a matrix structure to store the psa parameters
  v_names_param         <- df_param$parameter # extract the parameter names
  m_param_psa           <- matrix(data = NA, 
                                  nrow = n_iter, 
                                  ncol = length(v_names_param),
                                  dimnames = list(paste("iter", 1:n_iter),
                                               v_names_param)) # make matrix


  for(p in v_names_param){    # loop over all parameters
    # select the row with information about the parameter
    v_param <- df_param[df_param$parameter == p, ]
    # name the parameter
    v_shapes        <- v_param[, c("shape1",      "shape2",      "shape3")]
    names(v_shapes) <- v_param[, c("type.shape1", "type.shape2", "type.shape3")]
    v_shapes        <- v_shapes[c(!is.na(v_shapes))] # remove NA's 

    param <- cbind(v_param[, c("parameter", "unit", "distribution")], 
                   v_shapes)
    
    # Beta distribution
    if(param$distribution == "beta"){
      # if mean and sigma -> get the shapes and add
      if(c("sigma") %in% names(param)){
        l_shapes   <-  with(param, 
                                dampack::beta_params(mean  = mean,
                                                     sigma = sigma))
      param$shape1 <- l_shapes$alpha
      param$shape2 <- l_shapes$beta } 
      
      # now we have the shapes sample the values from the distribution
      # store in the matrix 
      m_param_psa[, p] <- with(param, 
                               rbeta(n = n_iter,
                                shape1 = shape1,
                                shape2 = shape2))
      
      } else if(param$distribution == "gamma"){
      # if mean and sigma -> get the shapes and add
      if(c("sigma") %in% names(param)){
        l_shapes   <-  with(param, 
                            dampack::gamma_params(mu  = mean,
                                                 sigma = sigma))
        param$shape  <- as.numeric(l_shapes$shape)
        param$scale  <- as.numeric(l_shapes$scale) }
        
        m_param_psa[, p] <- with(param, 
                                 rgamma(n = n_iter,
                                        shape = shape,
                                        scale = scale))

      } else if (param$distribution == "triangle"){
        m_param_psa[, p] <- with(param, 
                                rtriangle(n = n_iter,
                                          a = lower_limit, 
                                          b = upper_limit, 
                                          c = mean))
        
      } else if (param$distribution == "normal"){
        if(c("upper_limit") %in% names(param)){
          
          param$sd <- with(param, 
                           (upper_limit - mean)/1.96) }
        
        m_param_psa[, p]  <- with(param, 
                                  rnorm(   n = n_iter ,
                                        mean = mean,
                                          sd = sd))
      } else if (param$distribution == "lognormal"){
        if(c("upper_limit") %in% names(param)){
        
        m_param_psa[, p] <- exp(with(param,
                                     rnorm(n = n_iter ,
                                           mean = log(mean),
                                           sd = (log(upper_limit) - log(mean)) / 1.96)))
        } # if the sd and mean are given
        m_param_psa[, p] <- exp(with(param,
                                     rnorm(n = n_iter ,
                                           mean = log(mean),
                                           sd = (log(sd)))))
        
        
      } else if (param$distribution == "uniform"){
        m_param_psa[, p] <- with(param,
                                 runif(n = n_iter,
                                     min = lower_limit,
                                     max = upper_limit))
        
        } else if(param$distribution == "weibull"){
          m_param_psa[, p] <- with(param,
                                   rweibull(n = n_iter,
                                         shape = shape,
                                         scale = scale))
    } else if (param$distribution == "NA"){
        m_param_psa[, p] <- rep(param$mean, n_iter) 
      }
    
  } # close loop for the parameters
  
  return(m_param_psa)  # Return the parameter values for the PSA runs
}  





