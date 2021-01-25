#' Gibbs sampler for SIS model
#' @description update the agent states for the SIS model given observations such that the complete likelihood is not zero.
#' @param xx a matrix of size model_config$N by length(y), can be initialized from sis_init_hidden_popobservation.R 
#' @param y a list of population observations 
#' @param model_config a list that contains parameters and network structure
#' @return hiden states
#' @export

sis_xx_gibbs_singlesite_kernel <- function(xx, y , model_config){
  num_observations <- length(y)
  t <- 1
  for (n in 1: model_config$N){
    xxt1 <- xxt0 <- rep(NA, model_config$N)
    xxt1[n] <- TRUE
    xxt0[n] <- FALSE 
    xxt1[-n] <- xxt0[-n] <- xx[-n,t]
    sum_x_notn <- sum(xxt0)
    a1 <- a0 <- rep(NA, model_config$N)
    a1 <- sis_get_alpha(xxt1, model_config)
    a0 <- sis_get_alpha(xxt0, model_config)
    logp <- rep(NA, 2) ## logp[1] = logp(xn = 1) and logp[2] = logp(xn = 0)
    ## complete likelihood contains ft from previous step to current step, current step to next step and current observation density
    logp[1] <- log(model_config$alpha0[n]) + 
      logdbern_sum_cpp(a1, xx[,t+1]) + 
      dbinom(x = y[t], size = sum_x_notn +1 , prob = model_config$rho, log = TRUE)
    logp[2] <- log(1 - model_config$alpha0[n]) + 
      logdbern_sum_cpp(a0, xx[,t+1]) + 
      dbinom(x = y[t] , size = sum_x_notn, prob = model_config$rho, log = TRUE)
    logp <- logp - max(logp)
    if(runif(1) < exp(logp[1]) / sum(exp(logp))){
      xx[n, t] <- TRUE
    }else{
      xx[n, t] <- FALSE 
    }
  }
  for (t in 2 : (num_observations - 1 )){
    previous_alpha <- sis_get_alpha(xx[, t - 1], model_config)
    for (n in 1 : model_config$N){
      xxt1 <- xxt0 <- rep(NA, model_config$N)
      xxt1[n] <- 1
      xxt0[n] <- 0 
      xxt1[-n] <- xxt0[-n] <- xx[-n,t]
      sum_x_notn <- sum(xxt0)
      a1 <- a0 <- rep(NA, model_config$N)
      a1 <- sis_get_alpha(xxt1, model_config)
      a0 <- sis_get_alpha(xxt0, model_config)
      logp <- rep(NA, 2) ## logp[1] = logp(xn = 1) and logp[2] = logp(xn = 0)
      ## complete likelihood contains ft from previous step to current step, current step to next step and current observation density
      logp[1] <- log(previous_alpha[n]) + 
        logdbern_sum_cpp(a1, xx[,t+1]) + 
        dbinom(x = y[t], size = sum_x_notn + 1, prob = model_config$rho, log = TRUE)
      logp[2] <- log(1 - previous_alpha[n]) + 
        logdbern_sum_cpp(a0, xx[,t+1]) + 
        dbinom(x = y[t] , size = sum_x_notn , prob = model_config$rho, log = TRUE)
      logp <- logp - max(logp)
      if(runif(1) < exp(logp[1]) / sum(exp(logp))){
        xx[n, t] <- 1
      }else{
        xx[n, t] <- 0 
      }
    }
  }
  t <- num_observations 
  previous_alpha <- sis_get_alpha(xx[, t - 1], model_config)
  for (n in 1 : model_config$N){
    xxt1 <- xxt0 <- rep(NA, model_config$N)
    xxt1[n] <- 1
    xxt0[n] <- 0 
    xxt1[-n] <- xxt0[-n] <- xx[-n,t]
    a1 <- a0 <- rep(NA, model_config$N)
    a1 <- sis_get_alpha(xxt1, model_config)
    a0 <- sis_get_alpha(xxt0, model_config)
    sum_x_notn <- sum(xxt1[-n])
    logp <- rep(NA, 2) ## logp[1] = logp(xn = 1) and logp[2] = logp(xn = 0)
    ## complete likelihood contains ft from previous step to current step and current observation density
    logp[1] <- log(previous_alpha[n])  + dbinom(y[t], sum_x_notn + 1, model_config$rho, TRUE) 
    logp[2] <- log(1 - previous_alpha[n]) + dbinom(y[t], sum_x_notn, model_config$rho, TRUE) 
    logp <- logp - max(logp)
    if(runif(1) < exp(logp[1]) / sum(exp(logp))){
      xx[n, t] <- 1
    }else{
      xx[n, t] <- 0 
    }
  }
  return(xx)
}
