#' @title Static Model Log Likelihood
#' @description  compute log marginal likelihood of observing y in the static model
#' @param parameters a list containing beta and rho
#' @param model_config a list containing population size and agent features
#' @param method can be from the list c('tp','exact','mc')
#' @param particle_config if method == 'mc', then use this argument and it must contain num_particles
#' @export
#' 

static_loglikelihood_marginal <- function(y, parameters, model_config,  method = 'tp', particle_config = NULL){
  rho_ <- parameters$rho
  alpha_ <- get_rates_from_features(beta = parameters$beta, features = model_config$features)
  ## error messages
  method <- match.arg(method, choices = c('tp', 'exact', 'mc'))
  if (method == 'mc' & is.null(particle_config$num_particles)){
    stop("please specify the number of particles.")
  }
  if (method == 'tp'){
    loglik <- lw.logsum(logdtranspoisson_approx(alpha = alpha_,  y = y:length(alpha_)) + dbinom(x = y, size = y:length(alpha_), prob = rho_, log = TRUE))
  }
  if (method == 'exact'){
    loglik <- lw.logsum(logdpoisbinom_cpp(alpha = alpha_)[(y + 1):(length(alpha_) + 1)] + dbinom(x = y, size = y:length(alpha_), prob = rho_, log = TRUE))
  }
  if (method == 'mc'){
    particles <- replicate(n = particle_config$num_particles, 
                           expr = sum(runif(length(alpha_)) < alpha_))
    logweights <- dbinom(x = y, size = particles, prob = rho_, log = TRUE)
    loglik <- lw.logmean(logweights)
  }
  return(loglik)
}
