#' @title check  model configurations
#' @description this functionc checks if model_config is eligible as input for simulations, particle filters and gibbs samplers.
#' @param model_config a list that must contain: TODO
#' @export
check_model_config <- function(model_config){
  if(is.null(model_config$N)) stop("please specify population size");
  if(is.null(model_config$network_type)) stop("please specify network type. it can be full or general");
  match.arg(model_config$network_type, choices = c("full", "general"));
  if(model_config$network_type == "general"){
    if(is.null(model_config$adjacency_matrix_b)) stop("please specify the normalized adjacency matrix");
    if(any(colSums(model_config$adjacency_matrix_b) != 1)){
      warning("the adjacency matrix is not normalized");
      # adjacency_matrix_b <- apply(adjacency_matrix_b, 2, function(col) col / sum(col));
    } 
  }
  # if(is.null(model_config$features)) stop("please specify agent features");
  # if(dim(model_config$features)[1] != model_config$N) stop("feature matrix is not the right size");
  if(is.null(model_config$alpha0)) stop("please specify alpha0");
  if(is.null(model_config$lambda)) stop("lambda is null");
  if(is.null(model_config$gamma)) stop("gamma is null");
  if(is.null(model_config$rho)) stop("rho is null");
  model_config;
}


#' @title check particle configurations
#' @description  this function checks if particle_config is eligible as input for particle filters
#' @param particle_config a list that contains that must contain: TODO
#' @export
check_particle_config <- function(particle_config){
  if(is.null(particle_config$num_particles)) stop("please specify number of particles");
  if(is.null(particle_config$ess_threshold)) stop("please give ess threshold for resampling");
  if(is.null(particle_config$save_genealogy)) particle_config$save_genealogy <- FALSE;
  if(is.null(particle_config$clock)) particle_config$clock <- FALSE;
  if(is.null(particle_config$verbose)) particle_config$verbose <- FALSE;
  if(is.null(particle_config$exact)){
   stop("using exact Poisbinom densities or exact CondBern samplers?");
    particle_config$exact <- T;
  }
  particle_config;
}