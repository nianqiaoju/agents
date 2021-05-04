#include <Rcpp.h>
#include "sis_cpp.h"
#include "sampling.h"
using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
double sis_apf_exact_full_cpp(const IntegerVector &y,
                        const NumericVector &alpha0,
                        const NumericVector &lambda,
                        const NumericVector &gamma,
                        const double rho,
                        const int num_particles,
                        const double threshold,
                        const int population_size){
  /// begin variable declarations
  int num_observations = y.size();
  // incremental log likelihoods
  NumericVector loglikelihood (num_observations);
  std::fill(loglikelihood.begin(), loglikelihood.end(), R_NegInf);
  // particles
  LogicalMatrix xts (num_particles, population_size);
  IntegerVector its(num_particles);
  NumericMatrix alphats (num_particles, population_size);
  // weights
  // logw: cumulative weights, used for resampling
  // logW: log of normalized logw
  // logweights: potential at the current step, use these to compute normalizing constant
  NumericVector logweights (num_particles);
  NumericVector logw (num_particles);
  NumericVector logW (num_particles);
  std::fill(logW.begin(), logW.end(), log(1 / (double) num_particles));
  IntegerVector ancestors (num_particles);
  for(int iparticle = 0; iparticle < num_particles; iparticle++){
    ancestors[iparticle] = iparticle;
  }
  NumericVector weights (num_particles);
  NumericVector runifs_p (num_particles);
  NumericVector runifs_n (population_size);
  NumericVector ess (num_observations);
  NumericVector logv0 (1 + population_size);
  NumericVector v0 (1 + population_size);
  NumericMatrix logvt (num_particles, 1 + population_size);
  // NumericMatrix vt (num_particles, 1 + population_size);
  /// end variable declarations

  // treat t = 0 differently
  // int t = 0;
  ess[0] = 1;
  logv0 = logdpoisbinom_cpp(alpha0);
  for(int it = 0; it <= population_size; it++){
    logv0[it] += R::dbinom(y[0],it,rho,true);
  }
  loglikelihood[0] = lw_logsum_cpp(logv0);
  // normalize to get conditional density of i0 | y0
  v0 = lw_normalize_cpp(logv0);
  // sample i0 | y0
  runifs_p = get_runifs(num_particles);
  its = multinomial_resampling_cpp(v0, num_particles, runifs_p);
  // sample x0 | i0
  for(int iparticle = 0; iparticle < num_particles; iparticle++){
    runifs_n = get_runifs(population_size);
    xts(iparticle, _) = idchecking_cpp(its[iparticle], alpha0, runifs_n);
  }
  
  // Rcout << "finished t = 0";
  // Rcout << ", loglikelihood = " << loglikelihood[0] << "\n";

  // // t > 0
  for(int t = 1; t < num_observations; t++){
    // Rcout << "t = " << t << " y[t] = " << y[t];
    // calculate p(x[t,k]  = 1 | x(t-1));
    for(int iparticle = 0; iparticle < num_particles; iparticle++){
      alphats(iparticle, _) = sis_get_alpha_full_cpp(xts(iparticle, _), lambda, gamma);
      logvt(iparticle,_) = logdpoisbinom_cpp(alphats(iparticle, _));
      for(int it = 0; it <= population_size; it++){
        logvt(iparticle, it) += R::dbinom(y[t],it,rho,true);
      }
      logweights[iparticle] = lw_logsum_cpp(logvt(iparticle,_));
      // vt(iparticle, _ ) = lw_normalize_cpp(logvt(iparticle, _));
    }
    loglikelihood[t]= lw_logsum_cpp(logweights) - log(num_particles);
    // Rcout << ", loglikelihood = " << loglikelihood[t];
    // resample: normalize the logweights
    logw = logW + logweights;
    weights = lw_normalize_cpp(logw);
    ess[t] = 1 / std::inner_product(weights.begin(), weights.end(), weights.begin(), 0.0) / num_particles; 
    // Rcout << ", ess = " << ess[t] << "\n";
    // the inner product computes sums of squares of weights
    if(ess[t] < threshold){
      // Rcout <<"RESAMPLE!\n";
      ancestors = multinomial_resampling_cpp(weights, num_particles, get_runifs(num_particles));
      std::fill(logW.begin(), logW.end(), - log((double) num_particles));
      std::fill(logw.begin(), logw.end(), 0);
    }else{
      for(int iparticle = 0; iparticle < num_particles; iparticle++){
        ancestors[iparticle] = iparticle;
      }
      logW = log(weights);
    }
    // Rcout << "log weights " << logweights << "\n";
    // sample it given yt and sample xt given it
    if(y[t] == population_size){ // no need to sample in this case
      std::fill(its.begin(), its.end(), population_size);
      for(int iparticle = 0; iparticle < num_particles; iparticle++){
        std::fill(xts.row(iparticle).begin(), xts.row(iparticle).end(), true);
      }
    }else{
      for(int iparticle = 0; iparticle < num_particles; iparticle++){
        GetRNGstate();
        double onerunif = unif_rand();
        PutRNGstate();
        its[iparticle] = multinomial_cpp(logvt(ancestors[iparticle],_), onerunif);
        xts(iparticle, _) = idchecking_cpp(its[iparticle], alphats(ancestors[iparticle],_), get_runifs(population_size));
      }
    }
  }
  return(sum(loglikelihood));
}
