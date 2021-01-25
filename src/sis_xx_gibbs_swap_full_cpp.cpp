#include <Rcpp.h>
#include "logsumexp.h"
#include "logdensities.h"
#include "sis_cpp.h"
using namespace Rcpp;

// this function performs N swaps at each time
// current implementation assumes a fully connected network
// it modifies the input xx directly
// [[Rcpp::export]]
void sis_xx_gibbs_swap_full_cpp(LogicalMatrix xx, 
                                const NumericVector &alpha0,
                                const NumericVector &lambda, 
                                const NumericVector &gamma){
  GetRNGstate();
  int num_swaps = xx.nrow();
  int num_agents = xx.nrow();
  int num_steps = xx.ncol();
  
  int infection_t;
  int i1;
  int i0;
  int i1_index;
  int i0_index;

  NumericVector alphaprev2t(num_agents);
  NumericVector logoddsprev2t(num_agents);
  NumericVector alphat2next(num_agents);
  
  double logpswap;
  double i1_temp; // use to store i1 during swap
  
  // set of indices n such that X^n = 0 
  vector<int> s0;
  // set of indices n such that X^n = 1 
  vector<int> s1;
  
  
  // first step updates x0
  int time_step = 0;
  // find current s0 and s1
  s0.clear(); 
  s1.clear();
  s0.reserve(num_agents); 
  s1.reserve(num_agents);
  for(int n = 0; n < num_agents; n++){
    if(xx(n, time_step)){
      s1.push_back(n);
    }else{
      s0.push_back(n);
    }
  }
  infection_t = sum(xx(_,time_step));
  if(infection_t >  0 && infection_t < num_agents){
    alphaprev2t = alpha0; 
    logoddsprev2t = log(alphaprev2t) - log(1 - alphaprev2t);
    alphat2next = sis_get_alpha_full_cpp(xx(_, time_step), lambda,gamma);
    for(int iswap = 0; iswap < num_swaps; iswap++){
      // choose i1 and i0, indices to swap
      i1_index = floor(unif_rand() * (infection_t));
      i0_index = floor(unif_rand() * (num_agents - infection_t));
      i1 = s1[i1_index];
      i0 = s0[i0_index];
      logpswap = logoddsprev2t[i0] - logoddsprev2t[i1] ;
      logpswap += ((xx(i1, time_step + 1))? log(lambda[i1] * infection_t / num_agents) : log(1 - lambda[i1] * infection_t / num_agents));
      logpswap += ((xx(i0, time_step + 1))? log(1 - gamma[i0]): log(gamma[i0]));
      logpswap -= ((xx(i1, time_step + 1))? log(alphat2next[i1]) : log(1 - alphat2next[i1]));
      logpswap -= ((xx(i0, time_step + 1))? log(alphat2next[i0]) : log(1 - alphat2next[i0]));
      if(log(unif_rand()) < logpswap){
        i1_temp = i1;
        s1[i1_index] = i0;
        s0[i0_index] = i1_temp;
        xx(i1, time_step) = false;
        xx(i0, time_step) = true;
        alphat2next[i0] = (1 - gamma[i0]);
        alphat2next[i1] = lambda[i1] * infection_t / num_agents;
      }
    }
  }
  // swap each time step N times
  for(int time_step = 1; (time_step + 1) <  num_steps; time_step++){
    // find current s0 and s1
    s0.clear(); 
    s1.clear();
    s0.reserve(num_agents); 
    s1.reserve(num_agents);
    for(int n = 0; n < num_agents; n++){
      if(xx(n, time_step)){
        s1.push_back(n);
      }else{
        s0.push_back(n);
      }
    }
    infection_t = sum(xx(_,time_step));
    if(infection_t < num_agents && infection_t > 0){
      alphat2next = sis_get_alpha_full_cpp(xx(_, time_step), lambda,gamma);
      alphaprev2t = sis_get_alpha_full_cpp(xx(_, time_step  - 1), lambda, gamma) ;
      logoddsprev2t = log(alphaprev2t) - log(1 - alphaprev2t);
      for(int iswap = 0; iswap < num_swaps; iswap++){
        // choose i1 and i0, indices to swap
        i1_index = floor(unif_rand() * (infection_t));
        i0_index = floor(unif_rand() * (num_agents - infection_t));
        i1 = s1[i1_index];
        i0 = s0[i0_index];
        logpswap = logoddsprev2t[i0] - logoddsprev2t[i1] ;
        logpswap += ((xx(i1, time_step + 1))? log(lambda[i1] * infection_t / num_agents) : log(1 - lambda[i1] * infection_t / num_agents));
        logpswap += ((xx(i0, time_step + 1))? log(1 - gamma[i0]): log(gamma[i0]));
        logpswap -= ((xx(i1, time_step + 1))? log(alphat2next[i1]) : log(1 - alphat2next[i1]));
        logpswap -= ((xx(i0, time_step + 1))? log(alphat2next[i0]) : log(1 - alphat2next[i0]));
        if(log(unif_rand()) < logpswap){
          i1_temp = i1;
          s1[i1_index] = i0;
          s0[i0_index] = i1_temp;
          xx(i1, time_step) = false;
          xx(i0, time_step) = true;
          alphat2next[i0] = (1 - gamma[i0]);
          alphat2next[i1] = lambda[i1] * infection_t / num_agents;
        }
      }
    }
  }
  // terminal step
  time_step = num_steps - 1 ;
  // find current s0 and s1
  s0.clear(); 
  s1.clear();
  s0.reserve(num_agents); 
  s1.reserve(num_agents);
  for(int n = 0; n < num_agents; n++){
    if(xx(n, time_step)){
      s1.push_back(n);
    }else{
      s0.push_back(n);
    }
  }
  infection_t = sum(xx(_,time_step));
  if(infection_t < num_agents && infection_t > 0){
    alphaprev2t = sis_get_alpha_full_cpp(xx(_, time_step  - 1), lambda, gamma) ;
    logoddsprev2t = log(alphaprev2t) - log(1 - alphaprev2t);
    for(int iswap = 0; iswap < num_swaps; iswap++){
      // choose i1 and i0, indices to swap
      i1_index = floor(unif_rand() * (infection_t));
      i0_index = floor(unif_rand() * (num_agents - infection_t));
      i1 = s1[i1_index];
      i0 = s0[i0_index];
      logpswap = logoddsprev2t[i0] - logoddsprev2t[i1] ;
      if(log(unif_rand()) < logpswap){
        i1_temp = i1;
        s1[i1_index] = i0;
        s0[i0_index] = i1_temp;
        xx(i1, time_step) = false;
        xx(i0, time_step) = true;
      }
    }
  }
  PutRNGstate();
}
