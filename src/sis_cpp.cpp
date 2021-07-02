#include <Rcpp.h>
#include "sis_cpp.h"
#include "logdensities.h"
#include "logsumexp.h"
using namespace Rcpp;
using namespace std;



/* 
 * compute the infection probabilities given previous state, lambda and gamma
 * for the fully connected network
 */

// [[Rcpp::export]]
NumericVector sis_get_alpha_full_cpp(const LogicalVector &xx,
                                     const NumericVector &lambda, 
                                     const NumericVector &gamma){
  NumericVector alpha (xx.length());
  double num_infections = sum(xx);
  // Rcout << "there are " << num_infections << "infections\n";
  double prop_infections = num_infections / xx.length();
  for(int n = 0; n < xx.length(); n++){
    alpha[n] = (xx[n]) ? (1 - gamma[n]) : (prop_infections * lambda[n]);
  }
  return alpha;
}

/*
 * for a sparse network, 
 * compute the infection probabilities given previous state, lambda and gamma.
 * network given in matrix neighbors, where neighbors(n, _) is a list of neighbors of agent n.
 * neighbors(n, _) = (1,2,4,-1,-1,-1) means that agent 1,2, and 4 are neighbors with agent n. 
 */
// [[Rcpp::export]]
NumericVector sis_get_alpha_sparse(const LogicalVector & xx,
                                   const NumericVector & lambda,
                                   const NumericVector & gamma,
                                   const IntegerMatrix & neighbors){
  NumericVector alpha(xx.length());
  int sumIn; // sum of infectious neighbors
  int Dn; // number of neighbors
  for(int n = 0; n< xx.length(); n++){
    sumIn = 0;
    Dn = 0;
    if(xx[n]){
      alpha[n] = 1 - gamma[n];
    }else{
      // alpha[n] = lambda[n] * number of infected neighors of n / number of neighbors of n 
      while(neighbors(n, Dn) > 0 & Dn < neighbors.ncol()){
        sumIn += xx[neighbors(n, Dn)];
        Dn ++;
        }
      if(sumIn > 0){alpha[n] = lambda[n] * sumIn / Dn;}
      }
    }
  return alpha;
  }

/*
 * given parameters and observations,performs blocked updates of the hidden states X;
 * it modifies the input xx directly;
 * in every block, all the agents states at all times are updated jointly, using a forward-backward algorithm;
 *  assumptions:
 *  1. fully connected network
 *  2. population size N is divisible by block_size
 *  3. state_space is 2**(block_size) by population size;
*/

// [[Rcpp::export]]
void sis_xx_gibbs_blocked_full_cpp(LogicalMatrix  xx,
                                   const IntegerVector &y,
                                   const NumericVector &alpha0,
                                   const NumericVector &lambda,
                                   const NumericVector &gamma,
                                   const double rho,
                                   const int block_size,
                                   const LogicalMatrix &state_space){
  int population_size = xx.nrow();
  int nblocks = population_size / block_size;
  int nstates = state_space.nrow();
  int nobs = y.size();
  
  IntegerVector sum_state_space = rowSums(state_space);  // I(x) for each block
  NumericMatrix m(nstates, nobs);
  // m(row , col) = log p(x[row] & y[0 : col]), instatiate it to -Inf
  for(int istate = 0; istate < nstates; istate++){
    for(int t = 0; t < nobs; t++){
      m(istate,t) = R_NegInf;
    }
  }
  
  NumericVector logweights(nstates);
  // weights to sample a state (using multinomial_cpp)
  
  LogicalVector xxnow(population_size); 
  LogicalVector xxprev(population_size);
  
  double f; // transition probabilities
  double maxl; // used for log sum exp operations
  double logfm;
  
  // record the agents in the block being updated
  // the first block is [01...B]
  // before the forward-backward procedure starts,
  // set this block to 0s
  IntegerVector block_members (block_size);
  for(int imember = 0; imember < block_size; imember++){
    block_members[imember] = imember;
    std::fill(xx.row(block_members[imember]).begin(), xx.row(block_members[imember]).end(), false);
  }
  // number of infections at each time outside of the block
  IntegerVector infection_outside (nobs);
  infection_outside = colSums(xx);
  
  // update each block
  for(int iblock = 0; iblock < nblocks ; iblock++){
    // forward algorithm, marginalizes over the hidden states
    // t = 0 uses alpha0
    for(int jstate = 0; jstate < nstates; jstate++){ // current state
      m(jstate, 0) = R::dbinom(y[0], sum_state_space[jstate] + infection_outside[0], rho, true);
      for(int imember = 0; imember < block_size; imember++){
        m(jstate, 0) += (state_space(jstate, imember)? log(alpha0[block_members[imember]]) : log(1- alpha0[block_members[imember]]));
      }
    }
    // t > 0, the transition depends on the agents not in the current block
    for(int t = 1; t < nobs; t++){
      // complete vectors xxprev and xxnow;
      for(int n = 0; n < population_size; n++){
        xxprev[n] = xx(n, t - 1);
        xxnow[n] = xx(n, t);
      }
      // compute the transition probability p(xxnow | xxprev) for each xxprev
      // and marginalize over xxprev
      for(int istate = 0; istate < nstates; istate++){ // previous state
        for(int imember = 0; imember < block_size; imember++){
          xxprev[block_members[imember]] = state_space(istate,imember);
        }
        NumericVector alpha_prev2now = sis_get_alpha_full_cpp(xxprev, lambda, gamma);
        for(int jstate = 0; jstate < nstates; jstate++){ // current state
          for(int jmember = 0; jmember < block_size; jmember++){
            xxnow[block_members[jmember]] = state_space(jstate,jmember);
          }
          f = logdbern_sum_cpp(alpha_prev2now, xxnow);
          // this is the probability of going from istate to jstate at time t
          // we re-compute f in the backward sampling step
          // we can also choose to store all the f's is a nstate ** 2 by nobs matrix
          logfm = f + m(istate, t - 1);
          maxl = max(m(jstate, t), logfm);
          if (traits::is_infinite<REALSXP>(maxl)){
            m(jstate, t) = R_NegInf;
          }
          else{
            m(jstate, t) = maxl + log(exp(m(jstate, t) - maxl) + exp(logfm - maxl));
          }
        }
      }
      for(int jstate = 0; jstate < nstates; jstate++){
        // then include current observation densities
        m(jstate,t) = m(jstate,t) + R::dbinom(y[t], sum_state_space[jstate] + infection_outside[t], rho, true);
      }
    }
    
    // print the log marginal likelihood for this block, for testing. 
    // Rcout << "log marginal likelihood = " << lw_logsum_cpp(m(_, nobs - 1)) <<"\n";
    
    // backward sampling
    // get random numbers
    NumericVector runifs = get_runifs(nobs);
    int t = nobs - 1;
    // sample an index
    int jstate = multinomial_cpp(m(_, t), runifs[t]);
    // assign the binary values to agents in xx
    for(int imember = 0; imember < block_size; imember++){
      xx(block_members[imember],t) = state_space(jstate,imember);
    }
    for(int t = (nobs - 2); t >= 0; t--){
      // sample x(t) given x(t+1)
      // weight(istate) propto m(istate,t) * f(x(t+1)| istate)
      for(int n = 0; n < population_size; n++){
        xxprev[n] = xx(n, t);
        xxnow[n] = xx(n, t + 1);
      }
      for(int istate = 0; istate < nstates; istate++){
        for(int imember = 0; imember < block_size; imember++){
          xxprev[block_members[imember]] = state_space(istate,imember);
        }
        logweights[istate] = logdbern_sum_cpp(sis_get_alpha_full_cpp(xxprev,lambda, gamma), xxnow) + m(istate,t);
      }
      // sample an index
      int istate = multinomial_cpp(logweights, runifs[t]);
      // Rcout << "block " << iblock << " at time t = " << t << " choose state " << istate << "\n";
      for(int imember = 0; imember < block_size; imember++){
        xx(block_members[imember],t) = state_space(istate,imember);
      }
    }
    // prepare next block
    if(iblock != (nblocks - 1)){
      for(int imember = 0; imember < block_size; imember++){
        block_members[imember] = block_members[imember] + block_size;
        std::fill(xx.row(block_members[imember]).begin(), xx.row(block_members[imember]).end(), false);
      }
      infection_outside = colSums(xx);
    }
  }
}


/*
 * a complete single-site scan for all t and all n 
 */

// [[Rcpp::export]]
void sis_xx_gibbs_singlesite_full_cpp(LogicalMatrix xx, 
                                      const IntegerVector &y,
                                      const NumericVector &alpha0,
                                      const NumericVector &lambda, 
                                      const NumericVector &gamma, 
                                      const double rho){
  GetRNGstate();
  
  int population_size = xx.ncol();
  LogicalVector xxt(population_size);
  NumericVector alphaprev2t(population_size);
  NumericVector alphat2next(population_size);
  LogicalVector xxnext(population_size);
  double lw0, lw1, p1, maxlw;
  
  int time_step = 0;
  // Rcout << "t = " << time_step << "\n";
  int infection_t = 0;
  alphaprev2t = alpha0;
  xxt = xx(_ ,time_step);
  xxnext = xx(_, time_step + 1);
  for(int n = 0; n < population_size; n++){
    xxt[n] = false;
    infection_t = sum(xxt);
    alphat2next = sis_get_alpha_full_cpp(xxt,lambda,gamma);
    lw0 = log(1- alphaprev2t[n]) + R::dbinom(y[time_step], infection_t, rho, true) + logdbern_sum_cpp(alphat2next,xxnext);
    xxt[n] = true;
    alphat2next = sis_get_alpha_full_cpp(xxt,lambda,gamma);
    lw1 = log(alphaprev2t[n]) + R::dbinom(y[time_step], infection_t + 1, rho, true) +  logdbern_sum_cpp(alphat2next, xxnext);
    // normalize probabilities
    maxlw = max(lw0,lw1);
    p1 = exp(lw1- maxlw) / (exp(lw1 - maxlw) + exp(lw0 - maxlw));
    // sample xt[n]
    if(unif_rand() < p1){
      xxt[n] = true;
    }else{
      xxt[n] = false;
    }
  }
  xx(_, time_step) = xxt;
  for(time_step = 1; (time_step + 1) < y.length(); time_step++){
    xxt = xx(_, time_step);
    alphaprev2t = sis_get_alpha_full_cpp(xx(_, time_step - 1), lambda, gamma);
    xxnext = xx(_, time_step + 1);
    // Rcout << "BEFORE: t = " << time_step << "y = " <<  y[time_step] << "infection count is " <<  sum(xxt) << "\n";
    for(int n = 0; n < population_size; n++){
      xxt[n] = false;
      infection_t = sum(xxt);
      alphat2next = sis_get_alpha_full_cpp(xxt,lambda,gamma);
      lw0 = log(1 - alphaprev2t[n]) + R::dbinom(y[time_step], infection_t, rho, true) + logdbern_sum_cpp(alphat2next, xxnext);
      xxt[n] = true;
      alphat2next = sis_get_alpha_full_cpp(xxt,lambda,gamma);
      lw1 = log(alphaprev2t[n]) + R::dbinom(y[time_step], infection_t + 1, rho, true) + logdbern_sum_cpp(alphat2next, xxnext);
      // normalize probabilities
      maxlw = max(lw0,lw1);
      p1 = exp(lw1- maxlw) / (exp(lw1 - maxlw) + exp(lw0 - maxlw));
      // sample xt[n]
      if(unif_rand() < p1){
        xxt[n] = true;
      }else{
        xxt[n] = false;
      }
    }
    // Rcout << "AFTER: t = " << time_step << "y = " <<  y[time_step] << "infection count is " <<  sum(xxt) << "\n";
    xx(_, time_step) = xxt;
  }
  // final step
  // Rcout << "t = " << time_step << "\n";
  xxt = xx(_, time_step);
  alphaprev2t = sis_get_alpha_full_cpp(xx(_, time_step - 1), lambda, gamma);
  for(int n = 0; n < population_size; n++){
    xxt[n] = false;
    infection_t = sum(xxt);
    lw0 = log(1 - alphaprev2t[n]) + R::dbinom(y[time_step], infection_t, rho, true);
    xxt[n] = true;
    lw1 = log(alphaprev2t[n]) + R::dbinom(y[time_step], infection_t + 1, rho, true);
    // normalize probabilities
    maxlw = max(lw0,lw1);
    p1 = exp(lw1- maxlw) / (exp(lw1 - maxlw) + exp(lw0 - maxlw));
    // sample xt[n]
    if(unif_rand() < p1){
      xxt[n] = true;
    }else{
      xxt[n] = false;
    }
  }
  xx(_, time_step) = xxt;
  PutRNGstate();
}

/*
 * performs N swap updates for each t
 * current implementation assumes a fully connected network and modifies input xx directly
 */

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


// [[Rcpp::export]]
NumericMatrix sis_forward_algorithm_cpp(NumericVector logf0, 
                                        NumericMatrix logdtransition,  
                                        NumericVector y, 
                                        NumericVector all_sum_x, 
                                        double rho) {
  int size_state_space = logdtransition.nrow(); // this number is in fact 2**N for SIS or 3**N for SIR
  int num_observations = y.size();
  NumericMatrix m(size_state_space, num_observations); // each row is an individual and each column is a day 
  for (int istate = 0; istate < size_state_space; istate ++){
    m(istate,0) = logf0(istate) + R::dbinom(y(0), all_sum_x(istate), rho, true); 
  }
  for (int t = 1; t <= num_observations - 1; t++){
    for (int istate = 0; istate <  size_state_space ; istate++){
      NumericVector logweights = m(_, t - 1) + logdtransition( _ , istate);
      if(all(is_infinite(logweights))){
        m(istate, t) = R_NegInf;
      }else{
        double maxlogweights = max(logweights);
        m(istate,t) = maxlogweights + log(sum(exp(logweights - maxlogweights))) + R::dbinom(y(t), all_sum_x(istate), rho, true);
      }
    }
  }
  return m;
}



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

