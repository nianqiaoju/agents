#include <Rcpp.h>
#include "sis_cpp.h"
using namespace Rcpp;

// performs blocked on the hidden states updates given parameters and observations
// it modifies the input xx directly
// in every block, all the agents states at all times are updated jointly
// using a forward-backward algorithm
// current implementation assumes a fully connected network
// this function assumes that population size is divisible by block_size 
// it also assumes that the input state_space is a matrix of size 2^(block_size) by population_size


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