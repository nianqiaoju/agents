#include <Rcpp.h>
#include "sir_cpp.h"
#include "logdensities.h"
#include "sampling.h"
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
double sir_logfbar(const int & snow,
                   const int & inow, 
                   const int & snext, 
                   const int & inext,
                   const double & lambda, 
                   const double & gamma, 
                   const int & N){
  return(R::dbinom(snow + inow - inext - snext, inow, gamma, true) + R::dbinom(snow - snext, snow, lambda * inow / N, true));
};


bool xory(bool x, bool y){
  return(x || y);
};

int sir_get_state_index(const int s, const int i, const int N){
  int index = (N - s) * (N- s + 1) / 2 + N - s - i;
  return(index);
}



/* perform the BIF approximation for the SIR model, that observes icount partially
 * after the approximate trasition kernel fbar has been created by create_fbar_matrix 
 * or modified in update_fbar_matrix
 */

// [[Rcpp::export]]
NumericVector sir_bif_create_cpp(const IntegerVector & y, 
                                 const IntegerMatrix & nexts,
                                 const IntegerMatrix & nexti,
                                 const NumericMatrix & logfbar,
                                 const double & rho, 
                                 const int & N) {
  int days = y.size();
  // Rcout << "observation length = " <<  days <<"\n";
  // initialize a vector to -Inf
  NumericVector logpolicy((N+1) * (N+1)  * days);
  std::fill(logpolicy.begin(), logpolicy.end(), R_NegInf);
  // deal with the terminal time first
  int t = days - 1;
  // Rprintf("at day t = %i\n", t);
  for(int inow = N; inow >= y[t]; inow--){
    for(int snow = 0; snow + inow <= N; snow++){
      // Rprintf("s = %i, i = %i\n", snow, inow);
      logpolicy[t * (N+1) * (N+1) + (N+1) * inow + snow] = R::dbinom(y[t], inow, rho, true);
    }
  }
  NumericVector logg(N + 1);
  double curr_sum, add_item, curr_max;
  int cntstatenext;
  int istatenow;
  for(int t = days - 2; t >= 0; t--){
    // Rprintf("at day t = %i\n", t);
    // update observation densities at time t
    std::fill(logg.begin(), logg.end(), R_NegInf);
    for(int i = y[t]; i <= N; i++){
      logg[i] = R::dbinom(y[t], i, rho, true);
    }
    
    // compute the conditional expectation for each value of st and it
    for(int snow = N - y[t]; snow >= 0; snow--){
      istatenow = sir_get_state_index(snow, N - snow, N); // first eligible state for snow;
      for(int inow = N - snow; inow >= y[t]; inow--){
        curr_sum = R_NegInf; // initialize the conditional expectation to 0
        // can also set curr_sum to a very small value to satisfy the sufficient support condition
        cntstatenext = 0;
        while(cntstatenext < logfbar.ncol() && nexti(istatenow, cntstatenext) >= 0){
          add_item = logfbar(istatenow, cntstatenext) + logpolicy[(t + 1) * (N+1) * (N+1) + (N+1) * nexti(istatenow, cntstatenext) + nexts(istatenow, cntstatenext)];
          if(!Rcpp::traits::is_infinite<REALSXP>(add_item)){
            curr_max = max(curr_sum, add_item);
            curr_sum = log(exp(curr_sum - curr_max) + exp(add_item - curr_max)) + curr_max;
          }
          cntstatenext++;
        }
        // psi_t(st,it) is observation density * conditional expectation
        logpolicy[t * (N+1) * (N+1) + (N+1) * (inow) + snow] = curr_sum + logg[inow];
        istatenow++;
      }
    }
  }
  // convert vector to  3d-array
  IntegerVector dim = {N+1, N+1, days};
  logpolicy.attr("dim") = dim;
  return logpolicy;
}

/* for the sir dynamics with icount observed partially,
 * update the backward information filter given new values of lambda, gamma and rho.
 * This function should be run after sir_bif_create_cpp and create_fbar_matrix, 
 * and it modifies logpolicy and logfbar matrix directly.
 */

// [[Rcpp::export]]
void sir_bif_update_cpp(NumericVector & logpolicy, 
                        IntegerMatrix & nexts,
                        IntegerMatrix & nexti, 
                        NumericMatrix & logfbar, 
                        const double & lambda,
                        const double & gamma, 
                        const double & rho,
                        const IntegerVector & y, 
                        const int & N,
                        const double & c){
  update_fbar_matrix(nexts, nexti, logfbar, lambda, gamma, N, c);
  std::fill(logpolicy.begin(), logpolicy.end(), R_NegInf);
  // deal with the terminal time first
  int days = y.size();
  int t = days - 1;
  // Rprintf("at day t = %i\n", t);
  for(int inow = N; inow >= y[t]; inow--){
    for(int snow = 0; snow + inow <= N; snow++){
      // Rprintf("s = %i, i = %i\n", snow, inow);
      logpolicy[t * (N+1) * (N+1) + (N+1) * inow + snow] = R::dbinom(y[t], inow, rho, true);
    }
  }
  NumericVector logg(N + 1);
  double curr_sum, add_item, curr_max;
  int cntstatenext;
  int istatenow;
  for(int t = days - 2; t >= 0; t--){
    // Rprintf("at day t = %i\n", t);
    // update observation densities at time t
    std::fill(logg.begin(), logg.end(), R_NegInf);
    for(int i = y[t]; i <= N; i++){
      logg[i] = R::dbinom(y[t], i, rho, true);
    }
    
    // compute the conditional expecation for each value of st and it
    for(int snow = N - y[t]; snow >= 0; snow--){
      istatenow = sir_get_state_index(snow, N - snow, N); 
      for(int inow = N - snow; inow >= y[t]; inow--){
        curr_sum = R_NegInf; // initialize the conditonal expectation to 0
        // can also set curr_sum to a very small value to satisfy the sufficient support condition
        cntstatenext = 0;
        while(cntstatenext < logfbar.ncol() && nexti(istatenow, cntstatenext) >= 0){
          add_item = logfbar(istatenow, cntstatenext) + logpolicy[(t + 1) * (N+1) * (N+1) + (N+1) * nexti(istatenow, cntstatenext) + nexts(istatenow, cntstatenext)];
          if(!Rcpp::traits::is_infinite<REALSXP>(add_item)){
            curr_max = max(curr_sum, add_item);
            curr_sum = log(exp(curr_sum - curr_max) + exp(add_item - curr_max)) + curr_max;
          }
          cntstatenext++;
        }
        // psi_t(st,it) is observation density * conditional expectation
        logpolicy[t * (N+1) * (N+1) + (N+1) * (inow) + snow] = curr_sum + logg[inow];
        istatenow++;
      }
    }
  }
}

/*
 * create (subset of) the transition matrix and associated index vectors for the SIR dynamics
 * fbar matrix is the transition between aggregated states (st,it)
 * the factor c controls size of the subset
 */
// [[Rcpp::export]]
List create_fbar_matrix(const double lambda, 
                        const double gamma,
                        const int N,
                        const double c){
  // for every eligible state (st,it), compute the vector storing transition kernel
  // log f(s(t+1), i(t+1) | st,it)
  int ncurrstates = (N+1) * (N+2) / 2;
  int nnextstates = c * c * N + 4 + 4 * c * sqrt(N) + 1 + 1; // this is a crude upper bound
  nnextstates = min(nnextstates , ncurrstates); 
  // Rprintf("matrix dimension is %i by %i\n", ncurrstates, nnextstates);
  // row index correspond to current state
  // column index correspond to next state
  NumericMatrix logfbar(ncurrstates, nnextstates);
  IntegerMatrix nexts(ncurrstates, nnextstates);
  IntegerMatrix nexti(ncurrstates, nnextstates);
  std::fill(nexts.begin(), nexts.end(), -1);
  std::fill(nexti.begin(), nexti.end(), -1);
  std::fill(logfbar.begin(), logfbar.end(), 1.0);
  
  int snow, inow, snext, inext;
  int supper, slower, iupper, ilower;
  int istatenow;
  int cntstatenext;
  istatenow = 0;
  
  for(int snow = N; snow >= 0; snow--){
    for(int inow = N - snow; inow >= 0; inow--){
      // snow and inow defines a current state
      cntstatenext = 0;
      slower = snow * (1 - inow * lambda / N) - c * sqrt(N) / 2;
      supper = snow * (1 - inow * lambda / N) + c * sqrt(N) / 2 + 0.5;
      slower = max(0, slower);
      supper = min(snow, supper);
      for(int snext = supper; snext >= slower; snext--){
        ilower = inow + snow - snext - inow * gamma - c * sqrt(N) / 2;
        iupper = inow + snow - snext - inow * gamma + c * sqrt(N) / 2 + 0.5;
        ilower = max(snow - snext, ilower);
        iupper = min(inow + snow - snext, iupper);
        for(int inext = iupper; inext >= ilower; inext--){
          // snext and inext defines a future state
          logfbar(istatenow, cntstatenext) = sir_logfbar(snow, inow, snext, inext, lambda, gamma, N);
          nexts(istatenow, cntstatenext) = snext;
          nexti(istatenow, cntstatenext) = inext;
          cntstatenext++;
        }
      }
      istatenow++;
      // Rprintf("number of states = %i and analytic nnextstates = %i\n", cntstatenext, nnextstates);
    }
  }
  return List::create(Named("nexts") = nexts,
                      Named("nexti") = nexti,
                      // Named("nextstates") = next_states,
                      Named("logfbar") = logfbar);
}



/* assume that logfbar, nexts, nexti already exist
 * update these matrices given new values of lambda, gamma
 * this function assumes that c is the same as the factor used in create_fbar_matrix
 */

// [[Rcpp::export]]
void update_fbar_matrix(IntegerMatrix & nexts,
                        IntegerMatrix & nexti, 
                        NumericMatrix & logfbar, 
                        const double & lambda,
                        const double & gamma, 
                        const int & N,
                        const double & c){
  // row index correspond to current state
  // column index correspond to next state
  std::fill(nexts.begin(), nexts.end(), -1);
  std::fill(nexti.begin(), nexti.end(), -1);
  std::fill(logfbar.begin(), logfbar.end(), 1.0);
  
  int snow, inow, snext, inext;
  int supper, slower, iupper, ilower;
  int istatenow;
  int cntstatenext;
  istatenow = 0;
  
  for(int snow = N; snow >= 0; snow--){
    for(int inow = N - snow; inow >= 0; inow--){
      // snow and inow defines a current state
      // get lower and upper bounds for future states
      slower = snow * (1 - inow * lambda / N) - c * sqrt(N) / 2;
      supper = snow * (1 - inow * lambda / N) + c * sqrt(N) / 2 + 0.5;
      slower = max(0, slower);
      supper = min(snow, supper);
      cntstatenext = 0;
      for(int snext = supper; snext >= slower; snext--){
        ilower = inow + snow - snext - inow * gamma - c * sqrt(N) / 2;
        iupper = inow + snow - snext - inow * gamma + c * sqrt(N) / 2 + 0.5;
        ilower = max(snow - snext, ilower);
        iupper = min(inow + snow - snext, iupper);
        for(int inext = iupper; inext >= ilower; inext--){
          // snext and inext defines a future state
          logfbar(istatenow, cntstatenext) = sir_logfbar(snow, inow, snext, inext, lambda, gamma, N);
          nexts(istatenow, cntstatenext) = snext;
          nexti(istatenow, cntstatenext) = inext;
          cntstatenext++;
        }
      }
      istatenow++;
    }
  }
}

/* this function will run within the controlled particle filter
 * this is the cpp implementation of sir_logdpoismulti
 * compute the transition probability from x[t-1] to the aggregated state (st, it)
 * log f(st, it | x[t-1], theta)
 * it directly modifies the logf matrix
 * the current implementation assumes that the network is fully connected
 */

// [[Rcpp::export]]
void sir_csmc_update_f_matrix(NumericMatrix & logf,
                         const IntegerVector & xxprev,
                         const NumericVector & lambda_v,
                         const NumericVector & gamma_v,
                         const int & N){
  // assume that logf is a (N+1) by (N+1) matrix
  // xxprev is a length N vector representing agent states
  std::fill(logf.begin(), logf.end(), R_NegInf);
  int icount, scount;
  int inew, snew;
  NumericVector alphai2i(N);
  NumericVector alphas2i(N);
  NumericVector di2i(N);
  NumericVector ds2i(N);
  
  icount = sum(xxprev == 1);
  scount = sum(xxprev == 0);
  // Rprintf("iprev = %i, sprev = %i\n", icount, scount);
  for(int id = 0; id < N; id++){
    if(xxprev[id] == 0){
      alphas2i[id] = lambda_v[id] * icount / N;
    }else if(xxprev[id] == 1){
      alphai2i[id] = 1 - gamma_v[id];
    }
  }
  
  di2i = logdpoisbinom_cpp(alphai2i);
  ds2i = logdpoisbinom_cpp(alphas2i);
  // Rprintf("probability of no change is %.2f \n", ds2i[0] + di2i[icount]);
  
  for(int s2i = 0; s2i <= scount; s2i++){
    for(int i2i = 0; i2i <= icount; i2i++){
      inew = s2i + i2i;
      snew = scount - s2i;
      logf(snew, inew) = di2i[i2i] + ds2i[s2i];
      // Rprintf("st = %i, it = %i, logf = %.2f\n", snew, inew, logf(snew, inew));
    }
  }
}




/* For the sir model,
 * sample from the kernel f( * | xt) subject to 
 * the constraint that S(x[t+1]) = s and I(x[t+1]) = i
 * It modifies the matrixc xx directly, which is N by P
 */ 

// [[Rcpp::export]]
IntegerMatrix sir_sample_x_given_si(IntegerMatrix & xx,
                           const NumericVector & lambda,
                           const NumericVector & gamma, 
                           const IntegerVector & scount,
                           const IntegerVector & icount,
                           const int & N,
                           const int & P){
  NumericVector alphasi(N);
  NumericVector alphaii(N);
  LogicalVector xxsi(N);
  LogicalVector xxii(N);
  LogicalVector xxi(N);
  int snow, inow, rnow, rcount;
  int s2i, i2i;
  NumericVector rand(N);
  
  for(int p = 0; p < P; p++){// for every particle
    std::fill(alphasi.begin(), alphasi.end(), 0);
    std::fill(alphaii.begin(), alphaii.end(), 0);
    snow = sum(xx(_,p) == 0);
    inow = sum(xx(_,p) == 1);
    rnow = N - snow - inow;
    rcount = N - scount[p] - icount[p];
    s2i = snow - scount[p];
    i2i = icount[p] - s2i;
    // update alphasi and alphaii
    for(int n = 0; n < N; n++){
      if(xx(n,p) == 0){
        alphasi[n] = lambda[n] * inow / N;
      }else if (xx(n,p) == 1){
        alphaii[n] = 1 - gamma[n];
      }
    }
    GetRNGstate();
    rand = runif(N);
    PutRNGstate();
    xxsi = idchecking_cpp(s2i, alphasi, rand);
    GetRNGstate();
    rand = runif(N);
    PutRNGstate();
    xxii = idchecking_cpp(i2i, alphaii, rand);
    xxi = mapply(xxsi, xxii, xory);
    // compare xxi and xx(_,p), and change the values in xx(_,p)
    for(int n = 0; n < N; n++){
      if(xxi[n]){
        xx(n,p) = 1;
      }else{//xxi[n] == 0
        // 0 -> 0; 1 -> 2; 2 -> 2;
        xx(n,p) = (xx(n,p) == 0)? 0 : 2;
      }
    }
  }
  return xx;
};


