#include <Rcpp.h>
#include "boarding_cpp.h"
#include "sir_cpp.h"
#include "logdensities.h"
#include "sampling.h"
using namespace Rcpp;
using namespace std;

/* 
 * Given population size N, there are (N+1) * N / 2 low-dimensional states (s,i) possible. 
 * we create a matrix that stores all these states. 
 * 
 * 0-based index | s | i 
 * --------------|---|----
 *     0         | N | 0 
 *---------------|---|----
 *     1         |N-1| 1 
 *---------------|---|---- 
 */
// [[Rcpp::export]]
IntegerMatrix boarding_all_lowdim_states(const int N){
  int nlowdim_states = (N + 2) * (N + 1) / 2;
  IntegerMatrix allstates(nlowdim_states, 2);
  int istate = 0;
  for(int scnt = N; scnt>=0 ; scnt--){
    for(int icnt = N - scnt; icnt >= 0; icnt--){
      allstates(istate,0) = scnt;
      allstates(istate,1) = icnt;
      // Rcout << "state # " <<  istate <<" s = " << scnt << "i= " << icnt << "\n";
      istate++;
      }
    }
  return allstates;
  }

/*
 * given population size, lowdim summary (s,i)
 * return the index of this state in the matrix generated by boarding_all_lowdim_states
 */
// [[Rcpp::export]]
int boarding_lowdim2index(const int N, const int scnt, const int icnt){
  return int (N - scnt) * (N - scnt + 1) / 2  + N - scnt - icnt; 
  }

/*
 * this function is equivalent to fbar_create + bif_create_fast
 * 
 */
// [[Rcpp::export]]
NumericMatrix boarding_bif_create(const IntegerVector & y,
                                  const IntegerMatrix & all_lowdim_states,
                                  const double & lambda, 
                                  const double & gamma, 
                                  const double & rho,
                                  const int & N,
                                  const double & c){
  int days = y.size();
  // Rcout << "observation length = " <<  days <<"\n";
  
  // initialize the policy vector to -Inf
  NumericMatrix logpolicy(all_lowdim_states.nrow(),  days);
  std::fill(logpolicy.begin(), logpolicy.end(), R_NegInf);
  
  // start the recursion from t = days - 1. 
  int t = days - 1;
  for(int ilowdim = all_lowdim_states.nrow() - 1; ilowdim >= 0; ilowdim--){
    logpolicy(ilowdim, t) =  R::dbinom(y[t], all_lowdim_states(ilowdim, 1), rho, true);
    }
  
  int snow, inow, snext, inext;
  int supper, slower, iupper, ilower;
  double currlsum, ladd, currlmax;
  // run the recursion backwards
  for(int t = days - 2; t >= 0 ; t--){
    for(int ilowdimnow = all_lowdim_states.nrow() - 1; ilowdimnow >= 0; ilowdimnow--){ // for every state
      snow = all_lowdim_states(ilowdimnow, 0);
      inow = all_lowdim_states(ilowdimnow, 1);
      currlsum = R_NegInf;
      // consider a subset of states that can be reaches from (snow, inow);
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
          ladd = sir_logfbar(snow, inow, snext, inext, lambda, gamma, N) + logpolicy(boarding_lowdim2index(N, snext, inext), t+1);
          if(!Rcpp::traits::is_infinite<REALSXP>(ladd)){
            currlmax = max(currlsum, ladd);
            currlsum = log(exp(currlsum - currlmax) + exp(ladd - currlmax)) + currlmax;
            }
          }
        }
        logpolicy(ilowdimnow, t) = currlsum +  R::dbinom(y[t], inow, rho, true);
      }
    }
  return logpolicy;
} 

/*
 * given a policy matrix, update it given new values of y, lambda, and gamma.
 * this function assumes that rho, N, and c stays the same.
 */

// [[Rcpp::export]]
void boarding_bif_update(NumericMatrix logpolicy, 
                         const IntegerVector & y, 
                         const IntegerMatrix & all_lowdim_states,
                         const double & lambda, 
                         const double & gamma, 
                         const double & rho,
                         const int & N,
                         const double &c){
  
  int snow, inow, snext, inext;
  int supper, slower, iupper, ilower;
  double currlsum, ladd, currlmax;
  // run the recursion backwards
  for(int t = y.size()- 2; t >= 0 ; t--){
    for(int ilowdimnow = all_lowdim_states.nrow() - 1; ilowdimnow >= 0; ilowdimnow--){ // for every state
      snow = all_lowdim_states(ilowdimnow, 0);
      inow = all_lowdim_states(ilowdimnow, 1);
      currlsum = R_NegInf;
      // consider a subset of states that can be reaches from (snow, inow);
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
          ladd = sir_logfbar(snow, inow, snext, inext, lambda, gamma, N) + logpolicy(boarding_lowdim2index(N, snext, inext), t+1);
          if(!Rcpp::traits::is_infinite<REALSXP>(ladd)){
            currlmax = max(currlsum, ladd);
            currlsum = log(exp(currlsum - currlmax) + exp(ladd - currlmax)) + currlmax;
          }
        }
      }
      logpolicy(ilowdimnow, t) = currlsum +  R::dbinom(y[t], inow, rho, true);
    }
  }
}

/*
 * given a factor of c, 
 * for each lowdimensional state(st,it),
 * gives a list of states (s(t+1),i(t+1))
 * and compute the transition density from t to t+1
 */
// [[Rcpp::export]]
List boarding_fbar_create(const double lambda, 
                          const double gamma,
                          const int N,
                          const double c){
  // for every eligible state (st,it), compute the vector storing transition kernel
  // log f(s(t+1), i(t+1) | st,it)
  int ncurrstates = (N+1) * (N+2) / 2;
  int nnextstates = c * c * N + 4 + 4 * c * sqrt(N) + 1 + 1; // this is a crude upper bound
  nnextstates = min(nnextstates , ncurrstates); 
  NumericMatrix logfbar(ncurrstates, nnextstates);
  IntegerMatrix nextsi(ncurrstates, nnextstates);
  std::fill(nextsi.begin(), nextsi.end(), -1);
  std::fill(logfbar.begin(), logfbar.end(), 1.0);
  
  int snow, inow, snext, inext;
  int supper, slower, iupper, ilower;
  int istatenow;
  int cntstatenext;
  istatenow = 0;
      
  for(int snow = N; snow >= 0; snow--){
    for(int inow = N - snow; inow >= 0; inow--){ // istatenow  = boarding_lowdim2index(N,snow,inow);
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
          nextsi(istatenow, cntstatenext) = boarding_lowdim2index(N, snext, inext);
          cntstatenext++;
        }
      }
      istatenow++;
    }
  }
  
  return List::create(Named("nextsi") = nextsi,
                      Named("logfbar") = logfbar);
  }

/*
 * if logfbar and nextsi already exsit for factor c,
 * update them given new values of lambda and gamma
 */

// [[Rcpp::export]]
void boarding_fbar_update(NumericMatrix logfbar,
                          IntegerMatrix nextsi,
                          const double lambda,
                          const double gamma,
                          const double N,
                          const double c){
  std::fill(nextsi.begin(), nextsi.end(), -1);
  std::fill(logfbar.begin(), logfbar.end(), 1.0);
  
  int snow, inow, snext, inext;
  int supper, slower, iupper, ilower;
  int istatenow;
  int cntstatenext;
  istatenow = 0;
  
  for(int snow = N; snow >= 0; snow--){
    for(int inow = N - snow; inow >= 0; inow--){ // istatenow  = boarding_lowdim2index(N,snow,inow);
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
          nextsi(istatenow, cntstatenext) = boarding_lowdim2index(N, snext, inext);
          cntstatenext++;
        }
      }
      istatenow++;
    }
  }
}

/* 
 * given the result of pre-computed logfbar matrix, 
 * create the logpolicy matrix using approximate BIF
 */

// [[Rcpp::export]]
NumericMatrix boarding_bif_create_fast(const IntegerVector & y,
                                       const IntegerMatrix & nextsi,
                                       const NumericMatrix & logfbar,
                                       const double & rho,
                                       const int & N,
                                       const double & c){
  int days = y.size();
  // Rcout << "observation length = " <<  days <<"\n";
  int nlowdim_states = (N + 2) * (N + 1) / 2;
  // initialize the policy vector to -Inf
  NumericMatrix logpolicy(nlowdim_states,  days);
  std::fill(logpolicy.begin(), logpolicy.end(), R_NegInf);
  
  // start the recursion from t = days - 1. 
  int t = days - 1;
  int ilowdim = 0;
  for(int snow = N; snow >= 0; snow--){
    for(int inow = N - snow; inow >= 0; inow--){ // istatenow  = boarding_lowdim2index(N,snow,inow);
      logpolicy(ilowdim, t) =  R::dbinom(y[t], inow, rho, true);
      ilowdim++;
    }
  }
  
  int snow, inow;
  double currlsum, ladd, currlmax;
  int icolnext;
  int ilowdimnow = 0;
  
  // run the recursion backwards
  for(int t = days - 2; t >= 0 ; t--){
    ilowdimnow = 0;
    for(int snow = N; snow >= 0; snow--){
      for(int inow = N - snow; inow >= 0; inow--){ // istatenow  = boarding_lowdim2index(N,snow,inow);
        currlsum = R_NegInf;
        icolnext = 0;
        while(icolnext < logfbar.ncol() && nextsi(ilowdimnow, icolnext) >= 0){
          ladd = logfbar(ilowdimnow, icolnext) + logpolicy(nextsi(ilowdimnow, icolnext), t+1);
          if(!Rcpp::traits::is_infinite<REALSXP>(ladd)){
            currlmax = max(currlsum, ladd);
            currlsum = log(exp(currlsum - currlmax) + exp(ladd - currlmax)) + currlmax;
          }
          icolnext++;
        }
        logpolicy(ilowdimnow, t) = currlsum +  R::dbinom(y[t], inow, rho, true);
        ilowdimnow++;
      }
    }
  }
  return logpolicy;
}


/* given updated logfbar matrix,
 * update the logpolicy using approximate BIF
 */
// [[Rcpp::export]]
void boarding_bif_update_fast(NumericMatrix logpolicy,
                              const IntegerVector & y,
                              const NumericMatrix & logfbar,
                              const IntegerMatrix & nextsi,
                              const double & rho,
                              const int & N,
                              const double & c){
  int days = y.size();
  int snow, inow;
  double currlsum, ladd, currlmax;
  int icolnext;
  int ilowdimnow = 0;
  
  for(int t = days - 2; t >= 0 ; t--){
    ilowdimnow = 0;
    for(int snow = N; snow >= 0; snow--){
      for(int inow = N - snow; inow >= 0; inow--){ // istatenow  = boarding_lowdim2index(N,snow,inow);
        currlsum = R_NegInf;
        icolnext = 0;
        while(icolnext < logfbar.ncol() && nextsi(ilowdimnow, icolnext) >= 0){
          ladd = logfbar(ilowdimnow, icolnext) + logpolicy(nextsi(ilowdimnow, icolnext), t+1);
          if(!Rcpp::traits::is_infinite<REALSXP>(ladd)){
            currlmax = max(currlsum, ladd);
            currlsum = log(exp(currlsum - currlmax) + exp(ladd - currlmax)) + currlmax;
          }
          icolnext++;
        }
        logpolicy(ilowdimnow, t) = currlsum +  R::dbinom(y[t], inow, rho, true);
        ilowdimnow++;
      }
    }
  }
}


/*
 * given all particles (xts), update the following items:
 * (1) logf(s(t+1), i(t+1) | xt)
 * (2) alpha_{s to i}(xt)
 * (3) alpha_{i to i}(xt)
 * given the paramters lambda, gamma and the social network in neighbors
 */

// [[Rcpp::export]]

void boarding_logf_update_sparse(NumericMatrix logf,
                          NumericMatrix alphas2i,
                          NumericMatrix alphai2i,
                          const IntegerMatrix  & xts,
                          const NumericVector & lambda,
                          const NumericVector & gamma,
                          const IntegerMatrix & neighbors,
                          const int & N){
  std::fill(alphas2i.begin(), alphas2i.end(), 0);
  std::fill(alphai2i.begin(), alphai2i.end(), 0);
  std::fill(logf.begin(), logf.end(), R_NegInf);
  int ineighbors, mindex;
  int scnt, icnt;
  int snew, inew;
  NumericVector di2i(N + 1);
  NumericVector ds2i(N + 1);
  
  for (int p = 0; p < xts.ncol(); p++){
    // prepare for next particle
    scnt = 0;
    icnt = 0;
    for(int n = 0; n < xts.nrow(); n++){
      if(xts(n,p) == 0){
        mindex = ineighbors = 0;
        while(mindex < neighbors.ncol() & neighbors(n,mindex) >= 0){
          ineighbors += ((xts(neighbors(n,mindex),p) == 1)? 1 : 0);
          mindex++;
        }
        // Rprintf("agent %i has %i neighbors\n", n, mindex);
        // mindex is the number of neighbors of agent n 
        alphas2i(n,p) = lambda[n] * ineighbors / mindex;
        scnt++;
      }else if(xts(n,p) == 1){
        alphai2i(n,p) = 1 - gamma[n];
        icnt++;
      }
    }
    // alpha updates finished
    di2i = logdpoisbinom_cpp(alphai2i(_,p));
    ds2i = logdpoisbinom_cpp(alphas2i(_,p));
    // Rprintf("probability of no change is %.2f \n", ds2i[0] + di2i[icnt]);
    // Rprintf("(s,i) = (%i, %i)\n", scnt, icnt);
    for(int s2i = 0; s2i <= scnt; s2i++){
      for(int i2i = 0; i2i <= icnt; i2i++){
        inew = s2i + i2i;
        snew = scnt - s2i;
        logf(boarding_lowdim2index(N, snew, inew), p) = di2i[i2i] + ds2i[s2i];
      }
    }
  }
}

/*
 * given all particles (xts) and assuming a fully connected social network structure,
 * update the following items:
 * (1) logf(s(t+1), i(t+1) | xt)
 * (2) alpha_{s to i}(xt)
 * (3) alpha_{i to i}(xt)
 * given the paramters lambda, gamma.
 */

// [[Rcpp::export]]

void boarding_logf_update_full(NumericMatrix logf,
                                 NumericMatrix alphas2i,
                                 NumericMatrix alphai2i,
                                 const IntegerMatrix  & xts,
                                 const NumericVector & lambda,
                                 const NumericVector & gamma,
                                 const int & N){
  std::fill(alphas2i.begin(), alphas2i.end(), 0);
  std::fill(alphai2i.begin(), alphai2i.end(), 0);
  std::fill(logf.begin(), logf.end(), R_NegInf);
  int ineighbors;
  int scnt, icnt;
  int snew, inew;
  NumericVector di2i(N + 1);
  NumericVector ds2i(N + 1);
  
  for (int p = 0; p < xts.ncol(); p++){
    // prepare for next particle
    scnt = 0;
    icnt = 0;
    ineighbors = sum(xts(_,p) == 1); // the number of infections in the pth particle
    for(int n = 0; n < xts.nrow(); n++){
      if(xts(n,p) == 0){
        alphas2i(n,p) = lambda[n] * ineighbors / N;
        scnt++;
      }else if(xts(n,p) == 1){
        alphai2i(n,p) = 1 - gamma[n];
        icnt++;
      }
    }
    // alpha updates finished
    di2i = logdpoisbinom_cpp(alphai2i(_,p));
    ds2i = logdpoisbinom_cpp(alphas2i(_,p));
    // Rprintf("probability of no change is %.2f \n", ds2i[0] + di2i[icnt]);
    // Rprintf("(s,i) = (%i, %i)\n", scnt, icnt);
    for(int s2i = 0; s2i <= scnt; s2i++){
      for(int i2i = 0; i2i <= icnt; i2i++){
        inew = s2i + i2i;
        snew = scnt - s2i;
        logf(boarding_lowdim2index(N, snew, inew), p) = di2i[i2i] + ds2i[s2i];
      }
    }
  }
}


/* 
 * given the particles xts, and the updates matrices alphas2i and alphai2i, 
 * sample from the conditional kernel f(x(t+1) | xt, s(t+1), i(t+1)) for each particle
 */

// [[Rcpp::export]]
void boarding_sample_x_given_si(IntegerMatrix & xts,
                                       const NumericMatrix & alphas2i,
                                       const NumericMatrix & alphai2i,
                                       const NumericVector & lambda,
                                       const NumericVector & gamma,
                                       const IntegerVector & snext,
                                       const IntegerVector & inext,
                                       const int & N,
                                       const int & P){
  int snow, inow;
  int s2i, i2i;
  LogicalVector xxsi(N);
  LogicalVector xxii(N);
  LogicalVector xxi(N);
  NumericVector rand(N);
  
  
  for(int p = 0; p < P; p++){
    snow = sum(xts(_,p) == 0);
    inow = sum(xts(_,p) == 1);
    s2i = snow - snext[p];
    i2i = inext[p] - s2i;
    
    GetRNGstate();
    rand = runif(N);
    PutRNGstate();
    xxsi = idchecking_cpp(s2i, alphas2i(_,p), rand);
    GetRNGstate();
    rand = runif(N);
    PutRNGstate();
    
    xxii = idchecking_cpp(i2i, alphai2i(_,p), rand);
    xxi = mapply(xxsi, xxii, xory);
    
    for(int n = 0; n < N; n++){
      if(xxi[n]){
        xts(n,p) = 1;
      }else{//xxi[n] == 0
        // 0 -> 0; 1 -> 2; 2 -> 2;
        xts(n,p) = (xts(n,p) == 0)? 0 : 2;
      }
    }
  }
  // return xts;
}
