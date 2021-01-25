#include <Rcpp.h>
#include "logsumexp.h"
using namespace Rcpp;
using namespace std;

// this function maps the vector (i,r) 
// to its index in the vector representation of state space

// [[Rcpp::export]]
int misc_ir2index_cpp(int i, int r, int N) {
  int index = (2 * N + 3 - i) * i / 2 + r;
  return index;
}

// this function computes the approximation transition kernel fbar
// [[Rcpp::export]]
double logfbar_cpp(double lambda_bar, double gamma_bar, int N, int i0, int r0, int i1, int r1){
	// check all the constraints first
	if(i0 + r0 > N) return(R_NegInf);
	if(i1 + r0 > N) return(R_NegInf);
	if(r1 < r0) return(R_NegInf);
	if(N - i0 - r0 < N - i1 - r1) return(R_NegInf);
	if(r1 - r0 > i0) return(R_NegInf);
	if(N - i0 - r0 < i1 - i0 + r1 - r0) return(R_NegInf);// [[Rcpp::export]]
	double infection_prob =  std::min<double>(1,lambda_bar * i0 / N);
	double x = R::dbinom(r1 - r0, i0, gamma_bar,TRUE) + R::dbinom(i1 - i0 + r1 - r0, N - i0 - r0, infection_prob , TRUE);
	// double x = R::dbinom(r1 - r0, i0, gamma_bar,TRUE) + R::dbinom(i1 - i0 + r1 - r0, N - i0 - r0, lambda_bar * i0 / N , TRUE);
	return(x);
}

// [[Rcpp::export]]
NumericMatrix sir_bif_policy_matrix_cpp(double lambda_bar, double gamma_bar,
	IntegerVector y, int N, double rho){
	int num_observations = y.length();
	int irspace_size = (N + 1) * (N + 2) / 2 ;
	// instantize the policy matrix at -Inf
	NumericMatrix policy_matrix(irspace_size,num_observations);
	std::fill(policy_matrix.begin(), policy_matrix.end(), R_NegInf);
	// // use vectors to track index2i and index2r;
	// IntegerVector index2i(irspace_size);
	// IntegerVector index2r(irspace_size);
	// for(int i = 0; i <=N; i++){
	// 	for(int r = 0; r<= (N - i); r++){
	// 		int ir2index = misc_ir2index_cpp(i,r,N);
	// 		index2r[ir2index] = r;
	// 		index2i[ir2index] = i;
	// 	}
	// }
	// instantiate the logfhar matrix
	NumericMatrix logfbar_matrix(irspace_size,irspace_size);
	std::fill(logfbar_matrix.begin(), logfbar_matrix.end(), R_NegInf);
	// compute the logfbar matrix
	// for(int prev_state = 0; prev_state < irspace_size; prev_state++){
	// 	for(int next_state = 0; next_state < irspace_size; next_state++){
	// 		logfbar_matrix(prev_state,next_state) = logfbar_cpp(lambda_bar,gamma_bar,N,index2i[prev_state],index2r[prev_state],index2i[next_state],index2r[next_state]);
	// 	}
	// }
	// can be very careful with the constraints
	for(int iprev = 0; iprev <= N; iprev++){
		for(int rprev = 0; rprev <= (N-iprev); rprev++){
			int prev_state = misc_ir2index_cpp(iprev, rprev, N);
			for(int rnext = rprev; rnext <= (iprev + rprev); rnext++){
				for(int inext = (iprev - rnext + rprev); inext <= (N - rnext) ; inext++){
					int next_state = misc_ir2index_cpp(inext,rnext, N);
					logfbar_matrix(prev_state,next_state) = logfbar_cpp(lambda_bar,gamma_bar,N,iprev,rprev,inext,rnext);
				}
			}
		}
	}
	// approximate dynamic programming runs backwards
	for(int i = y[num_observations - 1]; i <= N; i++){
		for (int r = 0; r <= (N - i); r++){
			policy_matrix(misc_ir2index_cpp(i,r,N), num_observations - 1) = R::dbinom(y[num_observations - 1], i, rho, TRUE);
		}
	}
	NumericVector loggt(N + 1);
	for(int t = num_observations - 2; t>=0; t--){
		std::fill(loggt.begin(),loggt.end(), R_NegInf);
		for(int i = y[t]; i <= N; i++){
			loggt[i] = R::dbinom(y[t],i,rho,TRUE);
			for(int r = 0; r <= (N - i); r++){
				int ir_index = misc_ir2index_cpp(i,r,N);
				policy_matrix(ir_index, t) = loggt[i] + lw_logsum_cpp(policy_matrix(_, t + 1) + logfbar_matrix(ir_index,_));
			}
		}
	}
	// // current_logpolicy is not necessary and wastes space
	// NumericVector current_logpolicy(irspace_size);
	// std::fill(current_logpolicy.begin(), current_logpolicy.end(),R_NegInf);
	// for(int i = y[num_observations - 1]; i <= N; i++){
	// 	for (int r = 0; r <= (N - i); r++){
	// 		current_logpolicy[misc_ir2index_cpp(i,r,N)] = R::dbinom(y[num_observations - 1], i, rho, TRUE);
	// 	}
	// }
	// policy_matrix(_, num_observations - 1) = current_logpolicy;
	// NumericVector loggt(N + 1);
	// for(int t = num_observations - 2; t>=0; t--){
	// 	std::fill(current_logpolicy.begin(), current_logpolicy.end(),R_NegInf);
	// 	std::fill(loggt.begin(),loggt.end(), R_NegInf);
	// 	for(int i = y[t]; i <= N; i++){
	// 		loggt[i] = R::dbinom(y[t],i,rho,TRUE);
	// 		for(int r = 0; r <= (N - i); r++){
	// 			int ir_index = misc_ir2index_cpp(i,r,N);
	// 			current_logpolicy[ir_index] = loggt[i] + lw_logsum_cpp(policy_matrix(_, t + 1) + logfbar_matrix(ir_index,_));
	// 		}
	// 	}
	// 	policy_matrix(_,t) = current_logpolicy;
	// }
	return policy_matrix;
}