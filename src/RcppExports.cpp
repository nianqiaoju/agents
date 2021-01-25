// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// get_runifs
NumericVector get_runifs(const int n);
RcppExport SEXP _agents_get_runifs(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(get_runifs(n));
    return rcpp_result_gen;
END_RCPP
}
// idchecking_cpp
LogicalVector idchecking_cpp(const int sum_x, const NumericVector& alpha, const NumericVector& random_number);
RcppExport SEXP _agents_idchecking_cpp(SEXP sum_xSEXP, SEXP alphaSEXP, SEXP random_numberSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type sum_x(sum_xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type random_number(random_numberSEXP);
    rcpp_result_gen = Rcpp::wrap(idchecking_cpp(sum_x, alpha, random_number));
    return rcpp_result_gen;
END_RCPP
}
// logdbern_sum_cpp
double logdbern_sum_cpp(const NumericVector& alpha, const LogicalVector& x);
RcppExport SEXP _agents_logdbern_sum_cpp(SEXP alphaSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const LogicalVector& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(logdbern_sum_cpp(alpha, x));
    return rcpp_result_gen;
END_RCPP
}
// logdpoisbinom_cpp
NumericVector logdpoisbinom_cpp(const NumericVector& alpha);
RcppExport SEXP _agents_logdpoisbinom_cpp(SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(logdpoisbinom_cpp(alpha));
    return rcpp_result_gen;
END_RCPP
}
// lw_logsum_cpp
double lw_logsum_cpp(NumericVector lw);
RcppExport SEXP _agents_lw_logsum_cpp(SEXP lwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type lw(lwSEXP);
    rcpp_result_gen = Rcpp::wrap(lw_logsum_cpp(lw));
    return rcpp_result_gen;
END_RCPP
}
// lw_normalize_cpp
NumericVector lw_normalize_cpp(const NumericVector& lw);
RcppExport SEXP _agents_lw_normalize_cpp(SEXP lwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type lw(lwSEXP);
    rcpp_result_gen = Rcpp::wrap(lw_normalize_cpp(lw));
    return rcpp_result_gen;
END_RCPP
}
// metropolis_condbern_cpp
LogicalVector metropolis_condbern_cpp(int sum_x, NumericVector alpha, int num_mcmc);
RcppExport SEXP _agents_metropolis_condbern_cpp(SEXP sum_xSEXP, SEXP alphaSEXP, SEXP num_mcmcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type sum_x(sum_xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type num_mcmc(num_mcmcSEXP);
    rcpp_result_gen = Rcpp::wrap(metropolis_condbern_cpp(sum_x, alpha, num_mcmc));
    return rcpp_result_gen;
END_RCPP
}
// multinomial_cpp
int multinomial_cpp(const NumericVector& logweights, double uniform);
RcppExport SEXP _agents_multinomial_cpp(SEXP logweightsSEXP, SEXP uniformSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type logweights(logweightsSEXP);
    Rcpp::traits::input_parameter< double >::type uniform(uniformSEXP);
    rcpp_result_gen = Rcpp::wrap(multinomial_cpp(logweights, uniform));
    return rcpp_result_gen;
END_RCPP
}
// multinomial_resampling_cpp
IntegerVector multinomial_resampling_cpp(const NumericVector& weights, int ndraws, const NumericVector& rand);
RcppExport SEXP _agents_multinomial_resampling_cpp(SEXP weightsSEXP, SEXP ndrawsSEXP, SEXP randSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< int >::type ndraws(ndrawsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type rand(randSEXP);
    rcpp_result_gen = Rcpp::wrap(multinomial_resampling_cpp(weights, ndraws, rand));
    return rcpp_result_gen;
END_RCPP
}
// rcppeigen_hello_world
Eigen::MatrixXd rcppeigen_hello_world();
RcppExport SEXP _agents_rcppeigen_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcppeigen_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcppeigen_outerproduct
Eigen::MatrixXd rcppeigen_outerproduct(const Eigen::VectorXd& x);
RcppExport SEXP _agents_rcppeigen_outerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcppeigen_outerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcppeigen_innerproduct
double rcppeigen_innerproduct(const Eigen::VectorXd& x);
RcppExport SEXP _agents_rcppeigen_innerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcppeigen_innerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcppeigen_bothproducts
Rcpp::List rcppeigen_bothproducts(const Eigen::VectorXd& x);
RcppExport SEXP _agents_rcppeigen_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcppeigen_bothproducts(x));
    return rcpp_result_gen;
END_RCPP
}
// misc_ir2index_cpp
int misc_ir2index_cpp(int i, int r, int N);
RcppExport SEXP _agents_misc_ir2index_cpp(SEXP iSEXP, SEXP rSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(misc_ir2index_cpp(i, r, N));
    return rcpp_result_gen;
END_RCPP
}
// logfbar_cpp
double logfbar_cpp(double lambda_bar, double gamma_bar, int N, int i0, int r0, int i1, int r1);
RcppExport SEXP _agents_logfbar_cpp(SEXP lambda_barSEXP, SEXP gamma_barSEXP, SEXP NSEXP, SEXP i0SEXP, SEXP r0SEXP, SEXP i1SEXP, SEXP r1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type lambda_bar(lambda_barSEXP);
    Rcpp::traits::input_parameter< double >::type gamma_bar(gamma_barSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type i0(i0SEXP);
    Rcpp::traits::input_parameter< int >::type r0(r0SEXP);
    Rcpp::traits::input_parameter< int >::type i1(i1SEXP);
    Rcpp::traits::input_parameter< int >::type r1(r1SEXP);
    rcpp_result_gen = Rcpp::wrap(logfbar_cpp(lambda_bar, gamma_bar, N, i0, r0, i1, r1));
    return rcpp_result_gen;
END_RCPP
}
// sir_bif_policy_matrix_cpp
NumericMatrix sir_bif_policy_matrix_cpp(double lambda_bar, double gamma_bar, IntegerVector y, int N, double rho);
RcppExport SEXP _agents_sir_bif_policy_matrix_cpp(SEXP lambda_barSEXP, SEXP gamma_barSEXP, SEXP ySEXP, SEXP NSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type lambda_bar(lambda_barSEXP);
    Rcpp::traits::input_parameter< double >::type gamma_bar(gamma_barSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(sir_bif_policy_matrix_cpp(lambda_bar, gamma_bar, y, N, rho));
    return rcpp_result_gen;
END_RCPP
}
// sis_apf_exact_full_cpp
double sis_apf_exact_full_cpp(const IntegerVector& y, const NumericVector& alpha0, const NumericVector& lambda, const NumericVector& gamma, const double rho, const int num_particles, const double threshold, const int population_size);
RcppExport SEXP _agents_sis_apf_exact_full_cpp(SEXP ySEXP, SEXP alpha0SEXP, SEXP lambdaSEXP, SEXP gammaSEXP, SEXP rhoSEXP, SEXP num_particlesSEXP, SEXP thresholdSEXP, SEXP population_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type alpha0(alpha0SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< const double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const int >::type num_particles(num_particlesSEXP);
    Rcpp::traits::input_parameter< const double >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< const int >::type population_size(population_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(sis_apf_exact_full_cpp(y, alpha0, lambda, gamma, rho, num_particles, threshold, population_size));
    return rcpp_result_gen;
END_RCPP
}
// sis_forward_algorithm_cpp
NumericMatrix sis_forward_algorithm_cpp(NumericVector logf0, NumericMatrix logdtransition, NumericVector y, NumericVector all_sum_x, double rho);
RcppExport SEXP _agents_sis_forward_algorithm_cpp(SEXP logf0SEXP, SEXP logdtransitionSEXP, SEXP ySEXP, SEXP all_sum_xSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type logf0(logf0SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type logdtransition(logdtransitionSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type all_sum_x(all_sum_xSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(sis_forward_algorithm_cpp(logf0, logdtransition, y, all_sum_x, rho));
    return rcpp_result_gen;
END_RCPP
}
// sis_get_alpha_full_cpp
NumericVector sis_get_alpha_full_cpp(const LogicalVector& xx, const NumericVector& lambda, const NumericVector& gamma);
RcppExport SEXP _agents_sis_get_alpha_full_cpp(SEXP xxSEXP, SEXP lambdaSEXP, SEXP gammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const LogicalVector& >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type gamma(gammaSEXP);
    rcpp_result_gen = Rcpp::wrap(sis_get_alpha_full_cpp(xx, lambda, gamma));
    return rcpp_result_gen;
END_RCPP
}
// sis_xx_gibbs_blocked_full_cpp
void sis_xx_gibbs_blocked_full_cpp(LogicalMatrix xx, const IntegerVector& y, const NumericVector& alpha0, const NumericVector& lambda, const NumericVector& gamma, const double rho, const int block_size, const LogicalMatrix& state_space);
RcppExport SEXP _agents_sis_xx_gibbs_blocked_full_cpp(SEXP xxSEXP, SEXP ySEXP, SEXP alpha0SEXP, SEXP lambdaSEXP, SEXP gammaSEXP, SEXP rhoSEXP, SEXP block_sizeSEXP, SEXP state_spaceSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< LogicalMatrix >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type alpha0(alpha0SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< const double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const int >::type block_size(block_sizeSEXP);
    Rcpp::traits::input_parameter< const LogicalMatrix& >::type state_space(state_spaceSEXP);
    sis_xx_gibbs_blocked_full_cpp(xx, y, alpha0, lambda, gamma, rho, block_size, state_space);
    return R_NilValue;
END_RCPP
}
// sis_xx_gibbs_singlesite_full_cpp
void sis_xx_gibbs_singlesite_full_cpp(LogicalMatrix xx, const IntegerVector& y, const NumericVector& alpha0, const NumericVector& lambda, const NumericVector& gamma, const double rho);
RcppExport SEXP _agents_sis_xx_gibbs_singlesite_full_cpp(SEXP xxSEXP, SEXP ySEXP, SEXP alpha0SEXP, SEXP lambdaSEXP, SEXP gammaSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< LogicalMatrix >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type alpha0(alpha0SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< const double >::type rho(rhoSEXP);
    sis_xx_gibbs_singlesite_full_cpp(xx, y, alpha0, lambda, gamma, rho);
    return R_NilValue;
END_RCPP
}
// sis_xx_gibbs_swap_full_cpp
void sis_xx_gibbs_swap_full_cpp(LogicalMatrix xx, const NumericVector& alpha0, const NumericVector& lambda, const NumericVector& gamma);
RcppExport SEXP _agents_sis_xx_gibbs_swap_full_cpp(SEXP xxSEXP, SEXP alpha0SEXP, SEXP lambdaSEXP, SEXP gammaSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< LogicalMatrix >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type alpha0(alpha0SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type gamma(gammaSEXP);
    sis_xx_gibbs_swap_full_cpp(xx, alpha0, lambda, gamma);
    return R_NilValue;
END_RCPP
}
// static_xx_gibbs_cpp
LogicalVector static_xx_gibbs_cpp(LogicalVector xx_previous, NumericVector alpha, double rho, int y);
RcppExport SEXP _agents_static_xx_gibbs_cpp(SEXP xx_previousSEXP, SEXP alphaSEXP, SEXP rhoSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< LogicalVector >::type xx_previous(xx_previousSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< int >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(static_xx_gibbs_cpp(xx_previous, alpha, rho, y));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_condbernmcmc_module();

static const R_CallMethodDef CallEntries[] = {
    {"_agents_get_runifs", (DL_FUNC) &_agents_get_runifs, 1},
    {"_agents_idchecking_cpp", (DL_FUNC) &_agents_idchecking_cpp, 3},
    {"_agents_logdbern_sum_cpp", (DL_FUNC) &_agents_logdbern_sum_cpp, 2},
    {"_agents_logdpoisbinom_cpp", (DL_FUNC) &_agents_logdpoisbinom_cpp, 1},
    {"_agents_lw_logsum_cpp", (DL_FUNC) &_agents_lw_logsum_cpp, 1},
    {"_agents_lw_normalize_cpp", (DL_FUNC) &_agents_lw_normalize_cpp, 1},
    {"_agents_metropolis_condbern_cpp", (DL_FUNC) &_agents_metropolis_condbern_cpp, 3},
    {"_agents_multinomial_cpp", (DL_FUNC) &_agents_multinomial_cpp, 2},
    {"_agents_multinomial_resampling_cpp", (DL_FUNC) &_agents_multinomial_resampling_cpp, 3},
    {"_agents_rcppeigen_hello_world", (DL_FUNC) &_agents_rcppeigen_hello_world, 0},
    {"_agents_rcppeigen_outerproduct", (DL_FUNC) &_agents_rcppeigen_outerproduct, 1},
    {"_agents_rcppeigen_innerproduct", (DL_FUNC) &_agents_rcppeigen_innerproduct, 1},
    {"_agents_rcppeigen_bothproducts", (DL_FUNC) &_agents_rcppeigen_bothproducts, 1},
    {"_agents_misc_ir2index_cpp", (DL_FUNC) &_agents_misc_ir2index_cpp, 3},
    {"_agents_logfbar_cpp", (DL_FUNC) &_agents_logfbar_cpp, 7},
    {"_agents_sir_bif_policy_matrix_cpp", (DL_FUNC) &_agents_sir_bif_policy_matrix_cpp, 5},
    {"_agents_sis_apf_exact_full_cpp", (DL_FUNC) &_agents_sis_apf_exact_full_cpp, 8},
    {"_agents_sis_forward_algorithm_cpp", (DL_FUNC) &_agents_sis_forward_algorithm_cpp, 5},
    {"_agents_sis_get_alpha_full_cpp", (DL_FUNC) &_agents_sis_get_alpha_full_cpp, 3},
    {"_agents_sis_xx_gibbs_blocked_full_cpp", (DL_FUNC) &_agents_sis_xx_gibbs_blocked_full_cpp, 8},
    {"_agents_sis_xx_gibbs_singlesite_full_cpp", (DL_FUNC) &_agents_sis_xx_gibbs_singlesite_full_cpp, 6},
    {"_agents_sis_xx_gibbs_swap_full_cpp", (DL_FUNC) &_agents_sis_xx_gibbs_swap_full_cpp, 4},
    {"_agents_static_xx_gibbs_cpp", (DL_FUNC) &_agents_static_xx_gibbs_cpp, 4},
    {"_rcpp_module_boot_condbernmcmc_module", (DL_FUNC) &_rcpp_module_boot_condbernmcmc_module, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_agents(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
