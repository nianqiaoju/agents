/*
 * This implements an MCMC algorithm to sample conditional
 * Bernoulli distribution, which is the distribution of N Bernoulli variables X = (X^1,...,X^N)
 * each with a possibly different probability so that P(X^n = 1) = \alpha^n, 
 * conditioned on the sum \sum_{n=1}^N X^n being equal to some integer in {0,...,N}.
 * 
 * The input is the vector 'logOdds' of log-odds = log(\alpha^n/(1-\alpha^n)),
 * as well as the sum on which to condition 'sumx'
 * 
 * The implementation involves O(1) operations per iteration,
 * and also includes a coupling that costs O(1) operations per iteration.
 * 
 */
#include <Rcpp.h>
#include <vector>
using namespace std;
using namespace Rcpp;

// the following function is taken from
// https://stackoverflow.com/questions/9218724/get-random-element-and-remove-it
// and returns+removes an element from a vector, after swapping it with the last element,
// which makes it O(1) instead of the method 'erase' which costs O(v.size())
template <typename T>
T remove_at(std::vector<T>&v, typename std::vector<T>::size_type n)
{
  T ans = std::move_if_noexcept(v[n]);
  v[n] = std::move_if_noexcept(v.back());
  v.pop_back();
  return ans;
}

// sample uniformly in {0,...,n}
// should we add some rng scope GetRNGstate PutRNGstate etc?
inline int randWrapper(const int n) { return floor(unif_rand()*n); }

// class to generate Markov chains targeting the conditional Bernoulli
// distribution on {0,1}^N also contains function to implement coupling of chain
// in O(1) operations per iteration.

class CondBernMCMC
{
public:
  // constructor using log odds = log(\alpha/(1-\alpha)) where \alpha is the N-vector of probabilities
  // and sumx is the sum of the Bernoulli variables, to condition on
  CondBernMCMC(NumericVector logOdds, int sumx);
  // destructor
  ~CondBernMCMC();
  // sample initial configuration in O(N) operations
  void init();
  // perform one step of MCMC in O(1) operations
  void iteration();
  // return N-vector, with binary entries summing to 'sumx'
  IntegerVector getState();
  // for coupled chains, sample initial configurations in O(N) operations
  void coupledinit();
  // perform one step of coupled MCMC in O(1) operations
  // the argument specifies if both chains should be updated, or just the first
  // which can be useful for implementation of unbiased MCMC estimators
  void couplediteration(bool updateboth);
  // return two N-vectors as two columns of a binary matrix, each column summing so 'sumx'
  IntegerMatrix getTwoStates();
  // return a boolean indicating whether two states are identical
  bool identicalStates();
  // size of the vector X
  int N;
  // log odds of the N Bernoulli variables
  NumericVector logOdds;
  // sum of vector X on which to condition
  int sumx;
  // set of indices n such that X^n = 0 
  vector<int> s0;
  // set of indices n such that X^n = 1 
  vector<int> s1;
  // number of accepts & iterations in the single MCMC algorithm
  // incremented in this->iteration not by this->couplediteration
  int naccepts, niteration;
  // 'sij' is set of indices n such that X^n = i, \tilde{X}^n = j
  vector<int> s00, s11, s01, s10;
  // some temporary objects used in this->iteration and this->couplediteration
  // these refer to indices within s00,s11,s01,s10
  int i0, i1, i0tilde, i1tilde;
  // these indicate which sets the indices above are indexing
  int groupi0, groupi1, groupi0tilde, groupi1tilde;
  // these refer to indices in the N-vectors X and \tilde{X}
  int i0inx, i1inx, i0tildeinx, i1tildeinx;
  // this is a temporary integer used to perform 'swaps'
  int tmp;
  // these are log of acceptance rates and uniform(0,1) draws
  double logar, logartilde, loguar;
};

// constructor using logOdds and sumx
CondBernMCMC::CondBernMCMC(NumericVector logOdds, int sumx) : logOdds(logOdds), sumx(sumx) {
  // could check here that sumx is <= N 
  N = logOdds.length();
}

// destructor
CondBernMCMC::~CondBernMCMC(){
}

// initialize state before running MCMC using this->iteration
void CondBernMCMC::init(){
  GetRNGstate();
  // reset s0, s1
  s0.clear(); s1.clear();
  s0.reserve(N); s1.reserve(N);
  // random permutation of indices in 0,...,N-1
  IntegerVector indices = seq(0, N-1);
  std::random_shuffle(indices.begin(), indices.end(), randWrapper);
  // first 'sumx' shuffled indices are set to one in vector X
  for (int i = 0; i < N; i++){
    if (i < sumx){
      s1.push_back(indices(i));
    } else {
      s0.push_back(indices(i));
    }
  }
  // reset counters
  naccepts = 0;
  niteration = 0;
  PutRNGstate();
}

// perform one MCMC step in O(1) operations
void CondBernMCMC::iteration(){
  GetRNGstate();
  niteration ++;
  // sample i0 in s0 and i1 in s1
  i0 = randWrapper(N-sumx);
  i1 = randWrapper(sumx);
  // log accept ratio 
  logar = logOdds(s0.at(i0)) - logOdds(s1.at(i1));
  loguar = log(unif_rand());
  // comparison done on log scale
  if (loguar < logar){
    // accept: swap elements
    tmp = s0.at(i0);
    s0.at(i0) = s1.at(i1);
    s1.at(i1) = tmp;
    naccepts ++;
  } else {
    // do nothing
  }
  PutRNGstate();
}


// initialize coupled states before running MCMC using this->couplediteration
void CondBernMCMC::coupledinit(){
  GetRNGstate();
  // clear memory for s00, s01, s10, s11
  s00.clear();
  s11.clear();
  s01.clear();
  s10.clear();
  s00.reserve(N);
  s11.reserve(N);
  s01.reserve(N);
  s10.reserve(N);
  // 
  IntegerVector indices =  seq(0, N-1);
  IntegerVector indicestilde = seq(0, N-1);
  // random permutation of indices 0,...,M-1
  std::random_shuffle(indices.begin(), indices.end(), randWrapper);
  std::random_shuffle(indicestilde.begin(), indicestilde.end(), randWrapper);
  // first 'sumx' shuffled indices are set to one in X and \tilde{X}
  IntegerVector x = rep(0,N);
  IntegerVector xtilde = rep(0,N);
  for (int i = 0; i < sumx; i++){
    x(indices(i)) = 1;
    xtilde(indicestilde(i)) = 1;
  }
  // fill up s00, s11, s01, s10 based on x and xtilde
  for (int i = 0; i < N; i++){
    if (x(i) == 0 && xtilde(i) == 0){
      s00.push_back(i);
    }
    if (x(i) == 1 && xtilde(i) == 1){
      s11.push_back(i);
    }
    if (x(i) == 0 && xtilde(i) == 1){
      s01.push_back(i);
    }
    if (x(i) == 1 && xtilde(i) == 0){
      s10.push_back(i);
    }
  }
  PutRNGstate();
}


// perform one coupled MCMC step in O(1) operations; main idea of the
// implementation is to keep track of indices within 's' sets (s00, s11, s01,
// s10), indices in 'x' N-vector and 'group' associated with each index
// (1,2,3,4) respectively 
void CondBernMCMC::couplediteration(bool updateboth = true){
  GetRNGstate();
  // groups are 1,2,3,4 for {00,11,01,10}
  // sample i0, i0tilde from max coupling
  // probability of sampling identical index
  double zerocommonmass = (double) s00.size() / (N - sumx);
  double ucommonzero = unif_rand();
  bool commonzero = (ucommonzero < zerocommonmass);
  if (commonzero){
    // sample in s00
    i0 = floor(unif_rand() * s00.size());
    i0inx = s00.at(i0);
    groupi0 = 1;
    i0tildeinx = i0inx;
    groupi0tilde = 1;
    i0tilde = i0;
    // start computation of accept ratio
    logar = logOdds(s00.at(i0));
    logartilde = logOdds(s00.at(i0tilde));
  } else {
    // sample i0 in s01
    i0 = floor(unif_rand() * s01.size());
    i0inx = s01.at(i0);
    groupi0 = 3;
    // sample i0tilde in s10
    i0tilde = floor(unif_rand() * s10.size());
    i0tildeinx = s10.at(i0tilde);
    groupi0tilde = 4;
    logar = logOdds(s01.at(i0));
    logartilde = logOdds(s10.at(i0tilde));
  }
  // sample i1, i1tilde from max coupling
  double onecommonmass = (double) s11.size() / (sumx);
  double ucommonone = unif_rand();
  bool commonone = (ucommonone < onecommonmass);
  if (commonone){
    // sample in s11
    i1 = floor(unif_rand() * s11.size());
    groupi1 = 2;
    i1inx = s11.at(i1);
    i1tilde = i1;
    i1tildeinx = i1inx;
    groupi1tilde = 2;
    logar -= logOdds(s11.at(i1));
    logartilde -= logOdds(s11.at(i1tilde));
  } else {
    // sample i1 in s10 and i1tilde in s01
    i1 = floor(unif_rand() * s10.size());
    i1inx = s10.at(i1);
    groupi1 = 4;
    i1tilde = floor(unif_rand() * s01.size());
    i1tildeinx = s01.at(i1tilde);
    groupi1tilde = 3;
    logar -= logOdds(s10.at(i1));
    logartilde -= logOdds(s01.at(i1tilde));
  }
  // use common uniform to accept or not
  loguar = log(unif_rand());
  bool accept = (loguar < logar);
  bool accepttilde = (loguar < logartilde);
  if (accept){
    // update s00,s11,s01,s10
    // swap i0 in groupi0 with i1 in groupi1
    // needs to keep track carefully of the sets in which
    // i0tilde and i1tilde belong, 
    // otherwise we wouldn't be able to perform the 'accept' step
    // of the second chain in O(1) operations
    if (groupi0 == 1 && groupi1 == 2){
      // swap i0 and i1
      // taking s00(i0) from s00 and putting it in s01 
      // taking s11(i1) from s11 and putting it in s11 
      remove_at(s00, i0);
      remove_at(s11, i1);
      // s11.remove_at(s11.begin()+i1);
      s01.push_back(i1inx);
      s10.push_back(i0inx);
      i1tilde = s01.size()-1;
      i0tilde = s10.size()-1;
      groupi0tilde = 4;
      groupi1tilde = 3;
      // note i0tildeinx and i1tildeinx aren't changed 
    }
    if (groupi0 == 1 && groupi1 == 4){
      // taking s00(i0) from s00 and putting it in s10
      // taking s10(i1) from s10 and putting it in s00
      // s00.erase(s00.begin()+i0);
      // s10.erase(s10.begin()+i1);
      remove_at(s00, i0);
      remove_at(s10, i1);
      s00.push_back(i1inx);
      s10.push_back(i0inx);
      i0tilde = s10.size()-1;
      groupi0tilde = 4;
      // note i1tilde is necessarily in group 3 here and is unchanged
    }
    if (groupi0 == 3 && groupi1 == 2){
      // s01.erase(s01.begin()+i0);
      // s11.erase(s11.begin()+i1);
      remove_at(s01, i0);
      remove_at(s11, i1);
      s01.push_back(i1inx);
      s11.push_back(i0inx);
      i1tilde = s01.size()-1;
      groupi1tilde = 3;
      // i0tilde is necessarily in group 4 here and is unchanged
    }
    if (groupi0 == 3 && groupi1 == 4){
      // s01.erase(s01.begin()+i0);
      // s10.erase(s10.begin()+i1);
      remove_at(s01, i0);
      remove_at(s10, i1);
      s00.push_back(i1inx);
      s11.push_back(i0inx);
      // now i0tildeinx could be identical to i1inx, 
      // in which case we need to change groupi0tilde
      if (i0tildeinx == i1inx){
        i0tilde = s00.size()-1;
        groupi0tilde = 1;
      }
      // and i1tildeinx could be identical to i0inx, 
      // in which case we need to change groupi1tilde
      if (i1tildeinx == i0inx){
        i1tilde = s11.size()-1;
        groupi1tilde = 2;
      }      
    }
    // could implement this with 'switch case'...
  } else {
    // do nothing
  }
  // 
  if (updateboth && accepttilde){
    // update s00,s11,s01,s10
    // swap i0tilde in groupi0tilde with i1tilde in groupi1tild
    if (groupi0tilde == 1 && groupi1tilde == 2){
      // s00.erase(s00.begin()+i0tilde);
      // s11.erase(s11.begin()+i1tilde);
      remove_at(s00, i0tilde);
      remove_at(s11, i1tilde);
      
      s01.push_back(i0tildeinx);
      s10.push_back(i1tildeinx);
    }
    if (groupi0tilde == 1 && groupi1tilde == 3){
      // s00.erase(s00.begin()+i0tilde);
      // s01.erase(s01.begin()+i1tilde);
      remove_at(s00, i0tilde);
      remove_at(s01, i1tilde);
      
      s00.push_back(i1tildeinx);
      s01.push_back(i0tildeinx);
    }
    if (groupi0tilde == 4 && groupi1tilde == 2){
      // s11.erase(s11.begin()+i1tilde);
      // s10.erase(s10.begin()+i0tilde);
      remove_at(s11, i1tilde);
      remove_at(s10, i0tilde);
      
      s11.push_back(i0tildeinx);
      s10.push_back(i1tildeinx);
    }
    if (groupi0tilde == 4 && groupi1tilde == 3){
      // s10.erase(s10.begin()+i0tilde);
      // s01.erase(s01.begin()+i1tilde);
      remove_at(s10, i0tilde);
      remove_at(s01, i1tilde);
      
      s11.push_back(i0tildeinx);
      s00.push_back(i1tildeinx);
    }
  } else {
    // do nothing
  }
  PutRNGstate();
}


// retrieve N-state of X, necessarily in either N-I or I operations;
// here we choose I but we could take min(I, N-I) to speed things up 
IntegerVector CondBernMCMC::getState(){
  IntegerVector x = rep(0, N);
  for (int i = 0; i < s1.size(); i++){
    x(s1.at(i)) = 1;
  }
  return x;
}

// are the two states identical?
bool CondBernMCMC::identicalStates(){
  bool identical = false;
  if (s00.size() == N-sumx){
    // then in principle the states must be equal with necessarily
    // s11.size() == N, s01.size() == 0, s10.size() == 0
    identical = true;
  }
  return(identical);
}
  
// retrieve N-state of X and \tilde{X} in O(N) operations
// here we choose to do it in O(I) operations but it could
// be done in O(min(I,N-I)) operations
IntegerMatrix CondBernMCMC::getTwoStates(){
  IntegerMatrix xxtilde(N,2);
  // fill with zero
  std::fill(xxtilde.begin(), xxtilde.end(), 0);
  // and put one in adequate places
  for (int i = 0; i < s01.size(); i++){
    xxtilde(s01.at(i),1) = 1;
  }
  for (int i = 0; i < s10.size(); i++){
    xxtilde(s10.at(i),0) = 1;
  }
  for (int i = 0; i < s11.size(); i++){
    xxtilde(s11.at(i),0) = 1;
    xxtilde(s11.at(i),1) = 1;
  }
  return xxtilde;
}

// The following exports the class CondBernMCMC
// so that an R user can type:
// cbm <- new(CondBernMCMC, logOdds, sum_x)
// cbm$coupledinit()
// cbm$couplediteration(TRUE)
// cbm$getTwoStates()

RCPP_MODULE(condbernmcmc_module) {
  class_<CondBernMCMC>( "CondBernMCMC" )
  .constructor<NumericVector,int>()
  .field( "N", &CondBernMCMC::N )
  .field( "naccepts", &CondBernMCMC::naccepts)
  .field( "niteration", &CondBernMCMC::niteration)
  .method( "init", &CondBernMCMC::init )
  .method( "iteration", &CondBernMCMC::iteration )
  .method( "getState", &CondBernMCMC::getState )
  .method( "coupledinit", &CondBernMCMC::coupledinit )
  .method( "couplediteration", &CondBernMCMC::couplediteration )
  .method( "getTwoStates", &CondBernMCMC::getTwoStates )
  .method( "identicalStates", &CondBernMCMC::identicalStates )
  ;
}
