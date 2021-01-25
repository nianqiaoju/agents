#' @title a coupling kernel that leaves conditional-bernoullis invariant
#' @description discription of the algorithm in page 21 of Jacob et al. (2017).
#' for one chain, it samples index of a 0 entry in X uniformly and index of a 1 entry uniformly
#' it proposes a switch of the 1 and 0 and accepts the swap with probability ratio p(proposal) / p(current)
#' @param xx1 current state of chain1
#' @param xx2 current state of chain2
#' @param sum_x1 sum of xx1
#' @param sum_x2 sum of xx2
#' @param alpha1 vector of probabilities for chain1
#' @param alpha2 vector of probabilities for chain2
#' @return a matrix row 1 is new xx1 and row 2 is new xx2

couple_condbern_kernel<- function(xx1,xx2,sum_x1,sum_x2,alpha1,alpha2){
  ## propose i0 such that xx1[i0] = 0 
  ## propose i1 such that xx1[i1] = 1
  ## propose j0 such that xx2[j0] = 0 , j1 such that xx2[j1] = 1
  indices <- unbiasedmcmc::coupled_pairs01(xx1,xx2)
  ## indices is c(i0,j0,i1,j1)
  ## probability of accept swap on chain 1 
  accept1 <- (1 - alpha1[indices[3]]) * (alpha1[indices[1]]) / (1 - alpha1[indices[1]]) / alpha1[indices[3]] 
  ## probability of accept swap on chain 2 
  accept2 <- (1 - alpha2[indices[4]]) * (alpha2[indices[2]]) / (1 - alpha2[indices[2]]) / alpha2[indices[4]] 
  rU <- runif(1)
  if (rU < accept1){
    xx1[indices[3]] <- 0 
    xx1[indices[1]] <- 1
  }
  if (rU < accept2){
    xx2[indices[2]] <- 1
    xx2[indices[4]] <- 0 
  }
  return(rbind(xx1,xx2))
}
