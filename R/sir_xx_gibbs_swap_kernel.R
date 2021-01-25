#' @title Gibbs sampler for SIR model
#' @description given indices i and j , the proposed new state $Z$ is such that $Z^i \leftarrow X^j$ and $Z^j \leftarrow X^i$ and $Z^k \leftarrow X^k$ for all $k\not= i,j$. We accept the proposed state with probability \deqn{\min\{1, \frac{p_{\theta}(Z,y)}{p_{\theta}(X,y)}\}.}
#' @param xx current state x
#' @param ij a vector [i,j]
#' @param y observations
#' @param model_config

sir_xx_gibbs_swap_kernel <- function(xx, ij, y, model_config){
	ji <- rev(ij);
	logx <- sir_loglikelihood_complete(xx, y, model_config);
	zz <- xx;
	zz[ij,] <- xx[ji,];
	logz <- sir_loglikelihood_complete(zz, y, model_config);
	maxlog <- max(logx, logz);
	logx <- logx - maxlog;
	logz <- logz - maxlog;
	if(runif(1) < exp(logx)/(exp(logx) + exp(logz))){
		return(list(xx = xx, accept = FALSE));
	}else{
		return(list(xx = zz, accept = TRUE));
	}
}