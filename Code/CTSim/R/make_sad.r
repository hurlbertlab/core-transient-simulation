#' Abundance distribution generation
#'
#' Samples a given distribution for a set number of species abundances.
#' Implements several commonly used distributions (see details).
#'
#' The form of the distributition is controlled by a named list with
#' 	 the following components:
#'	\describe{
#'		\item{\code{type}}{(Required) indicates the form of the SAD.
#'			Can be: 'same', 'uniform', 'power', 'logseries', 'lognormal',
#' 			'poisson'}
#'		\item{\code{maxN}}{maximum abundance, used to find SAD parameters
#'			(if not supplied as p1 or p2)}
#'		\item{\code{P_maxN}}{probability of getting an abundance greater than maxN}
#'		\item{\code{p1}}{parameter controling the shape of the SAD
#'			(only certain types)}
#'		\item{\code{p2}}{parameter controling the shape of the SAD
#'			(only certain types)}
#'	}
#'
#' @note Code is set up so that abundance may be correlated with another
#' 	using Gaussian copulas, but this is currently not implemented. 
#' 	To add a correlation, generate x as a multivariate random normal 
#' 	variable with specified covariance matrix.
#' 
#' @param N_S number of species
#' @param distribution named list of  parameters describing the functional 
#' 	form of the distribution (SAD).
#' @return a vector of integers of length \code{N_S}
#'
#' @import poweRlaw
#' @import sads

make_sad = function(N_S, distribution){
	
	# Generate normal random variable
	x = rnorm(N_S)
	U = pnorm(x)

	# Convert to correct distribution
	if(distribution$type=='same'){
		use_p = ifelse(is.null(distribution$p1), 1, distribution$p1)
		abuns = rep(use_p, N_S)
	}

	if(distribution$type=='uniform'){
		if(is.null(distribution$maxN)){ 
			use_p = distribution$p1
		} else {
			use_p = distribution$maxN
		}
		abuns = qunif(U, 1, use_p)
	}
		
	if(distribution$type=='power'){
		if(is.null(distribution$maxN)){
			use_p = distribution$p1
		} else {
			use_p = 1/log(1/(distribution$maxN), base=distribution$P_maxN)
		}
		abuns = (1-U)^(-1/use_p)
	}
	if(distribution$type=='logseries'){
		# p2 = N, p1 = Fisher's alpha
		if(is.null(distribution$maxN)){
			use_N = distribution$p2
		} else {
			use_N = 1
			this_N = qls(1-distribution$P_maxN, use_N, distribution$p1)	
			while(this_N < distribution$maxN){
				use_N = use_N + 1
				this_N = sads::qls(1-distribution$P_maxN, use_N, distribution$p1)	
			}
			if(this_N > 10) use_N = use_N - 1
		}
		abuns = sads::qls(U, N = use_N, alpha = distribution$p1)		
	}
	if(distribution$type=='lognormal'){
		if(is.null(distribution$maxN)){
			use_mean = distribution$p1
			use_sd = distribution$p2
		} else {
			# Assume mean = 0
			use_mean = 0
			try_sig = 1
			this_N = qlnorm(1-distribution$P_maxN, use_mean, try_sig)
			i = 1
			while(abs(distribution$maxN - this_N) > 0.5){
				if(this_N > distribution$maxN){
					try_sig = try_sig - try_sig/2
					
				} else {
					try_sig = try_sig + try_sig/2
				}
				this_N = qlnorm(1-distribution$P_maxN, 0, try_sig)
			}
			use_sd = try_sig
		}

		# p1 = mean, p2 = sd
		abuns = ceiling(qlnorm(U, use_mean, use_sd))
	}
	if(distribution$type=='poisson'){
		# p1 = lambda
		if(is.null(distribution$maxN)){
			use_lamda = distribution$p1
		} else {
			use_maxN = distribution$maxN + 1
			try_lamda = use_maxN:1
			this_N = qpois(1-distribution$P_maxN, try_lamda)
			N_diffs = this_N - use_maxN
			i = 1
			while(!(0 %in% N_diffs)){
				new_start = which(abs(N_diffs)==min(abs(N_diffs)))
				try_lamda = seq(try_lamda[new_start+1], try_lamda[new_start-1], 10^(-i))
				this_N = qpois(1-distribution$P_maxN, try_lamda)
				N_diffs = this_N - use_maxN
				i = i + 1
			}
		
			use_lamda = min(try_lamda[which(N_diffs==0)])
		}

		abuns = qpois(U, use_lamda)+1
	}

	# Discretize so that it can be used as an individual fecundity rate
	abuns = round(abuns, 0)

	abuns
}
