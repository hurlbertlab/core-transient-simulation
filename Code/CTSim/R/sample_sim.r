#' Sample simulation with imperfect detection
#'
#' Stochastically samples a simulation using a fixed detection probability.
#'
#' The function samples the abundances of species at 
#' a set of locations and timepoints (given in \code{abuns}) using 
#' detection probabilities in \code{probs}. These are the probabilities
#' that an individual of a given species will be observed and can be set to
#' be the same for all species by only providing a single number. It returns
#' an array with the same dimensions as \code{abun} containing either the
#' observed abundances or presence/absence (depending on \code{return}). 
#'
#' @param abuns (required) array of abundance profiles (as returned by 
#' 	\code{\link{calc_abun_profile}})
#' @param probs probability of observing an individual. Can be either a 
#' 	single detection probability for all species, or a vector with different 
#' 	probabilies for each species. Defaults to 1.
#' @param return a character string indicating whether species abundances 
#' 	(\code{'abundance'}) or presence (\code{'presence'}) should be returned. 
#' 	Default is abundance.
#' @return an array with the same dimensions as \code{abuns}
#'
#' @import abind
#' @export

sample_sim = function(abuns, probs = NULL, return='abundance'){
	
	# Drop empty spaces from abuns, if present.
	abuns = abuns[,dimnames(abuns)[[2]]!='0',,drop=FALSE]

	# Determine number of species
	N_S = dim(abuns)[2]

	# If not specified, detection is 1 and the original abundance profiles are returned
	if(is.null(probs)){
		obs = abuns
		if(return=='presence') obs = abuns>0
	} else {
		# Get dimensions of abuns
		ndim = length(dim(abuns))  
		
		# Make vector of detectabilities if only one specified.
		if(length(probs)==1) probs = rep(probs, N_S)
		
		# If species presence to be returned		
		if(return=='presence'){

			# Calculate probability of observing species (1 - P(not observing))
			# Returns an array: [timepoints, sites, species]
			P = sapply(1:N_S, function(sp){
				apply(abuns[,sp,,drop=FALSE], 1:ndim, function(x) 1 - (1-probs[sp])^x )
			}, simplify='array')

			# Stochastically determine which species are observed
			rands = array(runif(length(P)), dim=dim(P))
			obs = rands <= P
		}
		
		if(return=='abundance'){
			
			# Stochasitically determine which individuals are observed
			obs = sapply(1:N_S, function(sp){
				apply(abuns[,sp,,drop=FALSE], 1:ndim, function(x){
					if(x>0){
						sum(runif(x) <= probs[sp])
					} else {
						0
					}
				})
			}, simplify='array')
		}
		
		# Rearrange dimensions to match abuns
		obs = abind::adrop(obs, 2)
		if(ndim==3) obs = aperm(obs, c(1,3,2))
	}

	# Return observed sample
	obs
} 
