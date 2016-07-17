#' Generate propagules
#'
#' Generates a pool of propagules produced by a community during one time step.
#'
#' \code{comm} is a vector of integers indicating which species are present.
#' The order of \code{b_rates} designates which species each birth rate
#' applies to. For example, \code{b_rates[1]} is the birth rate of species 1.
#'
#' @param comm (required) vector (or list) of species present
#' 	(may contain zeros for empty spaces)
#' @param b_rates (required) vector of birth rates for each species
#' @return a vector of integers, one for each propagule

reproduce = function(comm, b_rates){

	# Convert to vetor if not already
	comm = unlist(comm)
	
	# Make pool of propagules
	propagules = sapply(comm[comm>0], function(sp) rep(sp, b_rates[sp]))
	propagules = unlist(propagules)
	
	propagules
}