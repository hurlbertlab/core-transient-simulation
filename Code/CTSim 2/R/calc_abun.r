#' Calculate global relative abundance
#'
#' Calculates the relative abundance or occupancy of each species in a metacommunity
#'
#' The function counts the number of individuals of each species occuring
#' in a metacommunity and then divides by the total number of individuals.
#' If \code{only_species} is code{FALSE}, then empty spaces are included 
#' in the tally and the function returns the proportion of potential spaces
#' that each species occupies. Ranked abundances can be returned using the 
#' parameter \code{ranked}.
#'  
#' @param metacomm (required) matrix of lists defining metacommunity  
#' 	(as generated by \code{\link{populate_landscape}} or returned by
#' 	this function).
#' @param N_S (required) number of species expected in the metacommunity
#' @param only_species logical indicating whether the function should ignore
#' 	empty spaces. Defaults to \code{FALSE}
#' @param ranked logical indicating whether abundances should be returned in 
#' 	order of rank (\code{TRUE}) or in their original order (\code{FALSE}).
#' 	Default is \code{FALSE}.
#' @return a vector of length \code{N_S} or, if \code{only_species=FALSE}, 
#' 	a vector of length \code{N_S+1} where the first element corresponds to the number of
#' 	empty spaces 
#'
#' @export

calc_abun = function(metacomm, N_S, only_species=F, ranked=F){
	# Determine where to start species
	start = as.numeric(only_species)

	# Tally number of individuals of each species present, including empty spaces (0)
	abun = table(factor(unlist(metacomm), start:N_S))
	
	# Calculate relative abundance across all available spaces
	rel_abun = abun/sum(abun)

	# Rank abundances
	if(ranked) rel_abun = rel_abun[order(rel_abun, decreasing=T)]
	
	# Return abundances
	rel_abun
}