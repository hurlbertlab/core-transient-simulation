#' Death
#'
#' Determines which individuals in a community die during one time step.
#'
#' \code{comm} is a vector (or list) of integers indicating which species 
#' are present. Individuals in the community die with species-specific 
#' probabilities provided in \code{m_rates}. The order of \code{m_rates}  
#' designates which species each mortality rate applies to. 
#' For example, \code{m_rates[1]} is the mortality rate of species 1.
#'
#' @param comm (required) vector (or list) of species present
#' 	(may contain zeros for empty spaces)
#' @param m_rates (required) vector of mortality rates for each species
#' 	in this cell (may be habitat specific)
#' @return a vector of integers of the same lenth as \code{code}
#' 	where 0s replace individuals that have died

die = function(comm, m_rates){
	
	# Make comm into a vector if not already
	comm = unlist(comm)

	# Probabilistically kill each individual according to mortality rates in this habitat
	survived = sapply(comm, function(sp){
		new_sp = if(sp==0){ 0 } else { 
			ifelse(runif(1) <= m_rates[sp], 0, sp)
		}
	})

	# Return survivors
	survived
}