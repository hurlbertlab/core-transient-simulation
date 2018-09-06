#' Establishment
#'
#' Determines which individuals establish in the empty spaces in a community, 
#' given a pool of potenial recruits and probability that migrants will arrive
#' from outside the metacommunity. 
#'
#' For each empty space in the community \code{comm} (designated by a 0), 
#' the function first attempts to fill the space with a migrant from
#' outside the metacommunity with probability \code{m}. If successful, the species
#' identity is probabilistically chosen according to the relative abundances
#' provided in \code{gsad}. Spaces remaining empty after immigration are filled
#' stochastically from the pool of potential recruits provided in 
#' \code{propagules}. Individuals are chosen with equal probability.
#' Species recruitment rates define the probability that an individual 
#' selected for an empty space will actually establish.
#'
#' @param comm (required) vector (or list) of species present- may contain 
#' 	zeros for an empty spaces
#' @param propagules (required) vector of integers representing pool of 
#' 	potential recruits. Values indicate species identity of each individual. 
#' @param r_rates (required) vector of species' recruitment rates for this cell
#' 	(based on its habitat type). See details.
#' @param m immigration rate. Defaults to 0.
#' @param gsad vector defining the global relative abundance of each species.
#' 	Must be in the same order as \code{species}. Defaults to same abundance 
#'	for all species. See details.
#' @return a vector of the same length as \code{comm} defining the new community
#' 
#' @export

establish = function(comm, propagules, r_rates, m=0, gsad=NULL){
	
	# Make comm into a vector if not already
	comm = unlist(comm)

	# If m > 0, but gsad is unspecified, create gsad with equal abundance
	if(m > 0 & is.null(gsad)) gsad = make_sad(length(r_rates), distribution=list(type='same'))

	# For each empty space in the community
	recruits = comm[comm==0]

	# Select potential recruit from outside the landscape with probability m
	if(m>0){
		migrants = runif(length(recruits)) <= m
		recruits[migrants] = sample(length(gsad), size=sum(migrants), replace=T, prob=gsad)
	} else {
		migrants = rep(F, length(recruits))
	}
	
	# Select a potential recruit from the pool of propagules
	if(length(propagules)>1) recruits[!migrants] = sample(propagules)[1:sum(!migrants)]
	if(length(propagules)==1) recruits[!migrants] = propagules[1:sum(!migrants)]

	# Probabilistically determine whether propagules establish based on species recruitment rates for this habitat
	established = sapply(recruits, function(sp){
		ifelse(sp==0 | is.na(sp), 0, ifelse(runif(1) > r_rates[sp], 0, sp))
	})
	
	# Return new community
	comm[comm==0] = established
	comm
}
