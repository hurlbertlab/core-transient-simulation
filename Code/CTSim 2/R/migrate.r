#' Immigration
#'
#' Determine whether a species migrates into the metacommunity from outside 
#' the landscape.
#'
#' Species identity is determined probabilistically according to a
#' global species abundance distribution (\code{gsad}).
#'
#' @param m (required) probability that an individual will arrive from outside
#' 	the landscape.
#' @param gsad (required) global relative abundance of each species. Determines
#' probability that migrant will belong to a given species. 
#' @return integer representing the species that migrates. 0 = no migrant.
#'
#' @export

migrate = function(m, gsad){
	
	# Stochastically determine whether a migrant arrives
	arrived = runif(1) <= m

	# If a migrant arrived, choose its identity based on global abundance distribution
	migrant = ifelse(arrived, sample(length(gsad), 1, prob=gsad), 0)

	# Return the migrant
	migrant

}
