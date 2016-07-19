#' Determine species type
#'
#' Determines whether each species in a species pool is a specialist on
#' habitat A or B, or a generalist using birth rates.
#'
#' User should supply an array \code{species} in which the third dimension
#' 	has a column named 'b' that stores birth rates (such as returned by 
#' 	\code{\link{make_species}}.) 
#'
#' @param species (required) array of species vital rates 
#' @return a vector of characters indicating which habitat each species 
#' 	prefers: 'A', 'B', or 'AB' (for a generalist).
#'
#' @seealso \code{\link{make_species}}
#' @export

get_sptype = function(species){
	
	# For each species determine in which habitat type the species birth rate is > 0	
	apply(species[,,'b'], 1, function(x) paste(names(which(x>0)), collapse=''))
}
