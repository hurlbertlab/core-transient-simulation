#' Determine habitat type
#'
#' Determines the habitat type that corresponds to a numeric value
#' (as stored in a landscape).
#' 
#' @note Does not work directly on vectors. Use \code{sapply(x, get_habitat)}.
#'
#' @param x (required) a number
#' 
#' @return 'A' or 'B'
#'
#' @export

get_habitat = function(x){
	h = NA
	
	if(x<=0) h = 'A'
	if(x>0) h = 'B'
	

	# Return habitat types
	h
}