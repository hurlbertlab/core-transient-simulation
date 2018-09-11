#' Aggregate habitat types
#'
#' Determine the dominant habitat type in a set of spatial units
#'
#' The function determined which habitat type is present in the majority of
#' cells in each spatial unit specified in the list \code{locs}.
#' 
#' @param land (required) a raster landscape (as returned by 
#' 	\code{\link{make_landscape}}. 
#' @param locs (required) either a two-column matrix of cell locations defining
#' 	a single spatial unit or a list of matrices designating
#' 	multiple spatial units
#' @return a vector of habitat types corresponding to the spatial units
#'
#' @export

aggregate_hab_type = function(land, locs){

	# Convert matrix of cell locations to list if necessary.
	if(!is.list(locs)) locs = list(locs)
	
	# For each location, count the number of cells in each habitat
	hab_tab = sapply(locs, function(cell_block){
		
		# Get values and convert to habitat types
		hab_vals = land[cell_block]
		hab_types = sapply(hab_vals, get_habitat)
		
		# Table the types
		hab_tab = table(factor(hab_types, levels=c('A','B')))

		# Return counts and total
		c(hab_tab, n=sum(hab_tab))
	})

	# Calculate which habitat type dominates in each spatial unit
	habitats = apply(hab_tab, 2, function(x){ 
		dom = x[c('A','B')]/x['n'] > 0.5
		if(sum(dom)==0){
			sample(c('A','B'), 1)
		} else {
			names(which(dom))
		}
	})
	
	habitats
}
