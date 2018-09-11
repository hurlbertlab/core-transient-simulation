#' Average habitat values
#'
#' Calculate the habitat type that corresponds to the average
#' habitat value across a set of cells in a landscape.
#'
#' @param locs (required) a 2-column matrix of cell coordinates
#' @param land (required) a raster defining the landscape 
#' 	(see \code{\link{make_landscape}})
#' @return \code{'A'} or \code{'B'}
#'
#' @export

average_habitat = function(locs, land){

	# Get mean value across all cells
	vals  = apply(locs, 1, function(i) land[i[1], i[2]])

	# Determine habitat type
	# Note: equal parts 'A' and 'B' gets classified as 'A'
	get_habitat(mean(vals))
}