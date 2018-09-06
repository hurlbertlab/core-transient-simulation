#' Disperse individuals
#'
#' Disperses a set of individuals from a given cell and determines the locations 
#' of the new cells that each lands in. 
#'
#' For each dispersing individual, the function choses a random location 
#' within the cell of origin and then adds a dispersal vector, which is 
#' determined by the function \code{\link{get_dispersal_vec}}. The distance 
#' traveled is determined stochastically based on the dispersal kernel
#' defined in \code{form} and the species-specific dispersal rate defined in
#' the \code{d_rates} vector. Individuals are allowed to disperse off of the
#' landscape and the function may return cell locations outside of the 
#' intervals \code{[1,x] and [1,y]}.
#'
#' @param x (required) x coordinate of the origin cell or vector of length  
#' 2 with x and y coordinates
#' @param y y coordinate of the origin cell
#' @param propagules (required) vector of integers representing dispersing 
#' 	individuals. Values indicate which species each individual belongs to.
#' @param dim_land (required) vector with the dimensions of the landscape: 
#' 	(x, y)
#' @param  d_rates (required) vector of species dispersal rates from this cell
#' 	(based on its habitat type)
#' @param form list defining the dispersal kernel (see 
#' 	\code{\link{get_dispersal_vec}} for options). Defaults to Gaussian.
#' @return a matrix whose columns give the x and y coordinates of the cells 
#' where dispersing individuals end up. 
#'
#' @seealso \code{\link{get_dispersal_vec}} for calculation of dispersal vectors 
#' 	and for how to specify dispersal kernels.
#' @export

disperse = function(x, y=NA, propagules, dim_land, d_rates, form=list(type='gaussian')){

	# If only x specified, define y
	if(length(x)==2){	
		x = x[1]
		y = x[2]
	}

	# Define limits of origin cell
	origin_cell = rbind(c(x-1, x), c(y-1, y))	

	# For each propagule
	new_locs = sapply(propagules, function(i){

		# If the propagule has no species (i.e. and empty space)
		if(i == 0){

			destination_cell = c(NA, NA)

		# Otherwise determine dispersal based on species dispersal rate
		} else {

			# Get dispersal rate
			d = d_rates[i]	

			# Stochasticaly choose dispersal distance
			vec = get_dispersal_vec(d, form)

			# Stochastically choose an origin location from within origin cell
			origin_loc = apply(origin_cell, 1, function(cell) runif(1, cell[1], cell[2]))
	
			# Find cell that propagule lands in
			# Angles are counter-clockwise from east (as in mathematical notation)
			dx = vec[1]*sin(vec[2])
			dy = vec[1]*cos(vec[2])
	
			# Calculate destination and cell propagule lands in
			destination = origin_loc + c(dx, dy)
			destination_cell = ceiling(destination)	
		}

		destination_cell
	})

	# Return table of new locations
	t(new_locs)

}
