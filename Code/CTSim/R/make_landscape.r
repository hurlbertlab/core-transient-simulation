#' Creation of new landscape
#'
#' Creates a landscape as a grid with two habitat types represented by a raster
#' 	 of 1 and -1. 
#' 
#' @note Habitat heterogeneity increases as the ratio between \code{d} and grid 
#' 	dimensions decreases. WARNING: long compute time on grids > 50 x 50.
#' 
#' @param x (required) number of cells in the x dimension, or a vector of the x and y
#' 	dimensions
#' @param y number of cells in the y dimension
#' @param mod variogram model used to generate spatial autocorrelation 
#' 	(see  \code{\link{vgm}}).
#' 	Defaults to Exponential model with partial sill = 1.
#' @param d range of the variogram model. Defines the distance at which cells
#' 	are correlated. Defaults to 1/3 of the smallest dimension. 
#' @param prop proportion of cells assigned to have value = -1 (type 'A').
#' @param draw_plot logical indicating whether function should plot the
#' 	landscape
#' @return a raster object
#' @seealso \code{\link{vgm}} for how to specify variogram model
#'
#' @import gstat
#' @import sp
#' @import raster
#' @export

make_landscape = function(x=NULL, y=NA, mod=NULL, d=NA, prop=NA, draw_plot=F){

	# Catch error if no dimensions specified
	if(is.null(x)&is.na(y)){ stop('Must supply grid dimensions.') }
	if(is.na(y)&length(x)==2){
		y = x[2]
		x = x[1]
	} 
	if(is.na(y)) { stop('Must supply grid dimensions.') }
	
	# Define coordinates
	locs = expand.grid(1:x, 1:y)
	names(locs) = c('x','y')

	# Define spatial model with an expected value of 0 and no spatial trend (by universal kriging)
	# See vgm() in gstat package for explanation of spatial correlation models
	# Range defaults to 1/3 of the smallest grid dimension
	if(is.null(mod)){
		mod = gstat::vgm(psill=1, model='Exp', range=ifelse(is.na(d), min(x/3,y/3), d))
	}
	spatial_model = gstat::gstat(formula=z~1, locations=~x+y, dummy=T, beta=0, model=mod)

	# Define proportion habitat B if unspecified
	if(is.na(prop)) prop = 0.5

	# Simulate values at locations
	values = predict(spatial_model, locs, nsim=1)

	# Convert to spatial data and raster grid
	sp::gridded(values) = ~x+y
	values_grid = raster::raster(values)

	# Threshold to binary habitat based on proportion of habitat in each type
	threshold = raster::quantile(values_grid, prop)
	binary_grid = values_grid > threshold
	binary_grid[binary_grid==0] = -1

	# Draw landscape
	if(draw_plot) plot(binary_grid)

	# Return landscape
	binary_grid
}