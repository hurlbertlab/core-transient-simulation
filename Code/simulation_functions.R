## This script contains functions for running the Core-Transient Simulation

## Required packages
library(gstat)
library(sp)
library(raster)




####################################################################
### Functions for initialing a simulation


# Function that creates a landscape as a grid with two habitat types represented by a raster of 1 and -1
# Heterogeneity increases as the ratio between d and grid dimensions decreases.
# WARNING: THIS FUNCTION HAS A LONG COMPUTE TIME FOR LARGE GRIDS
# x : number of cells in the x dimension, or a vector of the x and y dimensions
# y : number of cells in the y dimension
# mod : variogram model used to generate spatial autocorrelation (defined by vgm() function in gstat package). Defaults to Exponential model with partial sill = 1
# d : range of the variogram model. Defines the distance at which cells are correlated. Defaults to 1/3 of the smallest dimension. 
# prop : proportion of cells assigned to have value = 1.
# draw_plot : logical indicating whether function should plot the landscape
make_landscape = function(x=NULL, y=NULL, mod=NULL, d=NULL, prop=0.5, draw_plot=F){
	
	# Catch error if no dimensions specified
	if(is.null(x)&is.null(y)){ stop('Must supply grid dimensions.') }
	if(is.null(y)&length(x)==2){
		y = x[2]
		x = x[1]
	} 
	if(is.null(y)) { stop('Must supply grid dimensions.') }
	
	# Define coordinates
	locs = expand.grid(1:x, 1:y)
	names(locs) = c('x','y')

	# Define spatial model with an expected value of 0 and no spatial trend (by universal kriging)
	# See vgm() in gstat package for explanation of spatial correlation models
	# Range defaults to 1/3 of the smallest grid dimension
	if(is.null(mod)){
		mod = vgm(psill=1, model='Exp', range=ifelse(is.null(d), min(x/3,y/3),d))
	}
	spatial_model = gstat(formula=z~1, locations=~x+y, dummy=T, beta=0, model=mod)

	# Simulate values at locations
	values = predict(spatial_model, locs, nsim=1)

	# Convert to spatial data and raster grid
	gridded(values) = ~x+y
	values_grid = raster(values)

	# Threshold to binary habitat based on proportion of habitat in each type
	threshold = quantile(values_grid, prop)
	binary_grid = values_grid > threshold
	binary_grid[binary_grid==0] = -1

	# Draw landscape
	if(draw_plot) plot(binary_grid)

	# Return landscape
	binary_grid
}

# Function that calculates degree of habitat aggregation






















####################################################################
### Functions for running a simulation


















####################################################################
### Functions for collecting data from a simulation