## This script contains functions for running the Core-Transient Simulation

## Required packages
library(gstat)
library(sp)
library(raster)
library(poweRlaw)




####################################################################
### Functions for initialing a simulation


# Function that creates a landscape as a grid with two habitat types represented by a raster of 1 and -1
# Heterogeneity increases as the ratio between d and grid dimensions decreases.
# WARNING: THIS FUNCTION HAS A LONG COMPUTE TIME FOR LARGE GRIDS
# 	x : number of cells in the x dimension, or a vector of the x and y dimensions
# 	y : number of cells in the y dimension
# 	mod : variogram model used to generate spatial autocorrelation (defined by vgm() function in gstat package). Defaults to Exponential model with partial sill = 1
# 	d : range of the variogram model. Defines the distance at which cells are correlated. Defaults to 1/3 of the smallest dimension. 
# 	prop : proportion of cells assigned to have value = 1.
# 	draw_plot : logical indicating whether function should plot the landscape
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






# A function that generates a species abundance distribution.
# Code is set up so that abundance may be correlated with another variable, but this is currently not implemented.
# To add a correlation, generate x as a multivariate random normal variable with specified covariance matrix.
# 	N_S : number of species
# 	distribution : list of parameters describing the functional form of the SAD
#		type = vector indicating the form of the SAD
#		maxN = maximum abundance, used to find SAD parameters if not supplied
#		P_maxN = probability of getting an abundance greater than maxN
#		p1 = parameter controling the shape of the SAD (only certain types)
#		p2 = parameter controling the shape of the SAD (only certain types)
#		corr = matrix giving correlations between global abundance and niche means and breadths (rows are env axes)
make_sad = function(N_S, distribution){
	
	# Generate normal random variable
	x = rnorm(N_S)
	U = pnorm(x)

	# Convert to correct distribution
	if(distribution$type=='same') abuns = rep(1, N_S)
	
	if(distribution$type=='uniform'){
		if(is.null(distribution$maxN)){ 
			use_p = distribution$p1
		} else {
			use_p = distribution$maxN
		}
		abuns = qunif(U, 1, use_p)
	}
		
	if(distribution$type=='power'){
		if(is.null(distribution$maxN)){
			use_p = distribution$p1
		} else {
			use_p = 1/log(1/(distribution$maxN), base=distribution$P_maxN)
		}
		abuns = (1-U)^(-1/use_p)
	}
	if(distribution$type=='logseries'){
		# p2 = N, p1 = Fisher's alpha
		if(is.null(distribution$maxN)){
			use_N = distribution$p2
		} else {
			use_N = 1
			this_N = qls(1-distribution$P_maxN, use_N, distribution$p1)	
			while(this_N < distribution$maxN){
				use_N = use_N + 1
				this_N = qls(1-distribution$P_maxN, use_N, distribution$p1)	
			}
			if(this_N > 10) use_N = use_N - 1

		}
		abuns = qls(U, N = use_N, alpha = distribution$p1)		
	}
	if(distribution$type=='lognormal'){
		if(is.null(distribution$maxN)){
			use_mean = distribution$p1
			use_sd = distribution$p2
		} else {
			# Assume mean = 0
			use_mean = 0
			try_sig = 1
			this_N = qlnorm(1-distribution$P_maxN, use_mean, try_sig)
			i = 1
			while(abs(distribution$maxN - this_N) > 0.5){
				if(this_N > distribution$maxN){
					try_sig = try_sig - try_sig/2
					
				} else {
					try_sig = try_sig + try_sig/2
				}
				this_N = qlnorm(1-distribution$P_maxN, 0, try_sig)
			}
			use_sd = try_sig
		}

		# p1 = mean, p2 = sd
		abuns = ceiling(qlnorm(U, use_mean, use_sd))
	}
	if(distribution$type=='poisson'){
		# p1 = lambda
		if(is.null(distribution$maxN)){
			use_lamda = distribution$p1
		} else {
			use_maxN = distribution$maxN + 1
			try_lamda = use_maxN:1
			this_N = qpois(1-distribution$P_maxN, try_lamda)
			N_diffs = this_N - use_maxN
			i = 1
			while(!(0 %in% N_diffs)){
				new_start = which(abs(N_diffs)==min(abs(N_diffs)))
				try_lamda = seq(try_lamda[new_start+1], try_lamda[new_start-1], 10^(-i))
				this_N = qpois(1-distribution$P_maxN, try_lamda)
				N_diffs = this_N - use_maxN
				i = i + 1
			}
		
			use_lamda = min(try_lamda[which(N_diffs==0)])
		}

		abuns = qpois(U, use_lamda)+1
	}

	# Discretize so that it can be used as an individual fecundity rate
	abuns = round(abuns, 0)

	abuns	
}



# Function that creates a dataframe defining species and matrix holding vital rates
# d : a vector of length 2 or 3 specifying specialist death rates in the preferred and not preferred habitat, and the generalist rate
make_species = function(S_A=NULL, S_B=NULL, S_AB=0, dist_b = list(maxN=10, type='uniform'), d){

	# Catch error if no species specified
	if(is.null(S_A)|is.null(S_B)) stop('Must specify number of species.')

	# Create array to hold species vital rates
	N_S = S_A + S_B + S_AB
	species_rates = array(NA, dim=c(N_S, 2, 2), dimnames=list(species=1:N_S, habitat=c('A','B'), rate=c('b','d')))

	# Specialist species birth rates in their preferred habitat ranged from 1 to maxN and are 0 in the unpreferred habitat
	species_rates[1:S_A, 'A', 'b'] = make_sad(S_A, dist_b)
	species_rates[1:S_A, 'B', 'b'] = 0
	species_rates[(S_A+1):(S_A+S_B), 'B', 'b'] = make_sad(S_B, dist_b)
	species_rates[(S_A+1):(S_A+S_B), 'A', 'b'] = 0

	# Generalist birth rates range from 1 to maxN in both habitats, but is not necessarily the same in both habitats
	if(S_AB > 0) species_rates[(N_S-S_AB+1):N_S, ,'b'] = make_sad(2*S_AB, dist_b)

	# Death rates are prespecified constants
	species_rates[1:S_A, 'A', 'd'] = d[1]
	species_rates[1:S_A, 'B', 'd'] = d[2]
	species_rates[(S_A+1):(S_A+S_B), 'B', 'd'] = d[1]
	species_rates[(S_A+1):(S_A+S_B), 'A', 'd'] = d[2]

	# Generalist death rates default to the death rate for the preferred habitat, but can be specified separately
	if(S_AB > 0) species_rates[(N_S-S_AB+1):N_S, ,'d'] = ifelse(length(d)==3, d[3], d[1])

	# TO DO: ADD IN DISPERSAL RATE FOR SPECIES AND ALLOW BIRTH AND DISPERSAL RATES TO BE DRAWN FROM CORRELATED DISTRIBUTIONS.




}
















####################################################################
### Functions for running a simulation


















####################################################################
### Functions for collecting data from a simulation