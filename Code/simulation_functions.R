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
# 	distribution : named list of parameters describing the functional form of the SAD
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
	if(distribution$type=='same'){
		use_p = ifelse(is.null(distribution$p1), 1, distribution$p1)
		abuns = rep(use_p, N_S)
	}

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
# 	S_A : number of specialist species on habitat type A
# 	S_B : number of specialist species on habitat type B
# 	S_AB : number of generalist species (defaults to 0)
# 	dist_b : nameed list describing the distributions from which birth rates should be generated (see make_sad() function). Defaults to uniform on [1,10].
# 	m : vector of length 2 or 3 specifying death rates [preferred habitat, not preferred habitat, generalist rate]
# 	dist_d : named list of parameters describing distribution from which dispersal kernals should be drawn
#		mu = mean dispersal distance
#		var = variance of distribution from which kernals drawn. Defaults to 0 for same dispersal kernal for all species.
make_species = function(S_A=NULL, S_B=NULL, S_AB=0, dist_b = list(maxN=10, type='uniform'), m, dist_d=list(mu=1, var=0)){

	# Catch error if no species specified
	if(is.null(S_A)|is.null(S_B)) stop('Must specify number of species.')

	# Create array to hold species vital rates (b = birth, m = death, d = dispersal)
	# Rates are per timestep
	N_S = S_A + S_B + S_AB
	species_rates = array(NA, dim=c(N_S, 2, 3), dimnames=list(species=1:N_S, habitat=c('A','B'), rate=c('b','m','d')))

	# Specialist species birth rates in their preferred habitat ranged from 1 to maxN and are 0 in the unpreferred habitat
	species_rates[1:S_A, 'A', 'b'] = make_sad(S_A, dist_b)
	species_rates[1:S_A, 'B', 'b'] = 0
	species_rates[(S_A+1):(S_A+S_B), 'B', 'b'] = make_sad(S_B, dist_b)
	species_rates[(S_A+1):(S_A+S_B), 'A', 'b'] = 0

	# Generalist birth rates range from 1 to maxN in both habitats, but is not necessarily the same in both habitats
	if(S_AB > 0) species_rates[(N_S-S_AB+1):N_S, ,'b'] = make_sad(2*S_AB, dist_b)

	# Death rates are prespecified constants
	# This could be modified
	species_rates[1:S_A, 'A', 'm'] = m[1]
	species_rates[1:S_A, 'B', 'm'] = m[2]
	species_rates[(S_A+1):(S_A+S_B), 'B', 'm'] = m[1]
	species_rates[(S_A+1):(S_A+S_B), 'A', 'm'] = m[2]

	# Generalist death rates default to the death rate for the preferred habitat, but can be specified separately
	if(S_AB > 0) species_rates[(N_S-S_AB+1):N_S, ,'m'] = ifelse(length(m)==3, m[3], m[1])

	# Mean dispersal distances for each species are drawn from a gamma distribution with mean 'mu' and variance 'var'
	if(dist_d$var==0){
		species_rates[,,'d'] = dist_d$mu
	} else {
		theta = dist_d$var / dist_d$mu
		k = dist_d$mu / theta
		species_rates[,,'d'] = rgamma(N_S, shape=k, scale=theta)
	}

	# Return rates
	species_rates
}


# A function that creates initial communities across a landscape
#	land : binary matrix of habitat types, as generated by make_landscape()
#	species : array of species rate, as generated by make_species()
#	gsad : vector of global relative abundance distribution for species. Must be in the same order as species. Defaults to same abundance for all species.
#	K : either a single number defining a carrying capacity for all cells or a matrix of carrying capacities for each landscape cell. Must match dimensions of land.
#	distribution : character describing how individuals should be disributed across the landscape. Defaults to 'same'.
#		'same' (all cells recieve same proportion of carrying capacity)
#		'uniform' (individuals uniformly distributed across landscape)
#		'designated' (only certain cells receive individuals)
#	p : proportion of carrying capacity to be filled. Defaults to 1 = full.
#	which_cells : two column matrix indicated the x and y coordinates of cells to receive individuals. May optionally contain 3rd column with different proportions for each cell.

populate_landscape = function(land, species, gsad=NULL, K, distribution='same', p=1, which_cells=NA){
	
	# Dimensions of landscape
	X = nrow(land)
	Y = ncol(land)

	# Make carrying capacity matrix, if not specified
	if(length(K)==1){
		K_mat = matrix(K, nrow=X, ncol=Y)
	} else {
		K_mat = K
	}

	# Catch error if landscape and K_mat do not match
	# Note that land may be a raster layer so dim() will not match
	if(nrow(K_mat)!=X | ncol(K_mat)!=Y) stop('Dimensions of land and K_mat must match.')

	# Number of species
	N_S = dim(species)[1]

	# Define global species abundance distribution
	if(is.null(gsad)) gsad = rep(1, N_S)
	gsad = gsad/sum(gsad)

	# Make matrix to hold metacommunity. Integers indicate which species is present in each space.
	metacomm = matrix(list(), nrow=X, ncol=Y)
	for(i in 1:X){
	for(j in 1:Y){
		metacomm[i,j] = list(rep(0,K_mat[i,j]))
	}}

	# All cells receive the same proportion of their carrying capacity
	if(distribution=='same'){
		for(i in 1:X){
		for(j in 1:Y){
			k = length(unlist(metacomm[i,j]))
			n = floor(p*k)
			if(n>0) metacomm[i,j][[1]][1:n] = sample(N_S, size=n, replace=T, prob=gsad) 
		}}
	}

	# Metacommunity receives fixed proportion of total carrying capacity with individuals distributed uniformly
	if(distribution=='uniform'){
		k = sum(K_mat)
		n = floor(k*p)
		inds = sample(N_S, size=n, replace=T, prob=gsad)
		
		# Keep track of which cells have been filled
		cell_full = matrix(F, nrow=X, ncol=Y)	
	
		for(x in inds){
			# Find cells with spaces available
			available = which(!cell_full, arr.ind=T)
			this_cell = available[sample(nrow(available), 1),]

			# Find an available space in this cell
			these_spots = metacomm[this_cell[1], this_cell[2]][[1]]
			this_spot = which(these_spots==0)[1]

			# Assign this individual to this cell
			metacomm[this_cell[1], this_cell[2]][[1]][this_spot] = x

			# Record whether the cell is now full
			if(this_spot==length(these_spots)) cell_full[this_cell[1], this_cell[2]] = T		
		}

	}

	# Individuals placed in designated cells only.
	if(distribution=='designated'){
	
		# Catch error if locations not specified
		if(is.na(which_cells)) stop("Must specify which_cells if distribution is 'designated'.")
		
		# Assign individuals to each specified cell
		# Currently adds fixed proportion 'p' as default if exact number of individuals for each cell are unspecified.
		for(r in 1:nrow(which_cells)){
			i = which_cells[r,1]
			j = which_cells[r,2]
			k = length(unlist(metacomm[i,j]))
			n = ifelse(ncol(which_cells)==3, which_cells[i,3], floor(k*p))

			if(n > k) stop(paste0('More individuals specified than carrying capacity of cell (',i,',',j,')'))
				
			metacomm[i,j][[1]][1:n] = sample(N_S, size=n, replace=T, prob=gsad)
		}		
	}

	# Return initial metacommunity
	metacomm
}


# TESTING
mycomm = populate_landscape(land,species,K=matrix(rpois(400,40),ncol=20), p=.3, distribution='uniform')


####################################################################
### Functions for running a simulation




disperse = function(mu, sigma, form='gaussian'){

	gaussian

	NE # negative exponential

	IP # inverse power

	ENE # extended negative exponential

	FT # fat tail

}














####################################################################
### Functions for collecting data from a simulation