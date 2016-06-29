## This script contains functions for running the Core-Transient Simulation

## Required packages
library(gstat, lib.loc='/nas02/home/j/r/jrcoyle/Rlibs/')
library(sp)
library(raster)
library(poweRlaw)
library(fdrtool) # half-normal distribution for gaussian dispersal kernal
library(sads) # qls()

### TO DO ###


### Conventions ###
# NA used for missing value with length 1
# NULL used for missing lists or other objects larger than length 1

####################################################################
### Functions for initiating a simulation


# Function that creates a landscape as a grid with two habitat types represented by a raster of 1 and -1
# Heterogeneity increases as the ratio between d and grid dimensions decreases.
# WARNING: THIS FUNCTION HAS A LONG COMPUTE TIME FOR LARGE GRIDS
# 	x : number of cells in the x dimension, or a vector of the x and y dimensions
# 	y : number of cells in the y dimension
# 	mod : variogram model used to generate spatial autocorrelation (defined by vgm() function in gstat package). Defaults to Exponential model with partial sill = 1
# 	d : range of the variogram model. Defines the distance at which cells are correlated. Defaults to 1/3 of the smallest dimension. 
# 	prop : proportion of cells assigned to have value = 1.
# 	draw_plot : logical indicating whether function should plot the landscape
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
		mod = vgm(psill=1, model='Exp', range=ifelse(is.na(d), min(x/3,y/3), d))
	}
	spatial_model = gstat(formula=z~1, locations=~x+y, dummy=T, beta=0, model=mod)

	# Define proportion habitat A if unspecified
	if(is.na(prop)) prop = 0.5

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
# ADD THIS LATER


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
# 	dist_b : named list describing the distributions from which birth rates should be generated (see make_sad() function). Defaults to uniform on [1,10].
# 	m : vector of length 2 or 3 specifying death rates [preferred habitat, not preferred habitat, generalist rate]
# 	r : vector of length 2 or 3 specifying recruitment rates [preferred habitat, not preferred habitat, generalist rate]
# 	dist_d : named list of parameters describing distribution from which propagule dispersal kernals should be drawn
#		mu = mean dispersal distance
#		var = variance of distribution from which kernals drawn. Defaults to 0 for same dispersal kernal for all species.
#	dist_v : named list of parameters describing distribution from which movement kernals should be drawn
#		mu = vector of length 2 or 3 specifying mean movement distance [preferred habitat, not preferred habitat, generalist rate]
#		var = vector of length 2 or 3 specifying variance in mean movement distance [preferred habitat, not preferred habitat, generalist rate]. Defaults to 0 (same dispersal kernal for all species).
make_species = function(S_A=NA, S_B=NA, S_AB=NA, dist_b = NULL, m, r, dist_d=NULL, dist_v=NULL){

	# Catch error if no species specified
	if(is.na(S_A)|is.na(S_B)) stop('Must specify number of species.')

	# Set generalists to 0 if undefined
	if(is.na(S_AB)) S_AB = 0

	# Set distributions to defaults if undefined
	if(is.null(dist_b)) dist_b = list(maxN=10, type='uniform')
	if(is.null(dist_d)) dist_d = list(mu=1, var=0)
	if(is.null(dist_v)) dist_v = list(mu=rep(0,3), var=rep(0,3))

	# Create array to hold species vital rates (b = birth, m = death, r = recruitment, d = dispersal, v = movement)
	# Rates are per timestep
	N_S = S_A + S_B + S_AB
	species_rates = array(NA, dim=c(N_S, 2, 5), dimnames=list(species=1:N_S, habitat=c('A','B'), rate=c('b','m','r','d','v')))

	# Specialist species birth rates in their preferred habitat ranged from 1 to maxN and are 0 in the unpreferred habitat
	species_rates[1:S_A, 'A', 'b'] = make_sad(S_A, dist_b)
	species_rates[1:S_A, 'B', 'b'] = 0
	species_rates[(S_A+1):(S_A+S_B), 'B', 'b'] = make_sad(S_B, dist_b)
	species_rates[(S_A+1):(S_A+S_B), 'A', 'b'] = 0

	# Death rates are prespecified constants
	# This could be modified to be drawn from a distribution
	species_rates[1:S_A, 'A', 'm'] = m[1]
	species_rates[1:S_A, 'B', 'm'] = m[2]
	species_rates[(S_A+1):(S_A+S_B), 'B', 'm'] = m[1]
	species_rates[(S_A+1):(S_A+S_B), 'A', 'm'] = m[2]

	# Specialist species recruitment rates are prespecified constants
	# This could be modified to be drawn from a distribution
	species_rates[1:S_A, 'A', 'r'] = r[1]
	species_rates[1:S_A, 'B', 'r'] = r[2]
	species_rates[(S_A+1):(S_A+S_B), 'B', 'r'] = r[1]
	species_rates[(S_A+1):(S_A+S_B), 'A', 'r'] = r[2]

	# Mean dispersal distances for each species are drawn from a gamma distribution with mean 'mu' and variance 'var'
	if(dist_d$var==0){
		species_rates[,,'d'] = dist_d$mu
	} else {
		theta = dist_d$var / dist_d$mu
		k = dist_d$mu / theta
		species_rates[,,'d'] = rgamma(N_S, shape=k, scale=theta)
	}

	# Mean movement distances for each species are drawn from a gamma distribution with mean 'mu' and variance 'var'.
	# Preferred habitat
	if(dist_v$var[1]==0){
		species_rates[1:S_A, 'A', 'v'] = dist_v$mu[1]
		species_rates[(S_A+1):(S_A+S_B), 'B', 'v'] = dist_v$mu[1]
	} else {
		theta = dist_v$var[1] / dist_v$mu[1]
		k = dist_v$mu[1] / theta
		species_rates[1:S_A, 'A', 'v'] = rgamma(S_A, shape=k, scale=theta)
		species_rates[(S_A+1):(S_A+S_B), 'B', 'v'] = rgamma(S_B, shape=k, scale=theta)
	}
	# Unpreferred habitat
	if(dist_v$var[2]==0){
		species_rates[1:S_A, 'B', 'v'] = dist_v$mu[2]
		species_rates[(S_A+1):(S_A+S_B), 'A', 'v'] = dist_v$mu[2]
	} else {
		theta = dist_v$var[2] / dist_v$mu[2]
		k = dist_v$mu[2] / theta
		species_rates[1:S_A, 'B', 'v'] = rgamma(S_A, shape=k, scale=theta)
		species_rates[(S_A+1):(S_A+S_B), 'A', 'v'] = rgamma(S_B, shape=k, scale=theta)
	}
	

	# Assign rates for generalists, if they exist
	if(S_AB > 0){
	
		# Generalist birth rates range from 1 to maxN in both habitats, but is not necessarily the same in both habitats
		species_rates[(N_S-S_AB+1):N_S, ,'b'] = make_sad(2*S_AB, dist_b)
	
		# Generalist death rates default to the death rate for the preferred habitat, but can be specified separately
		species_rates[(N_S-S_AB+1):N_S, ,'m'] = ifelse(length(m)==3, m[3], m[1])

		# Generalist recruitment rates default to the recruitment rate for the preferred habitat, but can be specified separately
		species_rates[(N_S-S_AB+1):N_S, ,'r'] = ifelse(length(r)==3, r[3], r[1])

		# Generalist movement rates are the same in both habitat types and default to the movement rate in the preferred habitat, but can be specificed separately
		if(length(dist_v$mu)==2) dist_v$mu = c(dist_v$mu, dist_v$mu[1])
		if(length(dist_v$var)==2) dist_v$var = c(dist_v$var, dist_v$var[1])

		if(dist_v$var[3]==0){
			species_rates[(N_S-S_AB+1):N_S, , 'v'] = dist_v$mu[3]
		} else {
			theta = dist_v$var[3] / dist_v$mu[3]
			k = dist_v$mu[3] / theta
			species_rates[(N_S-S_AB+1):N_S, ,'v'] = rgamma(S_AB, shape=k, scale=theta)
		}
	}

	# Return rates
	species_rates
}

# A function that returns a vector saying which type of habitat species are residents of.
# Determination based on species birth rates being greater than 0.
#	species : array of species vital rates.
get_sptype = function(species){
	
	# For each species determine in which habitat type the species birth rate is > 0	
	apply(species[,,'b'], 1, function(x) paste(names(which(x>0)), collapse=''))
}

# A function that creates initial communities across a landscape
#	land : binary matrix of habitat types, as generated by make_landscape()
#	species : array of species rate, as generated by make_species()
#	gsad : vector of global relative abundance distribution for species. Must be in the same order as species. Defaults to same abundance for all species.
#	K : either a single number defining a carrying capacity for all cells or a matrix of carrying capacities for each landscape cell. Must match dimensions of land.
#	distribution : character describing how individuals should be disributed across the landscape. Defaults to 'same' unless which_cells are specified.
#		'same' (all cells recieve same proportion of carrying capacity)
#		'uniform' (individuals uniformly distributed across landscape)
#		'designated' (only certain cells receive individuals)
#	p : proportion of carrying capacity to be filled. Defaults to 1 = full.
#	which_cells : two column matrix indicated the x and y coordinates of cells to receive individuals. May optionally contain 3rd column with different proportions for each cell.

populate_landscape = function(land, species, gsad=NULL, K, distribution=NA, p=NA, which_cells=NULL){
	
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

	# Define parameters to default values is unspecified
	if(is.na(distribution)&is.null(which_cells)) distribution = 'same'
	if(is.na(distribution)&!is.null(which_cells)) distribution = 'designated'
	if(is.na(p)) p = 1


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
		if(is.null(which_cells)) stop("Must specify which_cells if distribution is 'designated'.")
		
		# Assign individuals to each specified cell
		# Currently adds fixed proportion 'p' as default if exact number of individuals for each cell are unspecified.
		for(r in 1:nrow(which_cells)){
			i = which_cells[r,1]
			j = which_cells[r,2]
			k = length(unlist(metacomm[i,j]))
			n = ifelse(ncol(which_cells)==3, which_cells[r,3], floor(k*p)) # THROWS AN ERROR

			if(n > k) stop(paste0('More individuals specified than carrying capacity of cell (',i,',',j,')'))
				
			metacomm[i,j][[1]][1:n] = sample(N_S, size=n, replace=T, prob=gsad)
		}		
	}

	# Return initial metacommunity
	metacomm
}

####################################################################
### Functions for running a simulation

# Function that determines habitat type based on habitat value stored in landscape
# 	x : numeric habitat value stored in landscape
get_habitat = function(x){
	h = NA
	
	if(x<=0) h = 'A'
	if(x>0) h = 'B'
	

	# Return habitat types
	h
}


# Function that returns the pool of new propagules produce in a cell in a single time step
#	comm : vector (or list) of species present, may contain zeros for an empty space.
#	b_rates : vector of birth rates for species.
reproduce = function(comm, b_rates){

	# Convert to vetor if not already
	comm = unlist(comm)
	
	# Make pool of propagules
	propagules = sapply(comm[comm>0], function(sp) rep(sp, b_rates[sp]))
	propagules = unlist(propagules)
	
	propagules
}

# Function that probabilistically adds new immgrants to cell from global pool.
# 	m : probability that an individual will migrate in from the outside the landscape
#	gsad : relative abundance of each species globally. Could also represent long-distance dispersal ability. 
migrate = function(m, gsad){
	
	# Stochastically determine whether a migrant arrives
	arrived = runif(1) <= m

	# If a migrant arrived, choose its identity based on global abundance distribution
	migrant = ifelse(arrived, sample(length(gsad), 1, prob=gsad), 0)

	# Return the migrant
	migrant

}

# Function that disperses individual propagules from a given cell and returns a new cell location for each individual propagule.
# 	x : x coordinate of cell, or vector of length 2 with x and y coordinates
#	y : y coordinate of cell
#	propagules : vector of propagules of different species (denoted by numeric value).
#	dim_land : dimension fo landscape (x,y)
#	d_rates : vector of species dispersal rates from this cell, in same order as species. Dispersal rates are the mean distance for the given dispersal kernal form.
#	form : list describing the form of the dispersal kernel. See get_dispersal_dist() function 
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


# Function that stochastically generates a dispersal vector (distance, direction) for given a dispersal kernal
#	d : expected distance (mean)
#	form : a list defining the type of dispersal kernel. Currently only implements 'gaussian'.
#		type = character indicating the type
#	N : number of vectors to generate. Defaults to 1
get_dispersal_vec = function(d, form=list(type='gaussian'), N=1){
	
	# Gaussian kernel
	if(form$type=='gaussian'){
	
		# Calculate theta for mean = d
		theta = 1/d

		# Generate distance
		r = rhalfnorm(N, theta)	

		# Generate random angle (in radians)
		phi = runif(N, 0, 2*pi)	

	}

	# Adjacent cell dispersal
	# d is the probability that a propagule will leave its origin cell
	# form$moves = maximum number of moves that a propagule will complete in one timestep
	if(form$type=='adjacent'){

		# Generate random set of directional moves
		moves = sapply(1:N, function(n){
			sample(0:1, size=form$moves, prob= c(1-d, d), replace=T) * sample(1:4, size=form$moves, replace=T)
		})

		# Calculate distance and angle
		dX = apply(moves, 2, function(x) -1*sum(x==3) + 1*sum(x==1))
		dY = apply(moves, 2, function(x) -1*sum(x==4) + 1*sum(x==2))
		r = sqrt(dX^2 + dY^2)
		phi = atan2(dY,dX)
	}

	# Uniform: implements no disersal limitation
	# d is the maximum dispersal distance
	if(form$type=='uniform'){
	
		# Generate distance
		r = runif(N, 0, d)

		# Generate random angle		
		phi = runif(N, 0, 2*pi)	
	}


	# Return vector
	cbind(r, phi)
}


# A function that determines which propagules establish in a cell and returns a new community vector
# Includes potential recruitment from immigrants outside the landscape.
#	comm : vector (or list) of species present, may contain zeros for an empty space.
#	propagules : vector of propagules of different species (denoted by numeric value).
#	r_rates : vector of recruitment rates for each species in this habitat
# 	m : immigration rate. Probability that empty space colonized by propagule from outside the landscape. Defaults to 0.
#	gsad : global species abudance distribution. Defines probability that immigant from outside landscape will belong to a given species. Defaults to equal probability when unspecified.
establish = function(comm, propagules, r_rates, m=0, gsad=NULL){
	
	# Make comm into a vector if not already
	comm = unlist(comm)

	# If m > 0, but gsad is unspecified, create gsad with equal abundance
	if(m > 0 & is.null(gsad)) gsad = make_sad(length(r_rates), distribution=list(type='same'))

	# For each empty space in the community
	recruits = comm[comm==0]

	# Select potential recruit from outside the landscape with probability m
	if(m>0){
		migrants = runif(length(recruits)) <= m
		recruits[migrants] = sample(length(gsad), size=sum(migrants), replace=T, prob=gsad)
	} else {
		migrants = rep(F, length(recruits))
	}
	
	# Select a potential recruit from the pool of propagules
	if(length(propagules)>1) recruits[!migrants] = sample(propagules)[1:sum(!migrants)]
	if(length(propagules)==1) recruits[!migrants] = propagules[1:sum(!migrants)]

	# Probabilistically determine whether propagules establish based on species recruitment rates for this habitat
	established = sapply(recruits, function(sp){
		ifelse(sp==0 | is.na(sp), 0, ifelse(runif(1) > r_rates[sp], 0, sp))
	})
	
	# Return new community
	comm[comm==0] = established
	comm
}

# A function that determines which individuals die in a community
#	comm : vector (or list) of species present, may contain zeros for an empty space.
#	m_rates : vector of mortality rates (per time step) for each species in this habitat 
die = function(comm, m_rates){
	
	# Make comm into a vector if not already
	comm = unlist(comm)

	# Probabilistically kill each individual with according to death rates in this habitat
	survived = sapply(comm, function(sp){
		new_sp = if(sp==0){ 0 } else { 
			ifelse(runif(1) <= m_rates[sp], 0, sp)
		}
	})

	# Return survivors
	survived
}

# A function that runs a metacommunity on a landscape with a given species pool for one time step
# 	metacomm : a matrix of lists representing the community of individuals in each cell
#	land : a matrix of habitat values
#	species : an array containing species vital rates [species number, habitat type, rate type]
#	gsad : vector giving the global relative abundance of species in same order as species vital rates array
#	d_kernel : list defining the shape of the dispersal kernel for an individual (see disperse() function). Defaults to Gaussian.
#	imm_rate : probability that an empty space will be colonized by a propagule from outside the landscape (e.g. drawn from the global species abundance distribution). Defaults to 0.
run_timestep = function(metacomm, land, species, gsad, d_kernel=NULL, v_kernel=NULL, imm_rate=NA){

	# Define dimensions of landscape
	X = nrow(land)
	Y = ncol(land)

	# Catch error where metacomm and land are different dimensions
	if(nrow(metacomm)!=X | ncol(metacomm)!=Y) stop(paste0('Dimensions of metacommunity (',paste(dim(metacomm),collapse='x'),') must match dimensions of landscape (',X,'x',Y,').'))

	# Define individual dispersal parameters if unspecified
	if(is.null(d_kernel)) d_kernel = list(type='gaussian')
	if(is.null(v_kernel)) v_kernel = list(type='gaussian')
	if(is.na(imm_rate)) imm_rate = 0

	# Define array to hold new propagule pools
	propagule_pools = matrix(list(), nrow=X, ncol=Y)

	# For each cell, new individuals are born and disperse and existing individuals move around
	for(i in 1:X){
	for(j in 1:Y){

		# Define this community and habitat
		this_comm = metacomm[i,j][[1]]
		this_habitat = get_habitat(land[i,j])

		# Stochasitic mortality
		new_comm = die(this_comm, species[, this_habitat, 'm'])

		# New individuals born
		this_pool = reproduce(new_comm, species[,this_habitat,'b'])
	
		# Propagules disperse to new pools
		prop_locs = disperse(i, j, this_pool, c(X,Y), species[,this_habitat,'d'], form=d_kernel)
		if(length(prop_locs)>0){
			for(p in 1:nrow(prop_locs)){
				# Propagules that disperse off of the landscape are removed
				if(prop_locs[p,1]>0 & prop_locs[p,2]>0 & prop_locs[p,1]<=X & prop_locs[p,2]<=Y ) propagule_pools[prop_locs[p,1], prop_locs[p,2]] = list(unlist(c(propagule_pools[prop_locs[p,1], prop_locs[p,2]], this_pool[p])))
			}
		}

		# Previously established individuals have the opportunity to move 
		move_locs = disperse(i, j, new_comm, c(X,Y), species[,this_habitat,'v'], form=v_kernel)
		for(p in 1:length(new_comm)){
			
			# If the individual exists (i.e. this is not an empty place)
			if(new_comm[p]>0){	

			# If the individual moves to a different cell
			if(move_locs[p,1]!=i | move_locs[p,2]!=j){		

			# If if individual doesn't move off the landscape
			if(move_locs[p,1]>0 & move_locs[p,2]>0 & move_locs[p,1]<=X & move_locs[p,2]<=Y ){
				
				# Add individual to pool of propagules in the new cell
				propagule_pools[move_locs[p,1],move_locs[p,2]] = list(unlist(c(propagule_pools[move_locs[p,1],move_locs[p,2]], new_comm[p])))

				# Remove the individual from the current cell
				new_comm[p] = 0
		
			}}} # closes all three if clauses
		}

		# Save existing metacommunity
		metacomm[i,j] = list(as.numeric(new_comm))

	}}
	
	# Recruits establish from propagule pools (combination of new births and previously established individuals that moved)
	for(i in 1:X){
	for(j in 1:Y){

		# Define this community, habitat and pool of propagules
		this_comm = metacomm[i,j][[1]]
		this_habitat = get_habitat(land[i,j])
		this_pool = propagule_pools[i,j][[1]]
	
		# Determine which propagules estabish
		new_comm = establish(this_comm, this_pool, species[, this_habitat, 'r'], imm_rate, gsad)

		# Save results
		metacomm[i,j] = list(as.numeric(new_comm))
	}}

	# Return new metacommunity
	metacomm
}


# A function that runs a simulation for a given number of time steps.
# Returns an array of lists of individuals in each cell at each indicated time step.
#	steps : number of time steps to run the simulation.
# 	metacomm : a matrix of lists representing the community of individuals in each cell
#	land : a matrix of habitat values
#	species : an array containing species vital rates [species number, habitat type, rate type]
#	gsad : vector giving the global relative abundance of species in same order as species vital rates array
#	d_kernel = list defining the shape of the dispersal kernel for a new propagule (see disperse() function). Defaults to Gaussian.
#	v_kernel = list defining the shape of the movement kernel for an established individual (see disperse() function). Defaults to Gaussian.
#	imm_rate = probability that an empty space will be colonized by a propagule from outside the landscape (e.g. drawn from the global species abundance distribution).
#	save_steps : vector of timepoints at which to save the simulation. Defaults to all time steps.
# 	report : Number of time steps after which to report status. Defaults to 0, no reporting.
# 	ID : string defining this simulation run. Used when reporting progress. Defaults to NA.
run_sim = function(steps, metacomm, land, species, gsad, d_kernel=NULL, v_kernel=NULL, imm_rate=NA, save_steps = NULL, report=0, ID=NA){
	
	# Define simulation
	X = nrow(land)
	Y = ncol(land)

	# Create array to save simulations results
	if(is.null(save_steps)) save_steps = 0:steps
	sim_results = array(list(), dim=c(X, Y, length(save_steps)),
		dimnames=list(row=1:X, col=1:X, time=save_steps))
	
	# Save initial step
	if(0 %in% save_steps) sim_results[,,'0'] = metacomm

	# Run simulation
	new_metacomm = metacomm
	for(step in 1:steps){
		# Run for a timestep
		new_metacomm = run_timestep(new_metacomm, land, species, gsad, d_kernel, v_kernel, imm_rate)

		# Save results
		if(step %in% save_steps) sim_results[,,as.character(step)] = new_metacomm
		
		# Report progress
		if(report > 0) if(step %% report == 0) print(paste(Sys.time(), ': Finished', step, 'of', steps, 'from run', ID))
	}

	# Return results
	sim_results
}


# Function that runs a given number of replicates of a simulation under a given set of parameters
#	nruns : number of replicate simulations	
#	parms : list of parameters for running simulation, as would be created by read_parmfile() function
#	nparallel : number of cores to run the simulations on in parallel. Defaults to 1.
#	simID : string identifying this simulation run. Defaults to 'test'.
#	save_sim : file name under which to save simulation results. If none specified simulation results are returned but not saved.
#	sim_dir : directory where scripts with simulation functions are stored. Defaults to current directory.
# 	report : Number of time steps after which to report status. Defaults to 0, no reporting.
run_sim_N = function(nruns, parms, nparallel=1, simID='test', save_sim=NULL, sim_dir='./', report=0){

	if(nparallel > 1){
		require(parallel)
		cluster = makeCluster(nparallel, outfile=paste0(simID, '.Rout'))	
		
		# Send required functions and objects to each node
		clusterExport(cluster, c('parms','simID','save_sim','sim_dir','report'), envir=environment())
		clusterEvalQ(cluster, source(paste0(sim_dir, 'simulation_functions.R')))
		
		# Initialize simulation landscapes
		lands_N = parLapply(cluster, 1:nruns, function(j){
			with(parms, {
				x = dimX
				y = dimY
				if(!exists('vgm_mod')) vgm_mod = NULL
				d = ifelse(exists('vgm_dcorr'), vgm_dcorr, NA)
				prop = ifelse(exists('habA_prop'), habA_prop, 0.5)
				make_landscape(x, y, vgm_mod, d, prop, draw_plot=F)
			})
		})
	
		if(report>0) print(paste0(Sys.time(), ': Finished making landscapes.'))

		# Initialize species vital rates
		species_N = parLapply(cluster, 1:nruns, function(j){
			with(parms, {
				S_AB = ifelse(exists('S_AB'), S_AB, NA) 
				if(!exists('dist_b')) dist_b = NULL
				m = m_rates
				r = r_rates 
				if(!exists('dist_d')) dist_d = NULL
				if(!exists('dist_v')) dist_v = NULL
				make_species(S_A, S_B, S_AB, dist_b, m, r, dist_d, dist_v)
			})
		})

		if(report>0) print(paste0(Sys.time(), ': Finished making species pools.'))

		# Send initial landscapes and species to all cluster cores
		clusterExport(cluster, c('lands_N','species_N'), envir=environment())
		
		# Initialize global species abundance distribution
		# dist_gsad can be a generic distribution or 'b_rates' indicating that it should be the same as species birth rates
		gsad_N = parLapply(cluster, 1:nruns, function(j){
			with(parms, {
				N_S = dim(species_N[[j]])[1]
				if(exists('dist_gsad')){
					if(is.list(dist_gsad)){
						
						# Use specified distribution to generate abundances
	 					distribution = dist_gsad
						gsad_vec = make_sad(N_S, distribution)

					} else {

						# Make global abundaces equal to species birth rates	
						if(dist_gsad=='b_rates'){
			
							A_rates = species_N[[j]][1:S_A,'A','b']
							B_rates = species_N[[j]][(S_A+1):(S_A+S_B),'B','b']
							gsad_vec = c(A_rates, B_rates)
							if(exists('S_AB')) if(S_AB > 0) gsad_vec = c(gsad_vec, rowMeans(species_N[[j]][(S_A+S_B+1):(S_A+S_B+S_AB),,'b']))
					
						} else {
							stop('Unrecognized value for parameter dist_gsad.')
						}
					}

				# Defaults to same abundance for each species
				} else {	
					distribution = list(type='same')
					gsad_vec = make_sad(N_S, distribution)
				}
				
				# Return vector of global abundances
				gsad_vec
			})
		})

		if(report>0) print(paste0(Sys.time(), ': Finished making gsads.'))

		# Send global species abundance distributions to cluster
		clusterExport(cluster, 'gsad_N', envir=environment())

		# Run simulations
		sim_results = parLapply(cluster, 1:nruns, function(j){

			if(report>0) print(paste0(Sys.time(), ': Start run ', j))

			
			# Define this landscape and species pool
			this_land = lands_N[[j]]
			this_species = species_N[[j]]
			this_gsad = gsad_N[[j]]
			
			# Distribute species across landscape
			this_metacomm = with(parms, {
				p = ifelse(exists('prop_full'), prop_full, NA)
				distribution = ifelse(exists('init_distribute'), init_distribute, NA)
				if(exists('cells_distribute')){
					which_cells = cells_distribute
				} else {
					which_cells = NULL
				}
				populate_landscape(this_land, this_species, this_gsad, K, distribution, p, which_cells)
			})
			
			# Run simulation
			with(parms, {
				if(!exists('d_kernel')) d_kernel = NULL
				if(!exists('v_kernel')) v_kernel = NULL
				imm_rate = ifelse(exists('imm_rate'), imm_rate, NA)
				if(!exists('save_steps')) save_steps = NULL
				run_sim(nsteps, this_metacomm, this_land, this_species, this_gsad, d_kernel, v_kernel, imm_rate, save_steps, report, ID=j)
			})
		})

		stopCluster(cluster)	

	} else {

		# Initialize simulation landscapes
		lands_N = lapply(1:nruns, function(j){
			with(parms, {
				x = dimX
				y = dimY
				if(!exists('vgm_mod')) vgm_mod = NULL
				d = ifelse(exists('vgm_dcorr'), vgm_dcorr, NA)
				prop = ifelse(exists('habA_prop'), habA_prop, 0.5)
				make_landscape(x, y, vgm_mod, d, prop, draw_plot=F)
			})
		})

		if(report>0) print(paste0(Sys.time(), ': Finished making landscapes.'))

		# Initialize species vital rates
		species_N = lapply(1:nruns, function(j){
			with(parms, {
				S_AB = ifelse(exists('S_AB'), S_AB, NA) 
				if(!exists('dist_b')) dist_b = NULL
				m = m_rates
				r = r_rates 
				if(!exists('dist_d')) dist_d = NULL
				if(!exists('dist_v')) dist_v = NULL
				make_species(S_A, S_B, S_AB, dist_b, m, r, dist_d, dist_v)
			})
		})

		if(report>0) print(paste0(Sys.time(), ': Finished making species pools.'))

		# Initialize global species abundance distribution	
		gsad_N = lapply(1:nruns, function(j){
			with(parms, {
				N_S = dim(species_N[[j]])[1]
				if(exists('dist_gsad')){
					if(is.list(dist_gsad)){
						
						# Use specified distribution to generate abundances
	 					distribution = dist_gsad
						gsad_vec = make_sad(N_S, distribution)

					} else {

						# Make global abundaces equal to species birth rates	
						if(dist_gsad=='b_rates'){
			
							A_rates = species_N[[j]][1:S_A,'A','b']
							B_rates = species_N[[j]][(S_A+1):(S_A+S_B),'B','b']
							gsad_vec = c(A_rates, B_rates)
							if(exists('S_AB')) if(S_AB > 0) gsad_vec = c(gsad_vec, rowMeans(species_N[[j]][(S_A+S_B+1):(S_A+S_B+S_AB),,'b']))
					
						} else {
							stop('Unrecognized value for parameter dist_gsad.')
						}
					}

				# Defaults to same abundance for each species
				} else {	
					distribution = list(type='same')
					gsad_vec = make_sad(N_S, distribution)
				}
				
				# Return vector of global abundances
				gsad_vec
			})
		})

		if(report>0) print(paste0(Sys.time(), ': Finished making gsads.'))

		# Run simulations
		sim_results = lapply(1:nruns, function(j){
			
			if(report>0) print(paste0(Sys.time(), ': Start run ', j))

			# Define this landscape and species pool
			this_land = lands_N[[j]]
			this_species = species_N[[j]]
			this_gsad = gsad_N[[j]]
			
			# Distribute species across landscape
			this_metacomm = with(parms, {
				p = ifelse(exists('prop_full'), prop_full, NA)
				distribution = ifelse(exists('init_distribute'), init_distribute, NA)
				if(exists('cells_distribute')){
					which_cells = cells_distribute
				} else {
					which_cells = NULL
				}
				populate_landscape(this_land, this_species, this_gsad, K, distribution, p, which_cells)
			})
			
			# Run simulation
			with(parms, {
				if(!exists('d_kernel')) d_kernel = NULL
				if(!exists('v_kernel')) v_kernel = NULL
				imm_rate = ifelse(exists('imm_rate'), imm_rate, NA)
				if(!exists('save_steps')) save_steps = NULL
				run_sim(nsteps, this_metacomm, this_land, this_species, this_gsad, d_kernel, v_kernel, imm_rate, save_steps, report, ID=j)
			})
		})
	}

	# Save results
	if(!is.null(save_sim)) save(sim_results, lands_N, species_N, gsad_N, file=paste0(save_sim, simID, '.RData'))

	# Return results
	list(results = sim_results, species = species_N, lands = lands_N, gsads = gsad_N)
}




####################################################################
### Functions for collecting data from a simulation

# A function that groups cells together in a regular fashion and returns a list of cell locations.
#	X : number of rows or a vector of length 2 giving the dimensions of the grid
#	Y : number of columns
#	dX : number of rows to aggregate or a vector of length 2 giving the dimensions to aggregate
#	dY : number of columns to aggregate
#	form : string indicating how cell aggregations should be formed. Defaults to 'partition'.
#		'partition' = grid partitioned so that groups are non-overlapping
#		'window' = sliding window used. Results in overlapping groups. Default is 1 cell in each direction.
#		'origin' = cells aggregated around a focal cell or set of cells. Default is center cell.
#	slide : if form='window', the number of cells to slide the window in the x and y directions
#	locs : if form='origin', a matrix of cell locations that should be the center of aggregated groups.
aggregate_cells = function(X, Y=NA, dX, dY=NA, form=NA, slide=NULL, locs=NULL){
	
	# Catch errors where not enough parameter specified
	if(length(X)==1 & is.na(Y)) stop('Must specify dimensions of landscape.')
	if(length(dX)==1 & is.na(dY)) stop('Must specify aggregation dimensions dX x dY.')

	# Define missing parameters to defaults.
	if(length(X)==2){
		Y = X[2]
		X = X[1]
	}

	if(length(dX)==2){
		dY = dX[2]
		dX = dY[1]
	}

	if(is.na(form)) form = 'partition'
	
	if(form=='window'){
		if(is.null(slide)) slide = c(1,1)
		if(length(slide)==1) slide = rep(slide, 2)
	}

	if(form=='origin' & is.null(locs)){
		locs = matrix(c(ceiling(X/2), ceiling(Y/2)), nrow=1)		
	}
	
	# Aggregate cells
	if(form=='partition'){

		# Determine starting points
		Xs = seq(1, X, dX)
		Ys = seq(1, Y, dY)
		group_origins = as.matrix(expand.grid(Xs, Ys))

		# For each origin cell, find all cells within dX, dY
		groups = sapply(1:nrow(group_origins), function(i){
			o = group_origins[i,]
			Xcoords = o[1]:min(o[1]+dX-1, X)
			Ycoords = o[2]:min(o[2]+dY-1, Y)
			cells = as.matrix(expand.grid(Xcoords, Ycoords))
			colnames(cells) = c('x','y')
			list(cells)		
		})

	}
	if(form=='window'){
	
		# Determine starting points
		Xs = seq(1, X, slide[1])
		Ys = seq(1, Y, slide[2])
		group_origins = as.matrix(expand.grid(Xs, Ys))

		# For each origin cell, find all cells within dX, dY
		groups = sapply(1:nrow(group_origins), function(i){
			o = group_origins[i,]
			Xcoords = o[1]:min(o[1]+dX-1, X)
			Ycoords = o[2]:min(o[2]+dY-1, Y)
			cells = as.matrix(expand.grid(Xcoords, Ycoords))
			colnames(cells) = c('x','y')
			list(cells)		
		})
	}
	if(form=='origin'){

		# For each given origin cell, find all cells within dX, dY centered on that cell
		groups = sapply(1:nrow(locs), function(i){
			o = locs[i,]

			if(dX%%2==1){
				Xcoords =(o[1]-floor(dX/2)):(o[1] + floor(dX/2))
			} else {
				Xcoords = (o[1]-dX/2+1):(o[1]+dX/2)
			}

			if(dY%%2==1){
				Ycoords =(o[2]-floor(dY/2)):(o[2] + floor(dY/2))
			} else {
				Ycoords = (o[2]-dY/2+1):(o[2]+dY/2)
			}
			cells = as.matrix(expand.grid(Xcoords, Ycoords))
			colnames(cells) = c('x','y')
			list(cells)	
		})
	}


## DEBUG: CHECK GROUPS ##
#	image(0:X,0:Y,matrix(1, nrow=X, ncol=Y))
#	for(i in 1:length(groups)){
#		pts = groups[[i]]
#		text(pts[,1]-.5, pts[,2]-.5, labels=i)
#	}

	# Return list of locations
	groups
}


# A function that calculates species relative abundances in a landscape
#	metacomm : matrix of lists of species present in each cell, including 0 for empty spaces
#	N_S : number of species
#	only_species : whether empty spaces should be excluded from consideration. Defaults to FALSE.
calc_abun = function(metacomm, N_S, only_species=F){
	# Determine where to start species
	start = as.numeric(only_species)

	# Tally number of individuals of each species present, including empty spaces (0)
	abun = table(factor(unlist(metacomm), start:N_S))
	
	# Calculate relative abundance across all available spaces
	rel_abun = abun/sum(abun)

	# Return abundances
	rel_abun
}


# A function that calculates species abundance profiles during a given time window for a given set of locations.
# Returns an array of species abundances at each location through time [time, species, location]
#	locs : two-column matrix of cell locations where species occupancies should be calculated or a list of cell locations that should be aggregated.
#	t_window : either a list of start and stop times specifying all collected timepoints in a given interval or an explicit vector of timepoints to be considered
#	sim : an array of simulation results, as returned by run_sim() function
#	N_S : number of species
calc_abun_profile = function(locs, t_window, sim, N_S){
	
	# Determine which timepoints to evaluate
	timepoints = as.numeric(dimnames(sim)$time)
	if(is.list(t_window)){
		use_times = timepoints[timepoints >= t_window$start & timepoints <= t_window$stop]
	} else {
		# Catch error when a time is specified that was not recorded in simulation
		if(sum(t_window %in% timepoints) < length(t_window)){
			missing = t_window[!(t_window %in% timepoints)]
			stop(paste('Trying to measure occurance during timepoints not recorded in simulation. T =', paste(missing, collapse=' ')))
		}

		use_times = t_window
	}

	# Determine whether cells should be aggregated
	if(!is.list(locs)) locs = list(locs)
	
	# For each location:
	abun_profiles = sapply(locs, function(cell_block){
	
		# For each block of cells to be aggregated
		each_cell = sapply(1:nrow(cell_block), function(i){
			
			# Convert community profile to matrix
			x = as.numeric(cell_block[i,])
			comm_mat = simplify2array(sim[x[1],x[2],as.character(use_times)])

			# Calculate abundance of each species across timesteps
			abuns = sapply(0:N_S, function(sp) colSums(comm_mat==sp))		
			colnames(abuns) = 0:N_S

			# Return abundnaces
			abuns
		}, simplify='array')

		# Sum across cells and timepoints
		apply(each_cell, 1:2, sum)
		
	}, simplify='array')

	# Return abundance profiles
	abun_profiles
}


# A function that plots species abundance profiles
#	prof : a matrix of species abundance through time. Rows are timepoints and columns are species.
#	sp_type : a vector of chacters indicating whether species are specialists on 'A' or 'B' or generalists ('AB').
#	lcol : vector of length 2 or 3 giving line colors for specialists and generalists. Defaults to blue, red, purple
plot_abun = function(prof, lcol=c('royalblue','orchid','orangered')){
	
	# Total number of individuals in each time step
	N = rowSums(prof)

	# Number of timepoints and species
	times = 1:nrow(prof)
	sps = colnames(prof)
	
	# Convert profiles into cumulative profiles with empty spaces tallied at the end
	cum_prof = prof[,sps[sps!='0']]
	if('0' %in% sps) cum_prof = cbind(cum_prof, prof[,'0'])
	cum_prof = apply(cum_prof, 1, function(x){
		sapply(1:length(x), function(i) sum(x[1:i]) )
	})

	if('0' %in% sps){
		rownames(cum_prof) =  c(sps[sps!='0'], '0')
	} else {
		rownames(cum_prof) = sps
	}

	# Define colors
	names(lcol) = levels(factor(sp_type))

	# Set up plot and axes
	plot.new()
	plot.window(xlim=c(.5, nrow(prof)+.5), ylim=c(0,max(N))) 

	axis(1, at=times, labels=rownames(prof))
	mtext('Time', 1, 2.5)
	abline(h=par('usr')[3], lwd=3)
	axis(2, las=1)
	abline(v=par('usr')[1], lwd=3)
	mtext('Num. Individuals', 2, 3)

	# Add lines for each species colored by type
	for(i in rownames(cum_prof)){
		if(i=='0'){
			lines(times, cum_prof[i,], col='black')
		} else {
			lines(times, cum_prof[i,], col=lcol[sp_type[i]], lend=1)
		}	
	}
}


# A function that calculates species occupancies during a given time window for a given set of locations
#	locs : two-column matrix of cell locations where species occupancies should be calculated or a list of cell locations that should be aggregated.
#	t_window : either a list of start and stop times specifying all collected timepoints in a given interval or an explicit vector of timepoints to be considered
#	sim : an array of simulation results, as returned by run_sim() function
#	N_S : number of species
#	abuns : array of species abundance profiles, as generated by calc_abun_profiles() function.
# 	agg_times : either a single number specifing the number of timepoints that should be aggregated before calculating occupancy or a list defining exactly which timepoints should be aggregated. Timepoints are relative t_window. Defaults to no aggregation.
#	which_species : a vector indicating which species should be examined. Defaults to all species.
#	do_freq : logical indicating whether to return actual frequencies rather than occupancy (scaled by number of timepoints). Defaults to FALSE.
calc_occupancy = function(locs=NULL, t_window=NULL, sim=NULL, N_S=NULL, abuns=NULL, agg_times=NULL, which_species=NULL, do_freq=F){
	
	# Catch error where not enough information is specified
	if(is.null(abuns)&(is.null(locs)|is.null(t_window)|is.null(sim)|is.null(N_S))) stop('Must specify abuns or locs/t_window/sim/N_S.')	

	# Calculate species abundance profiles if not specified
	if(is.null(abuns)){
		abun_profiles = calc_abun_profile(locs, t_window, sim, N_S)
	} else {
		abun_profiles = abuns
	}

	# Number of timepoints measured
	N_t = dim(abun_profiles)[1]	

	# Determine which species to examine
	if(is.null(which_species)) which_species = colnames(abun_profiles)[2:ncol(abun_profiles)]

	# Aggregate timepoints, if specified
	if(!is.null(agg_times)){
		
		# Determine which times to aggregate
		if(is.numeric(agg_times)){
			if(length(agg_times)>1) stop('To specify non-uniform aggregation use list format.')
			int = agg_times
			agg_times = lapply(seq(1,N_t,int), function(x) x:(x+int-1))
		}

		# Calculate summed abundances across aggregation times
		# Note that new array will have dimensions  [species, locations, times]
		agg_abuns = sapply(agg_times, function(rows){
			rows = rows[rows %in% 1:N_t]
			if(length(rows)>1) new_abun = apply(abun_profiles[rows,,], 2:3, sum)
			if(length(rows)==1) new_abun = abun_profiles[rows,,]
			new_abun
		}, simplify='array')

		# Calculate number of timepoints that each specified species is present
		if(length(agg_times)==1){
			freqs = t(ifelse(agg_abuns[as.character(which_species),,]>0, 1, 0))
		} else {
			if(length(which_species)>1){
				freqs = apply(agg_abuns[as.character(which_species),,] > 0, 1:2, sum)
			} else {
				freqs = t(apply(agg_abuns[as.character(which_species),,] > 0, 1, sum))
			}
		}

		# Calculate occupancy
		if(!do_freq) freqs = freqs/length(agg_times)

	} else {

		# Calculate number of timepoints that each specified species is present
		if(length(which_species)>1){
			freqs = apply(abun_profiles[,as.character(which_species),] > 0, 2:3, sum)
		} else {
			freqs = t(apply(abun_profiles[,as.character(which_species),] > 0, 2, sum))
		}
	
		# Calculate occupancy
		if(!do_freq) freqs = freqs/N_t
	}

	# Return results as species x site matrix
	t(freqs)
}

# A function that cross-classifies species by their true vs. observed core/transient status.
# Returns a contingency table whose rows are the true classification based on birth rates and columns are the observed classes based on occupancy.
# Can specify a classification matrix or birth rates and habitat types
#	occupancy : matrix of species occupancies across spatial units (sites X species).
#	breaks : either a numeric vector of lengtth 1 or 2 giving occupancy breakpoints for transient vs. core status or a named list wuth specific intervals: list(trans=c(), core=c())
#	b_rates : matrix of birth rates for each species in each habitat types.
#	habitats : vector of habitat types ('A' or 'B') for each spatial unit
#	classification : matrix classifying each species on each site as 'core' or 'trans'. (sites X species)
# 	do_each : logical indicating whether one table should be returned across all sites (F) or an array of tables should be returned, one for each spatial unit (T).
#	return : character indicating whether a table of 'counts' or lists of species IDs ('ids') should be returned. Default is 'counts'.
cross_classify = function(occupancy, breaks, b_rates=NULL, habitats=NULL, classification=NULL, do_each=F, return='counts'){

	# Convert breaks to list of intervals defining core and transient occupancy levels if just specified as numeric
	if(!is.list(breaks)){
		if(length(breaks)==1) breaks = rep(breaks,2)
		breaks = list(trans=c(0,breaks[1]), core=c(breaks[2], 1))
	}

	# Classify species based on occupancy matrix
	
	if(breaks$trans[2]==breaks$core[1]){
		
		# Case when no occupancy levels excluded
		classes = matrix(cut(occupancy, breaks=c(breaks$trans, breaks$core[2]), include.lowest=T, labels=c('trans','core')), nrow=nrow(occupancy))
	
	} else {
		
		# Case when intermediate occupancy levels excluded
		classes = matrix(cut(occupancy, breaks=c(breaks$trans, breaks$core), include.lowest=T, labels=c('trans',NA,'core')), nrow=nrow(occupancy))
	}
	colnames(classes) = colnames(occupancy)

	# Get predefined classification or calculate from birth rates
	if(is.null(classification)){
		
		# Catch error where not enough information specified
		if(is.null(habitats)&is.null(b_rates)) stop('Must specify either a classification of each species in each spatial unit or habitat types and birth rates.')

		cores = t(sapply(habitats, function(h) b_rates[,h]>0))
		classification = apply(cores, 1:2, function(x) ifelse(x, 'core', 'trans'))		
	}

	# Check whether classification matrix matches occupancy matrix
	if(sum(dim(classification)==dim(classes))<2) stop('Incorrect dimensions for classification matrix. Should be sites x species.')

	# Calculate contingency tables
	if(do_each){
		# Calculate for each spatial unit
		if(return=='counts') cross_tab = sapply(1:nrow(classes), function(i) table(classification=classification[i,], occupancy=classes[i,]), simplify='array')
		if(return=='ids'){
			cross_tab = sapply(1:nrow(classes), function(i){
				this_tab = array(list(), dim=c(2,2), dimnames=list(classification=c('core','trans'), occupancy=c('core','trans')))
				this_tab['core','core'] = list(as.character(which(classes[i,]=='core' & classification[i,]=='core')))
				this_tab['core','trans'] = list(as.character(which(classes[i,]=='trans' & classification[i,]=='core')))
				this_tab['trans','core'] = list(as.character(which(classes[i,]=='core' & classification[i,]=='trans')))
				this_tab['trans','trans'] = list(as.character(which(classes[i,]=='trans' & classification[i,]=='trans')))
				this_tab
			}, simplify='array')
		}
	} else {
		# Sumarize across all spatial units
		if(return=='count') cross_tab = table(classification=classification, occupancy=classes)
		if(return=='ids'){
			cross_tab = array(list(), dim=c(2,2), dimnames=list(classification=c('core','trans'), occupancy=c('core','trans')))
			cross_tab['core','core'] = list(as.character(which(classification=='core'&classes=='core')))
			cross_tab['core','trans'] = list(as.character(which(classification=='core'&classes=='trans')))
			cross_tab['trans','core'] = list(as.character(which(classification=='trans'&classes=='core')))
			cross_tab['trans','trans'] = list(as.character(which(classification=='trans'&classes=='trans')))
		}
	}

	# Return tables
	cross_tab
}



# A function that calculates species richness profiles through time
#	locs : two-column matrix of cell locations where species occupancies should be calculated or a list of cell locations that should be aggregated.
#	t_window : either a list of start and stop times specifying all collected timepoints in a given interval or an explicit vector of timepoints to be considered
#	sim : an array of simulation results, as returned by run_sim() function
#	N_S : number of species
#	abuns : array of species abundance profiles, as generated by calc_abun_profiles() function.
# 	agg_times : either a single number specifying the number of timepoints that should be aggregated before calculating occupancy or a list defining exactly which timepoints should be aggregated. Timepoints are relative t_window. Defaults to no aggregation.
#	which_species : a vector indicating which species should be examined. Defaults to all species.
calc_rich = function(locs=NULL, t_window=NULL, sim=NULL, N_S=NULL, abuns=NULL, agg_times=NULL, which_species=NULL){
	
	# Catch error where not enough information is specified
	if(is.null(abuns)&(is.null(locs)|is.null(t_window)|is.null(sim)|is.null(N_S))) stop('Must specify abuns or locs/t_window/sim/N_S.')	

	# Calculate species abundance profiles if not specified
	if(is.null(abuns)){
		abun_profiles = calc_abun_profile(locs, t_window, sim, N_S)
	} else {
		abun_profiles = abuns
	}

	# Convert to array if necessary
	if(length(dim(abun_profiles))==2) abun_profiles = add_dim(abun_profiles)

	# Determine which species to examine
	if(is.null(which_species)) which_species = colnames(abun_profiles)[2:ncol(abun_profiles)]

	# Number of timepoints, species, and sites measured
	N_t = dim(abun_profiles)[1]	
	N_s = length(which_species)
	N_c = dim(abun_profiles)[3]

	# Aggregate timepoints, if specified
	if(!is.null(agg_times)){
		
		# Determine which times to aggregate
		if(is.numeric(agg_times)){
			if(length(agg_times)>1) stop('To specify non-uniform aggregation use list format.')
			int = agg_times
			agg_times = lapply(seq(1,N_t,int), function(x) x:(x+int-1))
		}

		# Calculate summed abundances across aggregation times
		# Note that new array will have dimensions  [species, locations, times]
		agg_abuns = sapply(agg_times, function(rows){
			rows = rows[rows %in% 1:N_t] 
			use_abun = abun_profiles[rows,,,drop=F]
			if(length(rows)>1){
				new_abun = apply(use_abun, 2:3, sum)
			}
			if(length(rows)==1) new_abun = use_abun
			new_abun
		}, simplify='array')
		if(N_c==1 & length(agg_times)>1 & length(dim(agg_abuns))==2) agg_abuns = add_dim(agg_abuns, 2)
		dimnames(agg_abuns) = list(dimnames(abun_profiles)[[2]], dimnames(abun_profiles)[[3]], 1:length(agg_times))
		

		# Calculate richness in each aggregation interval
		rich = apply(agg_abuns[as.character(which_species),,,drop=F] > 0, 2, colSums)

	} else {

		# Calculate richness at each timepoint
		rich = apply(abun_profiles[,as.character(which_species),,drop=F] > 0, 3, rowSums)
	}

	# Return results as species x site matrix
	rich
}

# A function that calculates richness of core and transient species
#	abuns : array of abundance profiles returned by calc_abun_profile() function
#	occupancy : matrix of species occupancies across the same spatial units represented in abuns, as returned by calc_occupancy(abuns, do_freq=F)
#	breaks : vector of ordered numbers between 0 and 1 denoting breakpoints for occupancy categories
# 	agg_times : either a single number specifying the number of timepoints that should be aggregated before calculating richness or a list defining exactly which timepoints should be aggregated. Timepoints are relative to times in abuns. Defaults to no aggregation.
calc_rich_CT = function(abuns, occupancy, breaks, agg_times=NULL){

	# Catch error when abuns and occupancy do not match
	if(nrow(occupancy)!=dim(abuns)[length(dim(abuns))]) stop('Abundance matrix and occupancy matrix must be supplied for same set of sites.')
	
	# Convert occupancy to factor based on breakpoints
	if(breaks[1]!=0) breaks = c(0,breaks)
	if(breaks[length(breaks)]!=1) breaks = c(breaks, 1)
	cats = matrix(cut(occupancy, breaks, include.lowest=F, labels=F), nrow=nrow(occupancy), ncol=ncol(occupancy))
	colnames(cats) = colnames(occupancy)
		
	# Number of sites and time windows
	N_c = dim(abuns)[3]
	N_t = ifelse(is.null(agg_times), dim(abuns)[1], length(agg_times))

	# Determine time aggregation windows
	if(is.numeric(agg_times)){
		if(length(agg_times)>1) stop('To specify non-uniform aggregation use list format.')
		int = agg_times
		agg_times = lapply(seq(1,N_t,int), function(x) x:(x+int-1))
	}

	# Richness of each category for each site during given aggregated time periods. Returns [times, categories, sites]
	rich_cat = sapply(1:N_c, function(i){
		sapply(1:(length(breaks)-1), function(cat){
			these_sp = which(cats[i,]==cat)
			if(length(these_sp)>0){
				calc_rich(abuns = abuns[,,i], agg_times = agg_times, which_species = these_sp)
			} else {
				rep(0, N_t)
			}
		})
	}, simplify='array')
	dimnames(rich_cat) = list(1:N_t, levels(cut(0, breaks)), 1:N_c)

	# Return richnes
	rich_cat
}

# A function that calculates the average habitat value for a set of cells
# 	locs : a matrix of cell locations
#	land : landscape matrix
average_habitat = function(locs, land){

	# Get mean value across all cells
	vals  = apply(locs, 1, function(i) land[i[1], i[2]])

	# Determine habitat type
	# Note: equal parts 'A' and 'B' gets classified as 'A'
	get_habitat(mean(vals))
}

# A function that returns an incomplete sample of a simulation based on species detectibility.
# 	abuns : species abundance profiles, as returned by calc_abun_profile(): [timepoints, species, sites]
# 	probs : probability that an individual is detected. Can be either a single detectability for all species, or a vector with different probabilies for each species. Defaults to 1.
#	return : a string indicating whether species abundances 'abundance' or presence 'presence' should be returned. Default is abundance.
sample_sim = function(abuns, probs = NULL, return='abundance'){
	
	# Drop empty spaces from abuns, if present.
	abuns = abuns[,dimnames(abuns)[[2]]!='0',]

	# Determine number of species
	N_S = dim(abuns)[2]

	# If not specified, detection is 1 and the original abundance profiles are returned
	if(is.null(probs)){
		obs = abuns
		if(return=='presence') obs = abuns>0
	} else {
		
		# Make vector of detectabilities if only one specified.
		if(length(probs)==1) probs = rep(probs, N_S)
		
		# If species presence to be returned		
		if(return=='presence'){

			# Calculate probability of observing species (1 - P(not observing))
			# Returns an array: [timepoints, sites, species]
			P = sapply(1:N_S, function(sp){
				apply(abuns[,sp,], 1:2, function(x) 1 - (1-probs[sp])^x )
			}, simplify='array')

			# Stochastically determine which species are observed
			rands = array(runif(length(P)), dim=dim(P))
			obs = rands <= P

			# Rearrange dimensions to match abuns
			obs = aperm(obs, c(1,3,2))	
		}
		
		if(return=='abundance'){
			
			# Stochasitically determine which individuals are observed
			obs = sapply(1:N_S, function(sp){
				apply(abuns[,sp,], 1:2, function(x){
					if(x>0){
						sum(runif(x) <= probs[sp])
					} else {
						0
					}

				})
			}, simplify='array')
			
			# Rearrange dimensions to match abuns
			obs = aperm(obs, c(1,3,2))
		}
	}

	# Return observed sample
	obs
} 


# DEBUGGING
#sim = sim_results$results[[1]]
#locs = scale_locs[[3]]
#t_window=list(start=76, stop=100)
#species = sim_results$species[[1]]
# habitats = sapply(locs, function(x) average_habitat(x, sim_results$lands[[1]]))

#t_window : timepoints to consider
#agg_times_occ : timepoints to aggregate before calculating occupancy
#agg_times_sum : timepoints to aggregate before calculating richness or abundance

summarize_sim = function(what, sim, species, land, gsad, locs, t_window, agg_times_occ, agg_times_sum, breaks, P_obs=NULL){
	
	# Get species birth rates
	b_rates = species[,,'b']
	
	# Get habitat types for spatial units
	habitats = sapply(locs, function(x) average_habitat(x, land))

	# Calculate classification of species based on birth rates
	cores = t(sapply(habitats, function(h) b_rates[,h]>0))
	classification = apply(cores, 1:2, function(x) ifelse(x, 'core', 'trans'))	

	# Calculate abundance profiles
	abuns_act = calc_abun_profile(locs, t_window, sim, dim(species)[1])
	
	# Calculate occupancy
	occ_act = calc_occupancy(abuns=abuns_act, agg_times=agg_times_occ, do_freq=F)

	# Cross-classify species in spatial units
	# May not need this
	cross_class_ids = cross_classify(occ_act, breaks, classification=classification, do_each=T, return='ids')

	# Calculate abundance of each class of species across spatial units
	

	# Calculate richness of each occupancy class across spatial units
	rich_occ_act = calc_rich_CT(abuns_act, occ_act, breaks, agg_times_sum)

	# Calculate abundance of each occupancy class across spatial units

	# Calculate richness of each biological classification across spatial units
	rich_bio_act = calc_rich_CT(abuns_act, apply(classification, 1:2, function(x) ifelse(x=='core', 1, 0)), breaks, agg_times_sum)

	# Calculate means across temporal units
	# THESE ARE NOT MATCHING UP_ DOES THIS MAKE SENSE??
	#apply(rich_occ_act, 2:3, mean)
	#apply(rich_bio_act, 2:3, mean)


	# If detection probability is given, sample the abundance profiles and return comparison between actual and observed
	if(!is.null(P_obs)){



	 abuns_obs = sample_sim(abuns_act, probs=P_obs, return='abundance')
	}


	


}



# Basic statistics:
# Abundance and richness of core and transient species, based on occupancy and based on birth rates
# Relationship between abundance and occupancy across species
# 

# Each statistic calculated at several spatial scales, including the entire grid.



# A function that summarizes results from multiple replicates of a simulation.
# Also works for a simgle simulation run
#	results : either a list of simulation results return by run_sim_N, a list of simulation results metacommunities from different runs, or a single simulation results metacommunity
summarize_sim_N = function(results, speciesN=NULL, landN=NULL, t_window=NULL, locs=NULL, agg_times=NULL, which_species=NULL, do_freq=F){

	# Check whether species are given
	if(!with(results, exists('species'))&is.null(speciesN)) stop('Must specify species vital rates lists if not included in results.')	
	if(!with(results, exists('lands'))&is.null(landN)) stop('Must specify landscapes lists if not included in results.')	

	# Assign species to internal objects
	if(with(results, exists('results'))){
		simN = results$results
		speciesN = results$species
		landN = results$lands
	} else {
		
		if(is.null(speciesN)) stop('Must specify species vital rates.')
		if(is.null(landN)) stop('Must specify landscapes.')

		# If this is a single simulation convert to a list
		if(!is.null(dim(results))){
			simN = list(results)
			speciesN = list(speciesN)
			landN = list(landN)
		} else {
			simN = results
		}
	}

	# For unspecified parameters define to default values
	if(is.null(t_window)){
		times = dimnames(simN[[1]])$time
		t_window = list(start=as.numeric(times[1]), stop=as.numeric(times[length(times)]))
	} 
	if(is.null(locs)){
		X = dim(results[[1]])[1]
		Y = dim(results[[1]])[2]
		locs = as.matrix(expand.grid(x=1:X, y=1:Y))
	}

	# Determine number of simulation replicates
	nsim = length(simN)

	# For each simulation replicate
	sapply(1:nsim, function(k){

		# Define this simulation
		this_sim = simN[[k]]
		this_species = speciesN	[[k]]
		this_land = landN[[k]]	
		N_S = dim(this_species)[1]
		
############ WORKING HERE #####################
		# Calculate species abundance profiles
		abun_profs = calc_abun_profile(locs, t_window, this_sim, N_S)

		# Calculate species occupancy
		occupancy = calc_occupancy(abuns=abun_profs, agg_times, which_species, do_freq=F)

		# Calculate species total and abundances across landscape
#		rel_abun = sapply(calc_abun(this_sim[,,step], N_S, only_species=T)
		
		# Calculate proportion of metacommunity that is full

	})



}


#######################################################
### Input/Output Functions

# A function that makes a list of parameters for writing to a file or passing to run_sim_N() from a given environment
# This function needs to be updated whenever new parameters are coded into the simulation
#	e : environment where parameter values can be found. Defaults to parent environment.
make_parmlist = function(e=parent.frame()){

	# Define list of required parameters
	parms = list(
		dimX = e$dimX,
		dimY = e$dimY,
		S_A = e$S_A, 
		S_B = e$S_B,
		m_rates = e$m_rates,
		r_rates = e$r_rates,
		K = e$K,								
		nsteps = e$nsteps,
		nruns = e$nruns
	)

	# Add on optional parameters
	if(exists('vgm_dcorr', e)) parms = c(parms, vgm_dcorr = e$vgm_dcorr)			
	if(exists('vgm_mod', e)) parms = c(parms, vgm_mod = list(e$vgm_mod)) 
	if(exists('habA_prop', e)) parms = c(parms, habA_prop = e$habA_prop)
	if(exists('S_AB', e)) parms = c(parms, S_AB = e$S_AB)
	if(exists('dist_b', e)) parms = c(parms, dist_b = list(e$dist_b))
	if(exists('dist_d', e)) parms = c(parms, dist_d = list(e$dist_d))
	if(exists('dist_v', e)) parms = c(parms, dist_v=list(e$dist_v))
	if(exists('dist_gsad',e)){
		if(is.list(e$dist_gsad)){
			parms = c(parms, dist_gsad = list(e$dist_gsad))
		} else {
			parms = c(parms, dist_gsad = e$dist_gsad)
		}
	}
	if(exists('prop_full', e)) parms = c(parms, prop_full = e$prop_full)
	if(exists('init_distribute', e)) parms = c(parms, init_distribute = e$init_distribute)
	if(exists('cells_distribute', e)) parms = c(parms, cells_distribute = list(e$cells_distribute))
	if(exists('d_kernel', e)) parms = c(parms, d_kernel = list(e$d_kernel))
	if(exists('v_kernel', e)) parms = c(parms, v_kernel = list(e$v_kernel))
	if(exists('imm_rate', e)) parms = c(parms, imm_rate = e$imm_rate)
	if(exists('save_steps', e)) parms = c(parms, save_steps = list(e$save_steps))
	if(exists('simID', e)) parms = c(parms, simID = e$simID)
}






#######################################################
#### Miscellaneous Functions


# Function that adds a dimension on to an array
add_dim = function(x, where=NULL){

	if(is.null(where)) where = length(dim(x)) + 1
	
	# Get original max dimension
	n = ifelse(is.null(dim(x)), 1, length(dim(x))) 
	
	# Determine whether original has dimnames
	if(n>1) has_names = !is.null(dimnames(x))
	if(n==1) has_names = !is.null(names(x))
	
	# Make new dimensions and names
	new_dim = as.numeric(1:(length(where)+n) %in% where)
	new_names = lapply(1:length(new_dim), function(i) ifelse(i %in% where, 1, NA))
	old_names = which(new_dim==0)
	
	if(n==1){
		if(has_names) new_names[[old_names]] = names(x) 
		new_dim[new_dim==0] = length(x)
	} else {
		if(has_names) for(i in 1:length(old_names)) new_names[[old_names[i]]] = dimnames(x)[[i]] 
		new_dim[new_dim==0] = dim(x)
	}

	if(has_names){
		new_array = array(x, dim=new_dim, dimnames=new_names)
	} else {
		new_array = array(x, dim=new_dim)
	}

	# Return new array
	new_array	
}












