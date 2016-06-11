## This script contains functions for running the Core-Transient Simulation

## Required packages
library(gstat)
library(sp)
library(raster)
library(poweRlaw)
library(fdrtool) # half-normal distribution for gaussian dispersal kernal



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
# 	r : vector of length 2 or 3 specifying recruitment rates [preferred habitat, not preferred habitat, generalist rate]
# 	dist_d : named list of parameters describing distribution from which dispersal kernals should be drawn
#		mu = mean dispersal distance
#		var = variance of distribution from which kernals drawn. Defaults to 0 for same dispersal kernal for all species.
make_species = function(S_A=NULL, S_B=NULL, S_AB=0, dist_b = list(maxN=10, type='uniform'), m, r, dist_d=list(mu=1, var=0)){

	# Catch error if no species specified
	if(is.null(S_A)|is.null(S_B)) stop('Must specify number of species.')

	# Create array to hold species vital rates (b = birth, m = death, r = recruitment, d = dispersal)
	# Rates are per timestep
	N_S = S_A + S_B + S_AB
	species_rates = array(NA, dim=c(N_S, 2, 4), dimnames=list(species=1:N_S, habitat=c('A','B'), rate=c('b','m','r','d')))

	# Specialist species birth rates in their preferred habitat ranged from 1 to maxN and are 0 in the unpreferred habitat
	species_rates[1:S_A, 'A', 'b'] = make_sad(S_A, dist_b)
	species_rates[1:S_A, 'B', 'b'] = 0
	species_rates[(S_A+1):(S_A+S_B), 'B', 'b'] = make_sad(S_B, dist_b)
	species_rates[(S_A+1):(S_A+S_B), 'A', 'b'] = 0

	# Generalist birth rates range from 1 to maxN in both habitats, but is not necessarily the same in both habitats
	if(S_AB > 0) species_rates[(N_S-S_AB+1):N_S, ,'b'] = make_sad(2*S_AB, dist_b)

	# Death rates are prespecified constants
	# This could be modified to be drawn from a distribution
	species_rates[1:S_A, 'A', 'm'] = m[1]
	species_rates[1:S_A, 'B', 'm'] = m[2]
	species_rates[(S_A+1):(S_A+S_B), 'B', 'm'] = m[1]
	species_rates[(S_A+1):(S_A+S_B), 'A', 'm'] = m[2]

	# Generalist death rates default to the death rate for the preferred habitat, but can be specified separately
	if(S_AB > 0) species_rates[(N_S-S_AB+1):N_S, ,'m'] = ifelse(length(m)==3, m[3], m[1])

	# Specialist species recruitment rates are prespecified constants
	# This could be modified to be drawn from a distribution
	species_rates[1:S_A, 'A', 'r'] = r[1]
	species_rates[1:S_A, 'B', 'r'] = r[2]
	species_rates[(S_A+1):(S_A+S_B), 'B', 'r'] = r[1]
	species_rates[(S_A+1):(S_A+S_B), 'A', 'r'] = r[2]

	# Generalist recruitment rates default to the recruitment rate for the preferred habitat, but can be specified separately
	if(S_AB > 0) species_rates[(N_S-S_AB+1):N_S, ,'r'] = ifelse(length(r)==3, r[3], r[1])

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
myland = make_landscape(30,30)
mysp = make_species(S_A=10, S_B=15, S_AB=3, dist_b = list(maxN=5, P_maxN=0.01, type='poisson'), m=c(0.1, 0.5), r =c(.9,.7), dist_d=list(mu=3, var=0.5))
mycomm = populate_landscape(myland, mysp, K=100, p=.5, distribution='uniform')
mygsad = make_sad(dim(mysp)[1], distribution=list(type='same'))

####################################################################
### Functions for running a simulation

# Function that determines habitat type based on habitat value stored in landscape
# 	x : numeric habitat value stored in landscape
get_habitat = function(x){
	h = NA
	
	if(x==-1) h = 'A'
	if(x==1) h = 'B'

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
disperse = function(x, y=NULL, propagules, dim_land, d_rates, form=list(type='gaussian')){

	# If only x specified, define y
	if(length(x)==2){	
		x = x[1]
		y = x[2]
	}

	# Define limits of origin cell
	origin_cell = rbind(c(x-1, x), c(y-1, y))	

	# For each propagule
	new_locs = sapply(propagules, function(i){
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

		destination_cell
	})

	# Return table of new locations
	t(new_locs)

}


# Function that stochastically generates a dispersal vector (distance, direction) for given a dispersal kernal
#	d : expected distance (mean)
#	form : a list defininting the type of dispersal kernel. Currently only implements 'gaussian'.
#		type = character indicating the type
#	N : number of vectors to generate. Defaults to 1
get_dispersal_vec = function(d, form=list(type='gaussian'), N=1){
	
	# Gaussian kernel
	if(form$type=='gaussian'){
	
		# Calculate theta for mean = d
		theta = 1/d

		# Generate distance
		r = rhalfnorm(N, theta)		
	}

	# Generate random angle (in radians)
	phi = runif(N, 0, 2*pi)

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
	}
	
	# Select a potential recruit from the pool of propagules
	recruits[!migrants] = sample(propagules)[1:sum(!migrants)]

	# Probabilistically determine whether propagules establish based on species recruitment rates for this habitat
	established = sapply(recruits, function(sp){
		ifelse(runif(1) > r_rates[sp] | is.na(sp), 0, sp)
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

run_timestep = function(metacomm, land, species, gsad, d_kernel=list(type='gaussian'), imm_rate=0){

	# Define dimensions of landscape
	X = nrow(land)
	Y = ncol(land)

	# Define array to hold new propagule pools
	propagule_pools = matrix(list(), nrow=X, ncol=Y)

	# For each cell, new individuals are born and disperse
	for(i in 1:X){
	for(j in 1:Y){

		# Define this community and habitat
		this_comm = metacomm[i,j][[1]]
		this_habitat = get_habitat(land[i,j])
		
		# New individuals born
		this_pool = reproduce(this_comm, species[,this_habitat,'b'])
	
		# Propagules disperse to new pools
		new_locs = disperse(i, j, this_pool, c(X,Y), species[,this_habitat,'d'], form=d_kernel)
		for(p in 1:nrow(new_locs)){
			# Propagules not that disperse off of the landscape are removed
			if(new_locs[p,1]>0 & new_locs[p,2]>0 & new_locs[p,1]<=X & new_locs[p,2]<=Y ) propagule_pools[new_locs[p,1], new_locs[p,2]] = list(unlist(c(propagule_pools[new_locs[p,1], new_locs[p,2]], this_pool[p])))
		}
	}}
	
	# Recruits establish from new propagule pools and then all community members experience stochastic mortality
	for(i in 1:X){
	for(j in 1:Y){

		# Define this community, habitat and pool of propagules
		this_comm = metacomm[i,j][[1]]
		this_habitat = get_habitat(land[i,j])
		this_pool = propagule_pools[i,j][[1]]
	
		# Determine which propagules estabish
		new_comm = list(establish(this_comm, this_pool, species[, this_habitat, 'r'], imm_rate, gsad))

		# Stochasitic mortality
		new_comm = die(new_comm, species[, this_habitat, 'm'])

		# Save results
		metacomm[i,j] = list(as.numeric(new_comm))
	}}

	# Return new metacommunity
	metacomm
}


# A function that runs a simulation for a given number of time steps
run_sim = function(steps, metacomm, land, species, gsad, parm_list, save_steps = NULL){
	
	# Define simulation
	X = nrow(land)
	Y = ncol(land)

	# Create array to save simulations results
	if(is.null(save_steps)) save_steps = 0:steps
	sim_results = array(list(), dim=c(X, Y, length(save_steps)),
		dimnames=list(row=1:X, col=1:X, time=save_steps))

	# Run simulation
	new_metacomm = metacomm
	for(step in 1:steps){
		new_metacomm = run_timestep(new_metacomm, land, species, gsad, parm_list$d_kernel, parm_list$imm_rate)
		if(step %in% save_steps) sim_results[,,as.character(step)] = new_metacomm
	}

	# Return results
	sim_results
}


parm_list = list(
	d_kernel=list(type='gaussian'),
	imm_rate = 0.2
)


# TESTING

test_run = run_sim(10, mycomm, myland, mysp, mygsad, parm_list, save_steps=seq(2,10,2))

####################################################################
### Functions for collecting data from a simulation


# A function that calculates species abundance profiles during a given time window for a given set of locations
#	locs : two-column matrix of cell locations where species occupancies should be calculated
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
	
	# For each location:
	abun_profiles = sapply(1:nrow(locs), function(i){
		# Convert community profile to matrix
		x = as.numeric(locs[i,])
		comm_mat = simplify2array(sim[x[1],x[2],])

		# Calculate abundance of each species across timesteps
		abuns = sapply(0:N_S, function(sp) colSums(comm_mat==sp))		
		colnames(abuns) = 0:N_S

		# Return abundnaces
		abuns
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
#	locs : two-column matrix of cell locations where species occupancies should be calculated
#	t_window : either a list of start and stop times specifying all collected timepoints in a given interval or an explicit vector of timepoints to be considered
#	sim : an array of simulation results, as returned by run_sim() function
#	N_S : number of species
calc_occupancy = function(locs, t_window, sim, N_S){
	
	# Calculate species abundance profiles
	abun_profiles = calc_abun_profile(locs, t_window, sim, N_S)

	### WORKING HERE


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

	apply(locs, 1, function(x){
		x = as.numeric(x)
		sim[x[1],x[2],]		

	})





}

# A function that calculates species relative abundances in a landscape
#	metacomm : matrix of lists of species present in each cell, including 0 for empty spaces
#	N_S : number of species
calc_abun = function(metacomm, N_S){
	# Tally number of individuals of each species present, including empty spaces (0)
	abun = table(factor(unlist(metacomm), 0:N_S))

	# Calculate relative abundance across all available spaces
	rel_abun = abun/sum(abun)

	# Return abundances
	rel_abun[as.character(1:N_S)]
}


# TESTING
abuns = apply(test_run, 3, function(x) calc_abun(x, 28))
b_rates = rowSums(species[,,'b'])
plot(b_rates, abuns[,'10'])
hist(abuns[,'10'], breaks=10)
 





for(i in 1:nrow(locs)){
	plot_abun(abun_profiles[,,i])
	Sys.sleep(2)
}



abun_profiles









