#' Summarize landscape habitats
#'
#' Summarizes habitat characteristics at a set of locations across a landscape.
#'
#' The function tabulates habitat types found in each spatial unit specified in
#' \code{locs} and each spatial unit is assigned to its dominant habitat type.
#' Then, the function calculates the proportion of spatial units in the 
#' dominant habitat type and the variance in habitat across spatial units
#' (assuming a Bernoulli distribution: \code{var = p(1-p)}). Finally, the function
#' calculates a measure of habitat heterogeneity from 
#' \href{http://dx.doi.org/10.1046/j.1440-1703.2000.00326.x}{Shiyomi 
#' and Yoshimura (2000), Ecological Research 15: 13-20},
#' which is a variance:mean ratio corrected for finite count data. 
#' For each spatial unit it counts the number of cells that are the dominant 
#' habitat typeand then compares the variance and mean of these frequencies. 
#' A value of zero indicates random distribution at that scale.
#'
#' @param sim (required) either a filename for a simulation run (as saved by
#' 	\code{\link{run_sim_N}} or a raster landscape (as returned by 
#' 	\code{\link{make_landscape}}. 
#' @param locs (required) either a two-column matrix of cell locations defining 
#' 	a single spatial unit or a list of matrices designating
#' 	multiple spatial units
#' @return a named vector with the dominant habitat proportion, habitat 
#' 	variance, and habitat heterogeneity (see details)
#' 
#' @export

summarize_land = function(sim, locs){

	# If sim is a file, then read in simulation run.
	if(is.character(sim)){
		load(sim)
		land = this_land
	}
	
	# If sim is a landscape
	if(class(sim) %in% c('RasterLayer', 'matrix')){
		land = sim
	}

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

	# Calculate proportion of dominant habitat cover
	hab_count = table(habitats)	
	hab_prop = hab_count[1]/sum(hab_count)
	if(hab_prop < 0.5) hab_prop = hab_count[2]/sum(hab_count)

	# Calculate variance- assumes Bernoulli distribution
	hab_var = hab_prop*(1-hab_prop)

	# DECIDED NOT TO USE SPATIALLY EXPLICIT MEASURES (E.G. REQUIRING DISTANCES) B/C LOCS ALLOWS FLEXIBLE SPECIFICATION OF SPATIAL UNITS ARE NOT NECESSARITY ON REGULAR GRID
	# INSTEAD USING MEASURE OF HETEROGENEITY FROM Shiyomi and Yoshimura 2000 Ecological Research 15: 13-20
	# Calculate heterogeneity if all spatial units sample the same number of cells
	if(length(unique(hab_tab['n',]))==1){
		# Dominant habitat type
		use_hab = names(hab_prop)
	
		# Number of cells in each spatial unit
		n = hab_tab['n',1]

		# Calculate mean and variance of frequency of dominant habitat type
		m = mean(hab_tab[use_hab,])
		v = var(hab_tab[use_hab,])

		# Calculate index based on variance to mean ratio
		hab_het = (v/(m*(1-(m/n)))) - 1
	} else {
		hab_het=NA	
	}

	c(prop=as.numeric(hab_prop), var=as.numeric(hab_var), het=as.numeric(hab_het))
}
