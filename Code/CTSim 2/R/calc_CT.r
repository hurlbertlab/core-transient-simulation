#' Calculate core/transient richness and abundance
#' 
#' Calculates richness or abundance of core and transient species through 
#' time for a simulation 
#'
#' These functions calculate species richness and abundance of
#' species in different categories of temporal occupancy.
#' The spatial and temporal scope is determined by an array of
#' species abundance profiles (\code{abuns}, as returned by
#' \code{\link{calc_abun_profile}} while the matrix of temporal occupancies 
#' used to categorize species are given in 
#' \code{occupancy} (see \code{\link{calc_occupancy}}). 
#' Species are assigned to one of two or more categories defined by
#' \code{breaks}, which is an ordered list of numbers in 
#' \code{[0,1]} defining occupancy category boundaries. Occupancy category 
#' intervals include their upper bound, but not their lower bound.
#' Thus, species with occupancy = \code{0} are not counted.
#' Specifying a large number of breaks may be useful for generating
#' occupancy distributions.
#' Timepoints may be aggregated prior to calculations using 
#' \code{agg_times}. This can be useful for evaluating the effects 
#' of sampling at different temporal scales on observed richness
#' and abundance. This aggregation is relative to the timepoints
#' stored in \code{abuns} and not to the actual 
#' timesteps in the simulation. For example, specifying 
#' \code{agg_times=10} on abundance profiles containing
#' every other timestep will aggregate abundances from \code{T=1:20,
#' T=21:40,...}.
#'
#' @param abuns (required) array of abundance profiles (as returned by 
#' 	\code{\link{calc_abun_profile}})
#' @param occupancy (required) matrix of species temporal occupancies at
#' 	a set of locations (as returned by \code{\link{calc_occupancy}})
#' @param breaks (required) a vector of numbers on \code{[0,1]} (in order) 
#' 	specifying boundaries of occupancy categories
#' @param agg_times either a single number specifing the number of timepoints 
#' 	that should be aggregated before calculating occupancy or a list of
#' 	vectors defining exactly which timepoints should be aggregated (see details).
#' 	 Defaults to no aggregation.
#' @return an array with dimensions \code{[temporal units, occupancy 
#' 	categories, spatial units]} 

#' @describeIn calc_CT Calculate core-transient species richness
#' @export
calc_rich_CT = function(abuns, occupancy, breaks, agg_times=NULL){

	# Catch error when abuns and occupancy do not match
	if(nrow(occupancy)!=dim(abuns)[length(dim(abuns))]) stop('Abundance matrix and occupancy matrix must be supplied for same set of sites.')
	
	# Convert occupancy to factor based on breakpoints
	if(breaks[1]!=0) breaks = c(0,breaks)
	if(breaks[length(breaks)]!=1) breaks = c(breaks, 1)
	cats = matrix(cut(occupancy, breaks, include.lowest=F, labels=F), nrow=nrow(occupancy), ncol=ncol(occupancy))
	colnames(cats) = colnames(occupancy)
		
	# If agg_times undefined, defaults to each timestep
	if(is.null(agg_times)) agg_times = 1

	# Number of sites and time points
	N_c = dim(abuns)[3]
	N_t = dim(abuns)[1]

	# Determine time aggregation windows
	if(is.numeric(agg_times)){
		if(length(agg_times)>1) stop('To specify non-uniform aggregation use list format.')
		int = agg_times
		agg_times = lapply(seq(1,N_t,int), function(x) x:(x+int-1))
	}
	
	# Names of aggregated times
	agg_time_names = sapply(agg_times, function(x) paste(dimnames(abuns)[[1]][x], collapse='-'))

	# Number of time windows
	N_tw = length(agg_times)

	# Richness of each category for each site during given aggregated time periods. Returns [times, categories, sites]
	rich_cat = sapply(1:N_c, function(i){
		sapply(1:(length(breaks)-1), function(cat){
			these_sp = which(cats[i,]==cat)
			if(length(these_sp)>0){
				calc_rich(abuns = abuns[,,i], agg_times = agg_times, which_species = these_sp)
			} else {
				rep(0, N_tw)
			}
		})
	}, simplify='array')
	dimnames(rich_cat) = list(agg_time_names, levels(cut(0, breaks)), 1:N_c)

	# Return richnes
	rich_cat
}

#' @describeIn calc_CT Calculate core-transient abundance
#' @export
calc_abun_CT = function(abuns, occupancy, breaks, agg_times=NULL){

	# Catch error when abuns and occupancy do not match
	if(nrow(occupancy)!=dim(abuns)[length(dim(abuns))]) stop('Abundance matrix and occupancy matrix must be supplied for same set of sites.')
	
	# Convert occupancy to factor based on breakpoints
	if(breaks[1]!=0) breaks = c(0,breaks)
	if(breaks[length(breaks)]!=1) breaks = c(breaks, 1)
	cats = matrix(cut(occupancy, breaks, include.lowest=F, labels=F), nrow=nrow(occupancy), ncol=ncol(occupancy))
	colnames(cats) = colnames(occupancy)
	
	# If agg_times undefined, defaults to each timestep
	if(is.null(agg_times)) agg_times = 1

	# Number of sites and time points
	N_c = dim(abuns)[3]
	N_t = dim(abuns)[1]

	# Determine time aggregation windows
	if(is.numeric(agg_times)){
		if(length(agg_times)>1) stop('To specify non-uniform aggregation use list format.')
		int = agg_times
		agg_times = lapply(seq(1,N_t,int), function(x) x:(x+int-1))
	}
	
	# Names of aggregated times
	agg_time_names = sapply(agg_times, function(x) paste(dimnames(abuns)[[1]][x], collapse='-'))

	# Number of time windows
	N_tw = length(agg_times)

	# Abundance of each category for each site during given aggregated time periods. Returns [times, categories, sites]
	abun_cat = sapply(1:N_c, function(i){
		sapply(1:(length(breaks)-1), function(cat){
			these_sp = which(cats[i,]==cat)
			if(length(these_sp)>0){
				sapply(agg_times, function(j){
					sum(abuns[j,as.character(these_sp),i])
				})
			} else {
				rep(0, N_tw)
			}
		})
	}, simplify='array')
	dimnames(abun_cat) = list(agg_time_names, levels(cut(0, breaks)), 1:N_c)

	# Return richnes
	abun_cat
}
