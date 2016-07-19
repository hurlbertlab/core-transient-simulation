#' Richness profiles
#'
#' Calculate species richness through time for a simulation
#'
#' The function calculates species richness for a set of locations
#' and timepoints over a given time interval. Users can either supply an 
#' array of abundance profiles (as returned by 
#' \code{\link{calc_abun_profile}}) or a set of spatial units
#' (\code{locs}), a time window (\code{t_window}), simulation
#' results (\code{sim}) and the number of species (\code{N_S}), 
#' from which abundance profiles can be calculated. These parameters
#' are explained in more detail in the documentation of 
#' \code{\link{calc_abun_profile}}. For each spatial and temporal unit in
#' the abundance profiles, the function calculates the number of species
#' present. Timepoints may be aggregated prior to this calculation using 
#' \code{agg_times}. This can be useful for evaluating the effects 
#' of sampling at different temporal scales on observed richness.
#' This aggregation is relative to the timepoints
#' stored in the abundance profiles array and not to the actual 
#' timesteps in the simulation. For example, specifying 
#' \code{agg_times=10} on abundance profiles containing
#' every other timestep will aggregate abundances from \code{T=1:20,
#' T=21:40,...}.
#' 
#' @inheritParams calc_abun_profile
#' @param abuns array of abundance profiles (as returned by 
#' 	\code{\link{calc_abun_profile}})
#' @param agg_times either a single number specifing the number of timepoints 
#' 	that should be aggregated before calculating occupancy or a list of
#' 	vectors defining exactly which timepoints should be aggregated (see details).
#' 	 Defaults to no aggregation.
#' @param which_species a vector of integers indicating which species should be
#' 	considered. Defaults to all species.
#' @return a matrix of species richness in each location at each timepoint 
#' (dimensions \code{[temporal unit, spatial unit]}) 
#'
#' @import abind
#' @export

calc_rich = function(locs=NULL, t_window=NULL, sim=NULL, N_S=NULL, abuns=NULL, agg_times=NULL, which_species=NULL){
	
	# Catch error where not enough information is specified
	if(is.null(abuns)&(is.null(locs)|is.null(t_window)|is.null(sim)|is.null(N_S))) stop('Must specify abuns or locs/t_window/sim/N_S.')	

	# Calculate species abundance profiles if not specified
	if(is.null(abuns)){
		abun_profiles = calc_abun_profile(locs, t_window, sim, N_S)
	} else {
		abun_profiles = abuns
	}

	# Convert to array if necessary (e.g., if this is a matrix of abundances from a single site)
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
				new_abun = add_dim(new_abun, 1)
			}
			if(length(rows)==1) new_abun = use_abun
			new_abun
		}, simplify='array')
		
		# Fix dimensions
		agg_abuns = abind::adrop(agg_abuns, drop=1)
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
