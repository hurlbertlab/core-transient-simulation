#' Summarize a simulation
#' 
#' Summarizes the output of a single simulation run.
#'
#' The function calculates summary statistics of four main quantities and 
#' returns them in a named list:
#' \describe{
#'	\item{bio}{richness and abundance of biologically core and transient 
#'		species (based on birth rates)}
#'	\item{occ}{richness and abundance of species in classes based on
#'		their temporal occupancy}
#'	\item{xclass}{number of species in each of four categories: biologically
#'		core - occupancy core, biologically core - occupancy transient,
#'		biologically transient - occupancy core, and biologically transient -
#'		occupancy transient (see\code{\link{cross_classify}})}
#'	\item{abun}{relative rank abundance of each species
#'		in the metacommunity}
#' }
#'
#' The function proceeds as follows:
#' \enumerate{
#'	\item Actual species' abundances are calculated for each spatial unit 
#'		in \code{locs} and temporal unit in \code{t_window}. See documentation
#'		of \code{\link{calc_abun_profile}} for more information on 
#'		these parameters.
#'	\item Observed species' abundance profiles are calculated by sampling
#'		actual abundance profiles using the detection probabilities 
#'		provided in \code{P_obs}.
#'	\item Observed abundance profiles are used to calculate species 
#'		temporal occupancy in each of the spatial units, after aggregating 
#'		timepoints according to \code{agg_times} (see 
#'		\code{\link{calc_occupancy}}).
#'	\item Richness and abundance of species in each of the temporal occupancy
#'		categories defined by \code{breaks} are calculated in each spatial 
#'		unit for each timepoint. Timepoints may
#'		be aggregated prior to calculations using the \code{agg_times}
#'		argument in \code{sum_parms}. See (see \code{\link{calc_rich_CT}}) for 
#'		details.
#'	\item For each spatial unit, habitat values are averaged across the
#'		cells that comprise it (see \code{\link{average_habitat}} and this 
#'		average habitat type is used to determine which species are 
#'		biologically core and transient (based on whether they have positive
#'		birth rates in that habitat). Richness and abundance of biologically
#'		core and transient species is then calculated for each spatial unit
#'		at each timepoint. As in the previous step, timepoints may be 
#'		aggregated prior to calculations using the \code{agg_times} 
#'		argument in \code{sum_parms}.
#'	\item For each spatial unit, species are classified as core or transient 
#'		based on whether they fall in the first or last occupancy interval defined
#'		by \code{breaks}. Then, the number of species classified as biologically
#'		core - occupancy core, biologically core - occupancy transient,
#'		biologically transient - occupancy core, and biologically transient -
#'		occupancy transient are counted for each spatial unit.
#'	\item The element \code{time_sum} in \code{sum_parms} defines which 
#'		timepoints are used for summarizing richness and abundance. If 
#'		\code{time_sum = 'none'} then all timepoints are summarized 
#'		individually. Other options are \code{'mean'} for averaging across all
#'		timepoints in \code{t_window} or \code{'last'} for only using the
#'		last timepoint in \code{t_window}.
#'	\item Summary statistics are calculated across all spatial units.
#'		Minimally, this is the mean and variance, but may also include any
#'		quantiles specified in the \code{quants} argument of \code{sum_parms}.
#'	\item The relative abundance of each species is calculated across the entire
#'		metacommunity for each set of observed abundances (from step 2).
#' }
#'
#' @param sim (required) either a filename for a simulation run (as saved by
#' 	\code{\link{run_sim_N}} or an array of simulation results (as returned by 
#' 	\code{\link{run_sim}}. If an array, must specify \code{species}, 
#' 	\code{land} and \code{gsad}.
#' @param breaks (required) a vector of numbers on \code{[0,1]} (in order) 
#' 	specifying boundaries of occupancy categories
#' @param locs (required) either a two-column matrix of cell locations where 
#' 	species' abundances should be summed or a list of matrices designating
#' 	multiple spatial units where abundance should be summed (see 
#' 	\code{\link{calc_abun_profile}})
#' @param t_window (required) either a list containing \code{start}
#' 	and \code{stop} specifying that all collected timepoints in that interval
#' 	should be considered or an explicit vector of timepoints
#' @param species array of species vital rates used in the simualation
#' 	(as generated by \code{\link{make_species}})
#' @param land matrix or raster of habitat types defining
#' 	the landscape used in the simulation
#' 	(as generated by \code{\link{make_landscape}})
#' @param gsad vector defining the global relative abundance of each species
#' 	used in the simulation. Must be in the same order as \code{species}. 
#' @param agg_times either a single number specifing the number of timepoints 
#' 	that should be aggregated before calculating occupancy or a list of
#' 	vectors defining exactly which timepoints should be aggregated (see details).
#' 	 Defaults to no aggregation.
#' @param P_obs vector of detection probabilities (see details and 
#' 	\code{\link{sample_sim}})
#' @param sum_parms list of parameters controlling how the simulation is 
#' 	summarized across spatial and temporal units. May contain:
#' 	\describe{
#' 		\item{agg_times}{specifies how time points should be aggregated 
#'			before calculations. See \code{\link{calc_rich_CT}} for details.}
#'		\item{time_sum}{character vector indicating which timepoint
#'			should be used in summary statistics (see details)}
#'		\item{quants}{vector of quantiles for summarizing across spatial units
#'			(see details)}
#'		\item{hab}{character string indicating whether summary should be performed
#'			across all habitat types (default) or for just one (indicated by 
#'			\code{'A'} or \code{'B'}.}
#'	}
#' @param sum_turn logical indicating whether turnover rates should be
#'  summarized. Defaults to \code{TRUE}.
#' @return a list of arrays defined as follows:
#' 	\describe{
#' 		\item{bio}{richness and abundance of biologically defined species,
#' 			has dimensions: \code{[richness or abundance, spatial 
#' 			summary statistic, core or transient, detection probability]}}
#' 		\item{occ}{richness and abundance in temporal occupancy categories,
#' 			has dimensions: \code{[richness or abundance, spatial 
#' 			summary statistic, occupancy category, detection probability]}}
#' 		\item{xclass}{core-transient cross-classification,
#' 			has dimensions: \code{[spatial summary statistic,
#' 			 cross-classification category, detection probability]}}
#' 		\item{abun}{species' relative abundances, has dimensions:
#' 			\code{[species, detection probability]}}
#' 	}
#' 	If \code{sum_parms$time_sum = 'none'} then each of these arrays will
#' 	have an additional dimension corresponding to timepoints.
#'
#' @import reshape2
#' @import abind
#' @export

summarize_sim = function(sim, breaks, locs, t_window, species=NULL, land=NULL, gsad=NULL, agg_times=NULL, P_obs=list(1), sum_parms=NULL, sum_turn=TRUE){
	
	# If sim is a file, then read in simulation run. Should have objects: results, this_land, this_species, this_gsad
	if(is.character(sim)){
		load(sim)
		
		# Check that required objects exist
		if(!exists('this_species')|!exists('this_land')|!exists('this_gsad')){
			stop(paste(sim, 'may not contain this_species, this_land, or this_gsad. Summary NOT run.'))		
		} else {
			# Assign objects
			species = this_species
			land = this_land
			gsad = this_gsad
		}
	} else {
		# If sim is an array of simulation results, then must specify species, land, gsad
		if(is.null(species)|is.null(land)|is.null(gsad)) stop('If sim is not a file, must supply species, land, and gsad.')
		results = sim
	}
	
	# If results is a list with two elements, then extract the elements corresponding to the metacommunity and species turnover
	if(length(results)==2){
		turnover = results$turnover
		results = results$sim
	}

	# Number of species
	N_S = dim(species)[1]

	# Calculate species abundance profiles at the spatial and temporal resolution given by locs and t_window
	abuns_act = calc_abun_profile(locs, t_window, results, N_S)

	# Apply observation bias
	abuns_obs = sapply(P_obs, function(p){
		sample_sim(abuns_act, probs = p, return='abundance')
	}, simplify='array')
	dimnames(abuns_obs)[[2]] = 1:N_S # Name columns with species names
	# dims are now [timepoint, species, spatial unit, P]

	# Calculate species occupancy
	occ = sapply(1:length(P_obs), function(i){
		use_abun = abuns_obs[,,,i, drop=FALSE]
		use_abun = abind::adrop(use_abun, 4) 
		calc_occupancy(abuns=use_abun, agg_times, do_freq=F)
	}, simplify='array')
	# dims are now: [spatial unit, species, P]

	# Calculate richness and in each occupancy category
	rich_ct = sapply(1:length(P_obs), function(i){
		use_abun = abuns_obs[,,,i, drop=FALSE]
		use_abun = abind::adrop(use_abun, 4) 
		use_occ = occ[,,i, drop=FALSE]
		use_occ = abind::adrop(use_occ, 3)
		calc_rich_CT(use_abun, use_occ, breaks, agg_times=sum_parms$agg_times)
	}, simplify='array')
	abun_ct = sapply(1:length(P_obs), function(i){
		use_abun = abuns_obs[,,,i, drop=FALSE]
		use_abun = abind::adrop(use_abun, 4) 
		use_occ = occ[,,i, drop=FALSE]
		use_occ = abind::adrop(use_occ, 3)
		calc_abun_CT(use_abun, use_occ, breaks, agg_times=sum_parms$agg_times)
	}, simplify='array')
	# dims are now [timepoint, category, spatial unit, P]

	# Get species birth rates
	b_rates = species[,,'b']
	
	# Get habitat types for spatial units
	habitats = sapply(locs, function(x) average_habitat(x, land))

	# Calculate classification of species based on birth rates
	cores = t(sapply(habitats, function(h) b_rates[,h]>0))
	classification = apply(cores, 1:2, function(x) ifelse(x, 'core', 'trans'))	

	# Calculate richness and abundnace of biologically core vs transient species
	occ_ab = apply(classification, 1:2, function(x) ifelse(x=='core', 1, .1))
	rich_ab = sapply(1:length(P_obs), function(i){
		use_abun = abuns_obs[,,,i, drop=FALSE]
		use_abun = abind::adrop(use_abun, 4) 
		calc_rich_CT(use_abun, occ_ab, 0.5, agg_times=sum_parms$agg_times)
	}, simplify='array')
	abun_ab = sapply(1:length(P_obs), function(i){
		use_abun = abuns_obs[,,,i, drop=FALSE]
		use_abun = abind::adrop(use_abun, 4) 
		calc_abun_CT(use_abun, occ_ab, 0.5, agg_times=sum_parms$agg_times)
	}, simplify='array')
	# dims are now [timepoint, category, spatial unit, P]

	# Calculate proportion mis-classified
	xclass = sapply(1:length(P_obs), function(i){
		use_occ = occ[,,i, drop=FALSE]
		use_occ = abind::adrop(use_occ, 3)
		tabs = cross_classify(use_occ, breaks, classification=classification, do_each=T, return='counts')
		reshape2::acast(reshape2::melt(tabs, varnames=c('bio','occ','sp_unit')), sp_unit ~ bio+occ)
	}, simplify='array')
	# dims are now [spatial unit, category, P]
	
	# Calculate species' relative abundances across the metacommunity
	abuns_global = sapply(1:length(P_obs), function(i){
		use_abun = abuns_obs[,,,i, drop=FALSE]
		use_abun = abind::adrop(use_abun, 4) 
		apply(use_abun, 1, function(x){
			# Sum across rows and order 
			rank_abun = rowSums(x)/sum(x)
			rank_abun = rank_abun[order(rank_abun, decreasing=T)]
		})
	}, simplify='array')
	# dims are now [species rank, timepoint, P]
	
	# Calculate species turnover rates, if present
	if(exists('turnover')&sum_turn){
		# Calculate for biologically core/transient
		# Currently does not aggregate time periods
		turn_ab = sapply(1:length(locs), function(i){
			calc_species_turnover(turnover, locs[i], t_window, which_species=classification[i,])
		}, simplify='array')
		# If sapply returns a list because core or transient species are missing from a spatial unit
		if(is.list(turn_ab)){
			turn_ab = sapply(turn_ab, function(x){
				dim_missing = c('core','trans')[!c('core','trans') %in% dimnames(x)[[4]]]
				if(length(dim_missing)==1){
					new_x = abind::abind(x, array(NA, dim=dim(x)[1:3]), along=4)
					dimnames(new_x)[[4]][2] = dim_missing
					new_x = new_x[,,,c('core','trans')]
				} else {
					new_x = x[,,,c('core','trans')]
				}
				new_x
			}, simplify='array')
		}
		if(length(dim(turn_ab))==5) turn_ab = abind::adrop(turn_ab, 2)
		turn_ab = aperm(turn_ab, c(1,3,4,2))
		# dimension are [time, category, spatial unit, rate type]
		
		# Calculate for each temporal occupancy class
		# Convert occupancy to factor based on breakpoints
		use_occ = calc_occupancy(abuns=abuns_act, agg_times, do_freq=F)
		if(breaks[1]!=0) breaks = c(0,breaks)
		if(breaks[length(breaks)]!=1) breaks = c(breaks, 1)
		cats = matrix(cut(use_occ, breaks, include.lowest=F, labels=F), nrow=nrow(use_occ), ncol=ncol(use_occ))

		turn_ct = sapply(1:length(locs), function(i){
			calc_species_turnover(turnover, locs[i], t_window, agg_times, cats[i,])
		}, simplify='array')
		
		## ACTUALLY DECIDED NOT TO OUTPUT THIS B/C INFO NOT NEEDED.
		## ALSO, TURN_CT IS USUALLY A LIST BECAUSE SOME OCCUPANCY CLASSES ARE USUALLY MISSING FROM MOST LOCATIONS
		## WOULD NEED TO CONVERT BACK FROM LIST TO USE BY FILLING IN MISSING ARRAY DIMENSIONS
	}

	# Determine which time window to use for summary
	# Defaults to mean
	if(is.null(sum_parms$time_sum)){time_sum = 'mean'} else { time_sum = sum_parms$time_sum }

	# Average across time windows
	if(time_sum=='mean'){
		rich_ab = apply(rich_ab, 2:4, mean)
		rich_ct = apply(rich_ct, 2:4, mean)
		abun_ab = apply(abun_ab, 2:4, mean)
		abun_ct = apply(abun_ct, 2:4, mean)
		abuns_global = apply(abuns_global,c(1,3), mean) 
		if(exists('turn_ab')) turn_ab = apply(turn_ab, 2:4, mean)

		# Define dimension that stores spatial units
		sp_dim = 2
	}
	
	# Only use the last time window
	if(time_sum=='last'){
		N_t = dim(rich_ab)[1]
		rich_ab = rich_ab[N_t,,,,drop=F]; rich_ab = abind::adrop(rich_ab, 1)
		rich_ct = rich_ct[N_t,,,,drop=F]; rich_ct = abind::adrop(rich_ct, 1)
		abun_ab = abun_ab[N_t,,,,drop=F]; abun_ab = abind::adrop(abun_ab, 1)
		abun_ct = abun_ct[N_t,,,,drop=F]; abun_ct = abind::adrop(abun_ct, 1)
		abuns_global = abuns_global[,N_t,,drop=F]; abuns_global = abind::adrop(abuns_global, 2)
		if(exists('turn_ab')){
			turn_ab = turn_ab[N_t,,,,drop=F]; turn_ab = abind::adrop(turn_ab, 1)
		}
	
		# Define dimension that stores spatial units
		sp_dim = 2
	}

	# Do not summarize across time windows
	if(time_sum=='none'){
		# Define dimension that stores spatial units
		sp_dim = 3
	}

	# Catch error if time_sum given does not match an available method.
	if(!(time_sum %in% c('none','last','mean'))) stop('Incorrect temporal summary method given in sum_parms.')
	
	# Calculate means and variances across spatial units
	# dim is the dimension where spatial units are stored
	# use_obs is a logical vector indicating which observations should be used 
	calc_stats = function(dat, sp_dim, quants=sum_parms$quants, use_obs=NULL){
		if(is.null(use_obs)) use_obs = rep(T, dim(dat)[sp_dim])
	
		keep_dims = 1:length(dim(dat))
		keep_dims = keep_dims[keep_dims!=sp_dim]
		
		apply(dat, keep_dims, function(x){
			stats = c(mean=mean(x[use_obs]), var=var(x[use_obs]))
			if(!is.null(quants)) stats = c(stats, quantile(x[use_obs], quants))
			stats
		})
	}

	# Determine which spatial units to used (based on habitat type)
	if(is.null(sum_parms$hab)){ hab='AB' } else {hab=sum_parms$hab}
	
	# If the summary is only for one habitat type
	if(hab!='AB'){
		# Check whether landscape is supplied
		if(is.null(land)) stop('Must supply landscape in order to summarize for only one habitat type.')
	
		# Determine the habitat type of the landscape
		land_types = aggregate_hab_type(land, locs)
		
		# Subset data to only the habitat type of interest
		sp_subset = land_types==hab	
	} else {
		sp_subset = rep(T, length(locs))
	}

	if(time_sum=='none'){
		bio_names = list(c('rich','abun'), NULL, NULL, c('trans','core'), P_obs)
		occ_names = list(c('rich','abun'), NULL, NULL, NULL, P_obs)
	} else{
		bio_names = list(c('rich','abun'), NULL, c('trans','core'), P_obs)
		occ_names = list(c('rich','abun'), NULL, NULL, P_obs)
	}
	
	bio_stats = abind::abind(calc_stats(rich_ab, sp_dim, use_obs=sp_subset), calc_stats(abun_ab, sp_dim, use_obs=sp_subset), along=0, new.names=bio_names)
	occ_stats = abind::abind(calc_stats(rich_ct, sp_dim, use_obs=sp_subset), calc_stats(abun_ct, sp_dim, use_obs=sp_subset), along=0, new.names=occ_names)
	xclass_stats = calc_stats(xclass, 1, use_obs=sp_subset)
	if(exists('turn_ab')) turn_stats = calc_stats(turn_ab, sp_dim, use_obs=sp_subset)
	dimnames(xclass_stats)[[3]] = P_obs
	dimnames(abuns_global)[[sp_dim]] = P_obs

	# Return list of statistics
	return_list = list(bio=bio_stats, occ=occ_stats, xclass=xclass_stats, abuns=abuns_global)
	if(exists('turn_stats')) return_list = c(return_list, turnover=list(turn_stats))
	return_list
}