#' Species turnover
#'
#' Calculate rates of species' gains and losses
#'
#' DOCUMENTATION NOT DONE
#' The function calculates the number of individuals in a group of 
#' species that are gained or lost per timestep. The function takes 
#' a turnover matrix, returned by \code{\link{run_sim}} when
#' code{calc_rates=TRUE} and then calculates the number of individuals
#' of each species that are gained or lost per timestep, summing across
#' cells in each spatial unit in \code{locs}. A parameter \code{agg_times}
#' can be used to aggregate timepoints prior to calculating rates. 
#' This aggregation is relative to the timepoints stored in the turnover
#' array and not to the actual timesteps in the simulation. 
#' For example, specifying \code{agg_times=10} on a turnover matrix containing
#' every other timestep will aggregate individuals from \code{T=1:20,
#' T=21:40,...} and then divide each by 20 timesteps to calculate the rates.
#' Once individual species rates have been calculated, the function sums
#' rates across species groups (\code{which_species}), given as a vector 
#' that can be interpreted as a factor identifying which species 
#' belong to which groups. If the user wishes to calculate
#' average rates within groups they will need to divide the results by
#' the number of species in each group. By default species are not grouped and 
#' rates for individual species are returned. 
#' 
#' @param turnover (required) an array of species rates of turnover, as
#' 	returned by \code{run_sim} when \code{calc_rates=TRUE}. Should have
#' dimensions \code{[time, row, column, species, gain/loss rate]}.
#' @param locs (required) either a two-column matrix of cell locations where 
#' 	species' rates should be summed or a list of matrices designating
#' 	multiple spatial units where rates should be summed (see details)
#' @param t_window (required) either a list containing \code{start}
#' 	and \code{stop} specifying that all collected timepoints in that interval
#' 	should be considered or an explicit vector of timepoints 
#' @param agg_times either a single number specifing the number of timepoints 
#' 	that should be aggregated before calculating rates or a list of
#' 	vectors defining exactly which timepoints should be aggregated (see details).
#' 	Defaults to no aggregation.
#' @param which_species a vector of integers indicating which species should be
#' 	grouped and averaged. Defaults to each species analyzed separately.
#' @return a array of rates of species gains and losses in each location at 
#' 	each timepoint (dimensions \code{[temporal unit, spatial unit, 
#' 	gain or loss, species group]}). 
#'
#' @import abind
#' @export

calc_species_turnover = function(turnover, locs, t_window, agg_times=NULL, which_species=NULL){

	# Determine which timepoints to evaluate
	timepoints = as.numeric(dimnames(turnover)[[1]])
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
	
	# Stop if no times
	if(length(use_times)==0) stop('Cannot calculate abundance on 0 timepoints. Check that t_window is in sim.')

	# Convert to actual rates (individuals per timestep) by calculating how many timesteps are included for each recorded timepoint
	ntimes = timepoints; names(ntimes) = dimnames(turnover)[[1]]
	ntimes[2:length(ntimes)] = ntimes[2:length(ntimes)] - ntimes[1:(length(ntimes)-1)]
	turnover_rates = apply(turnover, 2:5, function(x) x/ntimes)
	
	# Determine whether cells should be aggregated
	if(!is.list(locs)) locs = list(locs)

	# Aggregate cells
	turnover_summed = sapply(locs, function(cell_block){
	
		# For each block of cells to be aggregated
		each_cell = sapply(1:nrow(cell_block), function(i){
			
			# Convert community profile to matrix
			x = as.numeric(cell_block[i,])
			turn_mat = turnover_rates[as.character(use_times), x[1], x[2],,, drop=F]

		}, simplify='array')

		# Sum cells across species and timepoints
		summed_cells = apply(each_cell, 1:5, sum)
		
		abind::adrop(summed_cells, 2:3)
		# dims are now [timepoints, species, gain/loss]
	}, simplify='array')
	# dims are now [timepoints, species, gain/loss, spatial unit]
	names(dimnames(turnover_summed))[3:4] = c('rate','space')
	
	# Reduce ntimes to the focal timepoints
	ntimes = ntimes[as.character(use_times)]

	# Number of timepoints and sites measured
	N_t = dim(turnover_summed)[1]
	N_c = dim(turnover_summed)[4]

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
		agg_rates = sapply(agg_times, function(rows){
			rows = rows[rows %in% 1:N_t] 
			use_rates = turnover_summed[rows,,,,drop=F]
			if(length(rows)>1){
				new_rates = apply(use_rates, 2:4, function(x){
					# Must re-weight rates by number of timesteps included in these timepoints
					sum(x*ntimes[rows])/sum(ntimes[rows]) 
				})
				new_rates = add_dim(new_rates, 1)
			}
			if(length(rows)==1) new_rates = use_rates
			new_rates
		}, simplify='array')
		
		# Fix dimensions
		agg_rates = abind::adrop(agg_rates, drop=1)
		if(N_c==1 & length(agg_times)>1 & length(dim(agg_rates))==3) agg_rates = add_dim(agg_rates, 3) # add dimension for spatial units
		

		# Transpose and rename dimensions
		agg_rates = aperm(agg_rates, c(4,3,2,1))
		names(dimnames(agg_rates)) = c('time','space','rate','species')
		dimnames(agg_rates)[[1]] = paste('agg', 1:length(agg_times), sep='_')
		dimnames(agg_rates)[[2]] = 1:N_c

		return_results = agg_rates
	} else {
		return_results = aperm(turnover_summed, c(1, 4, 3, 2))
		dimnames(return_results)[[2]] = 1:N_c
	}
	
	# Calculate summed rates across species groups
	if(is.null(which_species)) which_species = dimnames(turnover)[[4]]
	
	return_results = sapply(levels(factor(which_species)), function(sp){
		apply(return_results[,,,which_species==sp,drop=F], c(1:3), sum)	
	}, simplify='array')

	# Return summed rates
	return_results
}
