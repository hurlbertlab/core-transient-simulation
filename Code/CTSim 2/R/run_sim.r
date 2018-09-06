#' Run simulation
#'
#' Runs the simulation on a metacommunity for a fixed number of timesteps.
#'
#' Given a metacommunity (\code{metacomm}, see 
#' \code{\link{populate_landscape}}), landscape (\code{land}, 
#' see \code{\link{make_landscape}}), species pool (\code{species}, see
#' \code{\link{make_species}}), and global species abundance distribution
#' (\code{gsad}), this function calls \code{\link{run_timestep}}
#' to progress a metacommunity through a fixed number of timesteps (given by 
#' \code{steps}). The user may specify which timesteps should be returned using the
#' \code{save_steps} parameter. Including \code{0} in \code{save_steps} will  
#' save the initial metacommunity. The function can also be told to calculate
#' species turnover rates (by specifying \code{calc_rates=TRUE}), which will
#' calculate the frequency at which each species establishes in or vacates a 
#' microsite (through death or movement). 
#' 
#' @param steps (required) number of timesteps to run the simulation
#' @inheritParams run_timestep
#' @param save_steps vector of timesteps at which to save the simulation. 
#' 	Defaults to all timesteps, including the initial metacommunity.
#' @param report integer giving the interval at which the function should 
#' 	report its progress to \code{stdout}. Defaults to 0 for no reporting.
#' @param ID character string identifying this simulation run when reporting progress
#' @param calc_rates logical indicating whether the function should calculate
#' 	rates of turnover for each species based on how frequently a species vacates
#' 	or establishes in a microsite. Defaults to \code{FALSE}.
#' @return if \code{calc_rates} is\code{FALSE} then the function returns 
#' a 3-dimensional array of lists defining the metacommunity at each 
#' 	timepoint in \code{save_steps} where time is in the 3rd dimension.
#' 	Otherwise, the function returns a list where the first element (\code{'sim'}
#' 	is the array recording the metacommunity through time and the second
#' 	element \code{'turnover'} is an array counting the number of species losses
#' 	and gains in each cell across time intervals specified in \code{save_steps}.
#' 	Dimensions of \code{turnover} are \code{[time step, row, column, species, 
#' 	gain/loss rate]}.
#'
#' @seealso \code{\link{run_timestep}} for how the simulation progresses 
#' 	through a single timestep. \cr
#' \code{\link{run_sim_N}} for running multiple independent simulations on 
#' 	the same set of parameters.
#' @export

run_sim = function(steps, metacomm, land, species, gsad, d_kernel=NULL, v_kernel=NULL, imm_rate=NA, save_steps = NULL, report=0, ID=NA, calc_rates=F){
	
	# Define simulation
	X = nrow(land)
	Y = ncol(land)
	N_S = dim(species)[1]

	# Create array to save simulations results
	if(is.null(save_steps)) save_steps = 0:steps
	sim_results = array(list(), dim=c(X, Y, length(save_steps)),
		dimnames=list(row=1:X, col=1:X, time=save_steps))

	# Create array for saving species turnover rates
	if(calc_rates){
		turnover = array(0, dim=c(X, Y, steps, N_S, 2), 
			dimnames=list(row=1:X, col=1:Y, time=1:steps, species=1:N_S, c('loss','gain'))
		)
	}
	
	# Save initial step
	if(0 %in% save_steps) sim_results[,,'0'] = metacomm

	# Run simulation
	new_metacomm = metacomm
	for(step in 1:steps){
		# Save old version of metacommunity
		old_metacomm = new_metacomm
	
		# Run for a timestep
		new_metacomm = run_timestep(old_metacomm, land, species, gsad, d_kernel, v_kernel, imm_rate, return_dead=calc_rates)

		# Calculate species turnover
		if(calc_rates){
			for(i in 1:X){
			for(j in 1:Y){
				
				# Calculate how frequently each species was lost
				sp_lost = old_metacomm[i,j][[1]][new_metacomm$comm_loss[i,j][[1]]]
				turnover[i,j,step,,'loss'] = table(factor(sp_lost, levels=1:N_S))
				
				# Calculate how frequently each species established
				empty_spots = old_metacomm[i,j][[1]]==0 | new_metacomm$comm_loss[i,j][[1]]
				sp_gain = new_metacomm$newcomm[i,j][[1]][empty_spots]
				turnover[i,j,step,,'gain'] = table(factor(sp_gain, levels=1:N_S))
			}}
			
			# Rename new metacommunity for next timestep
			new_metacomm = new_metacomm$newcomm
		}
		
		# Save results
		if(step %in% save_steps) sim_results[,,as.character(step)] = new_metacomm
		
		# Report progress
		if(report > 0) if(step %% report == 0) print(paste(Sys.time(), ': Finished', step, 'of', steps, 'from run', ID))
	}
	
	# Return results
	if(calc_rates){
	
		# Aggregate turnover rates to the steps that should be saved
		turnover_summed = apply(turnover, c(1,2,4,5), function(x){
			sums = sapply(save_steps[save_steps!=0], function(tp) sum(x[1:tp]))
			sums[2:length(sums)] = sums[2:length(sums)] - sums[1:(length(sums)-1)]
			sums
		})
		dimnames(turnover_summed)[[1]] = save_steps[save_steps!=0]
		names(dimnames(turnover_summed))[1] = 'time'
		
		list(sim=sim_results, turnover=turnover_summed)
	} else {
		sim_results
	}
}
