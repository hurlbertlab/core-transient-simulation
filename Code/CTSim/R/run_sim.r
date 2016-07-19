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
#' save the initial metacommunity.
#' 
#' @param steps (required) number of timesteps to run the simulation
#' @inheritParams run_timestep
#' @param save_steps vector of timesteps at which to save the simulation. 
#' 	Defaults to all timesteps, including the initial metacommunity.
#' @param report integer giving the interval at which the function should 
#' 	report its progress to \code{stdout}. Defaults to 0 for no reporting.
#' @param ID character string identifying this simulation run when reporting progress
#' @return a 3-dimensional array of lists defining the metacommunity at each 
#' 	timepoint in \code{save_steps}. Time is in the 3rd dimension.
#'
#' @seealso \code{\link{run_timestep}} for how the simulation progresses 
#' 	through a single timestep. \cr
#' \code{\link{run_sim_N}} for running multiple independent simulations on 
#' 	the same set of parameters.
#' @export

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
