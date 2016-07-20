#' Summarize with multiple parameter sets
#'
#' Summarize simulation runs with multiple sets of parameters 
#'
#' This function summarizes multiple simulation runs with multiple
#' sets of parameters, which are saved in the directory
#' given in \code{parm_dir}. Parameter filenames must start with 's_'.
#' Each parameter file should have a unique \code{sumID} defined within it
#' as well define objects that can be passed as parameters to 
#' \code{\link{summarize_sim_N}}. Parameters
#' requiring a value are \code{breaks}, \code{locs} and \code{t_window}.
#' The function will either summarize simulation results for a set time period
#' defined in \code{t_window} (default) or for multiple consecutive time windows
#'  (use \code{cross_time=TRUE}), in which case \code{t_window} 
#' defines the time interval and must be a list with named elements
#' \code{start} and \code{stop}. If \code{cross_time=FALSE} then 
#' two summary objects are saved to the .RData file: \code{sim_sum_ind}- 
#' includes a summary for each run, and \code{sim_sum}- summarizes quantities 
#' across runs using the function defined in \code{sum_func}. 
#' If \code{cross_time=TRUE} then only \code{sim_sum_ind} is saved and 
#' \code{T=<time>} is appended to the filename to denote the timestep at 
#' which the summary ends.
#' Results are saved in the directory \code{results_dir} and are not returned
#' by the function. Filenames follow the convention 
#' \code{<sumID>_summary.RData}.
#'
#' @param run_dir directory where simulation runs are located.
#' 	Defaults to the current directory.
#' @param parm_dir directory where parameter files are located.
#' 	Defaults to the current directory.
#' @param results_dir directory where results should be saved.
#' 	Defaults to 'Summaries' in the current directory.
#' @param cross_time logical indicating whether to use \code{t_window}
#' 	to summarize the simulation in windows of time
#' @return nothing
#'
#' @seealso \code{\link{summarize_simulation}} for command line execution
#' @export

summarize_sim_P = function(run_dir='./', parm_dir='./', results_dir='./Summaries/', cross_time=F){
	# Create summary directory if it does not exist
	if(!file.exists(results_dir)) dir.create(results_dir)
	
	# Read in parameter files
	# Parameter files are designated by starting with 'p_'
	file_list = list.files(parm_dir, '^s_')

	# Run set of simulations on each parameter file
	for(f in file_list){
		print(paste('Started',f))
	
		# Read in parameters
		parm_file = file.path(parm_dir, f)
		source(parm_file, local=TRUE)

		# Check that required parameters are present
		if(!exists('sumID')|!exists('t_window')|!exists('breaks')|!exists('locs')){
			warning(paste(f,'may not contain sumID, t_window, breaks, or locs. Summary NOT run.'))
		} else {
		
			# Assign missing parameters to NULL
			if(!exists('P_obs')) P_obs=NULL
			if(!exists('sum_parms')) sum_parms=NULL
			if(!exists('sum_func')) sum_func=NULL
			if(!exists('agg_times')) agg_times=NULL

			if(cross_time){
				# Define time windows
				dT = t_window$stop - t_window$start
				end_times = seq(t_window$stop, dT, -dT)

				# For each end time point, summarize the simulation runs
				for(endT in end_times){
					use_twindow = list(start=endT-dT+1, stop=endT)

					# Summaries for each individual run
					sim_sum_ind = summarize_sim_N(run_dir, breaks=breaks, locs=locs, t_window=use_twindow, agg_times=agg_times, P_obs=P_obs, sum_parms=sum_parms)

					# Save
					save(sim_sum_ind, file=file.path(results_dir, paste0(sumID,'-T', endT, '_summary.RData')))

					print(paste('Finished T =', endT))
				}
			} else {
				# Summaries for each individual run
				sim_sum_ind = summarize_sim_N(run_dir, breaks=breaks, locs=locs, t_window=t_window, agg_times=agg_times, P_obs=P_obs, sum_parms=sum_parms)
				
				# Summaries across runs
				sim_sum = summarize_sim_N(run_dir, breaks=breaks, locs=locs, t_window=t_window, agg_times=agg_times, P_obs=P_obs, sum_parms=sum_parms, sum_func=sum_func)

				# Save
				save(sim_sum_ind, sim_sum, file=file.path(results_dir, paste0(sumID,'_summary.RData')))
			}
		}
		
		# Remove parameters when finished with this file
		rm(sumID, t_window, breaks, locs, P_obs, sum_parms, sum_func)
	}
}
