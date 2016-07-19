#' Run multiple parameter sets
#'
#' Run simulations on multiple sets of parameters 
#'
#' This function runs simulations on parameter files stored in the directory
#' given in \code{parm_dir}. Parameter filenames must start with 'p_'.
#' Each parameter file should have a unique \code{simID} defined within it.
#' Results are saved in the directory \code{results_dir} and are not returned
#' by the function. 
#'
#' @param ncores number of cores to use for parallel simulation.
#' 	Defaults to 1 for serial simulation.
#' @param parm_dir directory where parameter files are located.
#' 	Defaults to the current directory.
#' @param results_dir directory where results should be saved.
#' 	Defaults to 'Results' in the current directory.
#' @param sim_dir directory where \code{CTSim} is installed, 
#' 	if not on the default search path
#' @param report timestep interval at which to report simulation status
#' 	to stdout. Defaults to 0 for no reporting.
#' @param restart logical indicating whether this a restart of a previously
#' 	existing set of runs. Defaults to \code{FALSE}
#' @return nothing
#'
#' @seealso \code{\link{run_simulation}} and \code{\link{restart_simulation}}
#' for command line execution
#' @export

run_sim_P = function(ncores=1, parm_dir='./', results_dir='./Results/', sim_dir=NULL, report=0, restart=F){
	# Read in parameter files
	# Parameter files are designated by starting with 'p_'
	file_list = list.files(parm_dir, '^p_')

	# Run set of simulations on each parameter file
	for(f in file_list){

		# Read in parameters
		parm_file = paste0(parm_dir, f)
		source(parm_file)
		parm_list = make_parmlist()
		
		# Check that required parameters are present
		if(!exists('nruns')|!exists('simID')){
			warning(paste(f,'may not contain simID or nruns. Simulation NOT run.'))
		} else {

			# Run simulations
			run_sim_N(nruns, parm_list, ncores, simID, save_sim=results_dir, 
				report=report, return_results=F, restart=restart, lib_loc=sim_dir)
			
			# Remove parameters
			rm(nruns, parm_list, simID)
		}
	}
}
