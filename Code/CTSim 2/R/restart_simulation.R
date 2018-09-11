#' Restart simulation 
#' 
#' Script used to restart simulations on a set of parameter files
#' from the command line using \code{R CMD BATCH}.
#'
#' This script is found in the exec/ directory where the package is 
#' installed. To run, move script to working directory and use:
#'
#' \code{R CMD BATCH '--args ncores parm_dir results_dir sim_dir report'  
#' 	restart_simulation.R outfile.Rout}
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
#'
#' @seealso \code{\link{run_sim_P}} for a description of parameters 
#' that can be passed through args
#' 
#' @name restart_simulation
NULL
