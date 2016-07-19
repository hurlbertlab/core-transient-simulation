#' Summarize simulations 
#' 
#' Script used to summarize simulations on a set of parameter files
#' from the command line using \code{R CMD BATCH}.
#'
#' This script is found in the exec/ directory where the package is 
#' installed. To run, move script to working directory and use:
#'
#' \code{R CMD BATCH '--args run_dir parm_dir results_dir sim_dir'  
#' 	summarize_simulation.R outfile.Rout}
#'
#' @param run_dir directory where simulation runs are located.
#' 	Defaults to the current directory.
#' @param parm_dir directory where parameter files are located.
#' 	Defaults to the current directory.
#' @param results_dir directory where results should be saved.
#' 	Defaults to 'Summaries' in the current directory.
#' @param sim_dir directory where \code{CTSim} is installed, 
#' 	if not on the default search path
#'
#' @seealso \code{\link{summarize_sim_P}} for a description of parameters 
#' that can be passed through args
#'
#' @name summarize_simulation
NULL