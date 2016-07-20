#' Restart simulation 
#' 
#' Script used to restart simulations on a set of parameter files
#' from the command line using \code{R CMD BATCH}.
#'
#' @name restart_simulation
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
#' @usage 
#' \code{R CMD BATCH '--args ncores parm_dir results_dir sim_dir report'  
#' 	restart_simulation.r outfile.rout}
#'
#' @seealso \code{\link{run_sim_P}} for a description of parameters 
#' that can be passed through args
#' 

# Set options and get arguments
options(stringsAsFactors=F)
args = commandArgs(trailingOnly=T)

# Set number of cores
ncores = ifelse(is.na(args[1]), 1, as.numeric(args[1]))

# Set directories
parm_dir = ifelse(is.na(args[2]), './', args[2])
results_dir = ifelse(is.na(args[3]), './Results/', args[3])
sim_dir = ifelse(is.na(args[4]), './', args[4])

# Set reporting interval
report = ifelse(is.na(args[5]), 0, args[5])

# Load CTSim package
tryCatch(library(CTSim), error=function(e) library(CTSim, lib.loc=sim_dir))

# Run function
run_sim_P(ncores, parm_dir, results_dir, sim_dir, report, restart=T)

# Quit R without saving workspace
quit('no')

