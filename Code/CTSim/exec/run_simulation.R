#' Run simulation 
#' 
#' Script used to run simulations on a set of parameter files
#' from the command line using \code{R CMD BATCH}.
#'
#' @name run_simulation
#' @inheritParams run_sim_P
#' @usage 
#' \code{R CMD BATCH '--args ncores parm_dir results_dir sim_dir report'  
#' 	run_simulation.r outfile.rout}
#'
#' @seealso \code{\link{run_sim_P}} for a description of parameters 
#' that can be passed through args

# Set options and get arguments
options(stringsAsFactors=F)
args = commandArgs(trailingOnly=T)

# Set number of cores
ncores = ifelse(is.na(args[1]), 1, as.numeric(args[1]))

# Set directories
parm_dir = ifelse(is.na(args[2]), './', args[2])
results_dir = ifelse(is.na(args[3]), './Results/', args[3])
sim_dir = ifelse(is.na(args[4]), NULL, args[4])

# Set reporting interval
report = ifelse(is.na(args[5]), 0, args[5])

# Load CTSim package
tryCatch(library(CTSim), error=function(e) library(CTSim, lib.loc=sim_dir))

# Run function
run_sim_P(ncores, parm_dir, results_dir, sim_dir, report)

# Quit R without saving workspace
quit('no')

