#' Summarize simulations 
#' 
#' Script used to summarize simulations on a set of parameter files
#' from the command line using \code{R CMD BATCH}.
#'
#' @name summarize_simulation
#' @inheritParams summarize_sim_P
#' @param sim_dir directory where \code{CTSim} is installed, 
#' 	if not on the default search path
#' @usage 
#' \code{R CMD BATCH '--args ncores run_dir parm_dir results_dir sim_dir'  
#' 	summarize_simulation_cross_time.R outfile.Rout}
#'
#' @seealso \code{\link{summarize_sim_P}} for a description of parameters 
#' that can be passed through args

# Set options and get arguments
options(stringsAsFactors=F)
args = commandArgs(trailingOnly=T)

# Set directories
run_dir = ifelse(is.na(args[1]), './', args[1])
parm_dir = ifelse(is.na(args[2]), './', args[2])
results_dir = ifelse(is.na(args[3]), './Summaries/', args[3])
sim_dir = ifelse(is.na(args[4]), './', args[4])

# Load package
tryCatch(library(CTSim), error=function(e) library(CTSim, lib.loc=sim_dir))

# Print versions of packages being used
sessionInfo()

# Run summary function
summarize_sim_P(run_dir, parm_dir, results_dir, cross_time=F)

# Quit R without saving workspace
q('no')