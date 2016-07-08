## This script is used to summarize multiple runs from a Core-Transient Simulation
## Accepts command line arguments for where to save results, but by default runs the simulation from the directory where it is called

options(stringsAsFactors=F)

# The following arguments should be passed (in order)
# directory where runs from a single simulation are saved : defaults to current directory
# directory with parameter files : defaults to current directory, parameter files begin with 's_'
# directory with simulation scripts : defaults to current directory
# directory in which to save results : defaults to './Results'
args = commandArgs(trailingOnly=T)

# Set directories
run_dir = ifelse(is.na(args[1]), './', args[1])
parm_dir = ifelse(is.na(args[2]), './', args[2])
sim_dir = ifelse(is.na(args[3]), './', args[3])
results_dir = ifelse(is.na(args[4]), './Summaries/', args[4])

# Load simulation scripts
source(file.path(sim_dir, 'simulation_functions.R'))

# Read in parameter files
# Parameter files are designated by starting with 's_'
file_list = list.files(parm_dir, '^s_')

# Run summaries for each parameter file
for(f in file_list){
	print(paste('Started', f))
	
	# Read in parameters
	parm_file = file.path(parm_dir, f)
	source(parm_file)

	# Summaries for each individual run
	sim_sum_ind = summarize_sim_N(run_dir, breaks=breaks, locs=locs, t_window=t_window, P_obs=P_obs, sum_parms=sum_parms)
	
	# Summaries across runs
	sim_sum = summarize_sim_N(run_dir, breaks=breaks, locs=locs, t_window=t_window, P_obs=P_obs, sum_parms=sum_parms, sum_func=sum_func)
	
	# Save
	save(sim_sum_ind, sim_sum, file=file.path(results_dir, paste0(sumID,'_summary.RData')))

}

q('no')