## This script is used to run a Core-Transient Simulation
## Accepts command line arguments for where to save results, but by default runs the simulation from the directory where it is called

options(stringsAsFactors=F)


# The following arguments should be passed (in order)
# number of cores allocated to this run : defaults to 1
# directory with parameter files : defaults to current directory
# directory with simulation scripts : defaults to current directory
# directory in which to save results : defaults to './Results'
args = commandArgs(trailingOnly=T)

# Set directories
parm_dir = ifelse(is.na(args[2]), './', args[2])
sim_dir = ifelse(is.na(args[3]), './', args[3])
results_dir = ifelse(is.na(args[4]), './Results/', args[4])

# Load simulation scripts
source(paste0(sim_dir, 'simulation_functions.R'))

# Read in parameter files
# Parameter files are designated by starting with 'p_'
file_list = list.files(parm_dir, '^p_')

# Set number of cores
ncores = ifelse(is.na(args[1]), 1, as.numeric(args[1]))

# Run set of simulations on each parameter file
for(f in file_list){
	
	# Read in parameters
	parm_file = paste0(parm_dir, f)
	source(parm_file)
	parm_list = make_parmlist()

	# Run CAMM
	sim_results = run_sim_N(nruns, parm_list, ncores, simID, sim_dir=sim_dir, save_sim=results_dir, report=100, return_results=F)

	# Summarize simulation

}

