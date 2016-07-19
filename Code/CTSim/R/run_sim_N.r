#' Run multiple simulations
#'
#' Runs multiple independent simulations on the same set of parameters.
#'
#' This function runs multiple independent simulations on a set of parameters 
#' given in \code{parms}. Each run generates a new landscape, species pool, 
#' global species abundance distribution and initial metacommunity.
#' Landscapes, species pools and gsads are saved as lists 
#' (object names: \code{lands_N}, \code{species_N}, \code{gsad_N})
#' in \code{<simID>_simobjects.RData}.
#' Simulations can be run in parallel by specifying 
#' \code{nparallel > 1}, which requires the \code{\link[doParallel]{doParallel}} and
#' \code{\link[foreach]{foreach}} packages.
#' By default, \code{nparallel = 1} and the simulations proceed serially.
#' Each run of the simulation is temporarily saved to the working directory
#' or permanently saved to the directory specified by \code{save_sim}. 
#' If this directory does not exist then it is created. Runs are saved as
#' \code{<simID>_run<i>.RData}. This file contains five objects:
#' \describe{
#'	\item{results}{an array of the metacommunity through time returned by 
#'		\code{\link{run_sim}}.}
#'	\item{this_land}{the landscape used for the simulation}
#'	\item{this_species}{the species pool used in the simulation}
#'	\item{this_gsad}{the global species abundance distribution used in
#'		the simulation}
#'	\item{this_metacomm}{the initial metacommunity}
#'}
#' If \code{return_results} is \code{TRUE}, then 
#' after all simulations are complete, all runs are read back into memory,  
#' compiled arrays, and returned as a list by the function. 
#' List components are \code{results}, \code{species}, \code{lands},
#' \code{gsads} and the first dimension of each array is the run to which it
#' corresponds. For example, \code{results[1,,,]} is the metacommunity 
#' from run 1. If \code{return_results} is \code{FALSE} then temporary
#' files from the simulation runs are NOT removed, even if \code{save_sim} is
#' not specified.
#'
#' This function can be used to restart a set of multiple simulation runs, 
#' but does not currently allow users to restart a simulation on an existing 
#' run. If \code{restart} is \code{TRUE}, the function navigates to the 
#' \code{save_sim} directory and searches for the first instance of 
#' code{i <= nruns} where \code{<simID>_run<i>.RData} does not exist.
#' It then starts simulations for all \code{i} that do not have saved files
#' using the objects saved in \code{[simID]_simobjects.RData}.
#'
#' @note Users should be cautious in specifying \code{return_results=TRUE}
#' for large simulations where memory requirements may not support large 
#' arrays.
#'
#' @param nruns (required) number of replicate simulations
#' @param parms (required) list of parameters for running simulation as 
#' 	created by \code{\link{make_parmlist}}
#' @param nparallel number of cores to run the simulations on in parallel. 
#' 	Defaults to 1.
#' @param simID character string identifying this simulation run. 
#' 	Defaults to 'test'.
#' @param save_sim directory in which to save simulation results. 
#' 	If none specified simulation results are not saved.
#' @param report number of timesteps after which to report status. 
#' 	Defaults to 0, no reporting.
#' @param return_results logical indicating whether function should
#' 	return simulation results and objects. See details. Defaults to TRUE.
#' @param restart logical indicating whether the simulation should continue 
#' 	from a set of saved runs. See details.
#' @param lib_loc location where \code{CTSim} is installed, if not on default
#' 	search path
#' @return nothing or a list of simulation results and objects. See details.
#'
#' @seealso \code{\link{run_sim}} for details on how each simulation runs \cr
#' \code{\link{make_parmlist}} for parameters that can be passed to the 
#' simulation
#' @export

run_sim_N = function(nruns, parms, nparallel=1, simID='test', save_sim=NULL, report=0, return_results=T, restart=F, lib_loc=NULL){
	sim_complete=F
	
	# Define working directory to save simulation runs
	save_dir = ifelse(is.null(save_sim), getwd(), save_sim)
	if(!file.exists(save_dir)) dir.create(save_dir)
	save_dir = file.path(save_dir, simID)
	if(!file.exists(save_dir)){
		dir.create(save_dir)
	} else {
		warning('Simulation run already exists and may be overwritten unless restart = TRUE')
	}
	
	# Simulate multiple runs in parallel
	if(nparallel > 1){
		# Attempt to load doParallel
		if(requireNamespace('doParallel', quietly=TRUE)&requireNamespace('foreach', quietly=TRUE)){
			# Attach functions in doParallel and foreach
			library(doParallel)
		
			# Make and register cluster
			cluster = makeCluster(nparallel, outfile=paste0(simID, '.Rout'))
			registerDoParallel(cluster)

			# Send required functions and objects to each node
			clusterExport(cluster, c('parms','simID','save_sim','report','save_dir','lib_loc'), envir=environment())
			clusterEvalQ(cluster, library(CTSim, lib.loc=lib_loc))
			
			# If this is not a restart of a previous run
			if(!restart){

				# Initialize simulation landscapes
				lands_N = parLapply(cluster, 1:nruns, function(j){
					with(parms, {
						x = dimX
						y = dimY
						if(!exists('vgm_mod')) vgm_mod = NULL
						d = ifelse(exists('vgm_dcorr'), vgm_dcorr, NA)
						prop = ifelse(exists('habA_prop'), 1-habA_prop, 0.5)
						make_landscape(x, y, vgm_mod, d, prop, draw_plot=F)
					})
				})
		
				# Report progress
				if(report>0) print(paste0(Sys.time(), ': Finished making landscapes.'))	

				# Initialize species vital rates
				species_N = parLapply(cluster, 1:nruns, function(j){
					with(parms, {
						S_AB = ifelse(exists('S_AB'), S_AB, NA) 
						if(!exists('dist_b')) dist_b = NULL
						m = m_rates
						r = r_rates
						if(!exists('dist_d')){
							if(exists('d_kernel')){
								dist_d = list(type=d_kernel$type)
							} else { dist_d = NULL }
						} else {
							if(exists('d_kernel')) dist_d = c(dist_d, type=d_kernel$type)
						}
						if(!exists('dist_v')){
							if(exists('v_kernel')){
								dist_v = list(type=v_kernel$type)
							} else { dist_v = NULL }
						} else {
							if(exists('v_kernel')) dist_v = c(dist_v, type=v_kernel$type)
						}
						make_species(S_A, S_B, S_AB, dist_b, m, r, dist_d, dist_v)
					})
				})

				# Report progress
				if(report>0) print(paste0(Sys.time(), ': Finished making species pools.'))

				# Send initial landscapes and species to all cluster cores
				clusterExport(cluster, c('lands_N','species_N'), envir=environment())
			
				# Initialize global species abundance distribution
				# dist_gsad can be a generic distribution or 'b_rates' indicating that it should be the same as species birth rates
				gsad_N = parLapply(cluster, 1:nruns, function(j){
					with(parms, {
						N_S = dim(species_N[[j]])[1]
						if(exists('dist_gsad')){
							if(is.list(dist_gsad)){
							
								# Use specified distribution to generate abundances
								distribution = dist_gsad
								gsad_vec = make_sad(N_S, distribution)

							} else {

								# Make global abundaces equal to species birth rates	
								if(dist_gsad=='b_rates'){
				
									A_rates = species_N[[j]][1:S_A,'A','b']
									B_rates = species_N[[j]][(S_A+1):(S_A+S_B),'B','b']
									gsad_vec = c(A_rates, B_rates)
									if(exists('S_AB')) if(S_AB > 0) gsad_vec = c(gsad_vec, rowMeans(species_N[[j]][(S_A+S_B+1):(S_A+S_B+S_AB),,'b']))
						
								} else {
									stop('Unrecognized value for parameter dist_gsad.')
								}
							}

						# Defaults to same abundance for each species
						} else {	
							distribution = list(type='same')
							gsad_vec = make_sad(N_S, distribution)
						}
					
						# Return vector of global abundances
						gsad_vec
					})
				})

				# Report progress
				if(report>0) print(paste0(Sys.time(), ': Finished making gsads.'))

				# Save for later restart
				save(lands_N, species_N, gsad_N, file=file.path(save_dir, 'sim_objects.RData'))

				# Send global species abundance distributions to cluster
				clusterExport(cluster, 'gsad_N', envir=environment())

			# If this is a restart of a previous run
			} else {
				
				# Read in lands, species, gsads from directory where simulation results saved
				load(file.path(save_dir, 'sim_objects.RData'))
				
				# Export objects to cluster
				clusterExport(cluster, c('lands_N','species_N','gsad_N'), envir=environment())
			}

			# Run simulations using foreach to reduce memory requirements
			foreach(j=1:nruns) %dopar% {

				# Define file to save results
				this_runfile = file.path(save_dir, paste0(simID, '_run', j, '.RData'))

				# Check whether this is a restart and whether this run has already be done
				if(restart & file.exists(this_runfile)){
					
					if(report>0) print(paste0(Sys.time(), ': Skipping run ', j))
				
				} else {

					if(report>0) print(paste0(Sys.time(), ': Start run ', j))
					
					# Define this landscape and species pool
					this_land = lands_N[[j]]
					this_species = species_N[[j]]
					this_gsad = gsad_N[[j]]
				
					# Distribute species across landscape
					this_metacomm = with(parms, {
						p = ifelse(exists('prop_full'), prop_full, NA)
						distribution = ifelse(exists('init_distribute'), init_distribute, NA)
						if(exists('cells_distribute')){
							which_cells = cells_distribute
						} else {
							which_cells = NULL
						}
						populate_landscape(this_land, this_species, this_gsad, K, distribution, p, which_cells)
					})
				
					# Run simulation
					results = with(parms, {
						if(!exists('d_kernel')) d_kernel = NULL
						if(!exists('v_kernel')) v_kernel = NULL
						imm_rate = ifelse(exists('imm_rate'), imm_rate, NA)
						if(!exists('save_steps')) save_steps = NULL
						run_sim(nsteps, this_metacomm, this_land, this_species, this_gsad, d_kernel, v_kernel, imm_rate, save_steps, report, ID=j)
					})

					# Save results
					save(results, this_species, this_land, this_metacomm, this_gsad, file=this_runfile)

					gc()
				}
			}
			
			sim_complete=T
			stopCluster(cluster)

		} else {
			stop('doParallel or foreach not found. Cannot run simulation in parallel without these package.')
		# Simulate runs sequentially
		}
	} else {

		# If this is not a restart of a previous simulation
		if(!restart){

			# Initialize simulation landscapes
			lands_N = lapply(1:nruns, function(j){
				with(parms, {
					x = dimX
					y = dimY
					if(!exists('vgm_mod')) vgm_mod = NULL
					d = ifelse(exists('vgm_dcorr'), vgm_dcorr, NA)
					prop = ifelse(exists('habA_prop'), 1-habA_prop, 0.5)
					make_landscape(x, y, vgm_mod, d, prop, draw_plot=F)
				})
			})
			
			# Report progress
			if(report>0) print(paste0(Sys.time(), ': Finished making landscapes.'))

			# Initialize species vital rates
			species_N = lapply(1:nruns, function(j){
				with(parms, {
					S_AB = ifelse(exists('S_AB'), S_AB, NA) 
					if(!exists('dist_b')) dist_b = NULL
					m = m_rates
					r = r_rates 
					if(!exists('dist_d')) dist_d = NULL
					if(!exists('dist_v')) dist_v = NULL
					make_species(S_A, S_B, S_AB, dist_b, m, r, dist_d, dist_v)
				})
			})
			
			# Report progress
			if(report>0) print(paste0(Sys.time(), ': Finished making species pools.'))

			# Initialize global species abundance distribution	
			gsad_N = lapply(1:nruns, function(j){
				with(parms, {
					N_S = dim(species_N[[j]])[1]
					if(exists('dist_gsad')){
						if(is.list(dist_gsad)){
						
							# Use specified distribution to generate abundances
	 						distribution = dist_gsad
							gsad_vec = make_sad(N_S, distribution)

						} else {

							# Make global abundaces equal to species birth rates	
							if(dist_gsad=='b_rates'){
			
								A_rates = species_N[[j]][1:S_A,'A','b']
								B_rates = species_N[[j]][(S_A+1):(S_A+S_B),'B','b']
								gsad_vec = c(A_rates, B_rates)
								if(exists('S_AB')) if(S_AB > 0) gsad_vec = c(gsad_vec, rowMeans(species_N[[j]][(S_A+S_B+1):(S_A+S_B+S_AB),,'b']))
					
							} else {
								stop('Unrecognized value for parameter dist_gsad.')
							}
						}

					# Defaults to same abundance for each species
					} else {	
						distribution = list(type='same')
						gsad_vec = make_sad(N_S, distribution)
					}
				
					# Return vector of global abundances
					gsad_vec
				})
			})

			# Report prgress
			if(report>0) print(paste0(Sys.time(), ': Finished making gsads.'))

			# Save simulation objects
			save(lands_N, species_N, gsad_N, file=file.path(save_dir, 'sim_objects.RData'))

		# If this is a restart of a previous run
		} else {
			
			# Read in lands, species, gsads from directory where simulation results saved
			load(file.path(save_dir, 'sim_objects.RData'))
		}

		# Run simulations
		for(j in 1:nruns){

			# Define file to save results
			this_runfile = file.path(save_dir, paste0(simID, '_run', j, '.RData'))

			# Check whether this is a restart and whether this run has already be done
			if(restart & file.exists(this_runfile)){
				
				if(report>0) print(paste0(Sys.time(), ': Skipping run ', j))
			
			} else {
				
				# Report progress
				if(report>0) print(paste0(Sys.time(), ': Start run ', j))

				# Define this landscape and species pool
				this_land = lands_N[[j]]
				this_species = species_N[[j]]
				this_gsad = gsad_N[[j]]
			
				# Distribute species across landscape
				this_metacomm = with(parms, {
					p = ifelse(exists('prop_full'), prop_full, NA)
					distribution = ifelse(exists('init_distribute'), init_distribute, NA)
					if(exists('cells_distribute')){
						which_cells = cells_distribute
					} else {
						which_cells = NULL
					}
					populate_landscape(this_land, this_species, this_gsad, K, distribution, p, which_cells)
				})
			
				# Run simulation
				results = with(parms, {
					if(!exists('d_kernel')) d_kernel = NULL
					if(!exists('v_kernel')) v_kernel = NULL
					imm_rate = ifelse(exists('imm_rate'), imm_rate, NA)
					if(!exists('save_steps')) save_steps = NULL
					run_sim(nsteps, this_metacomm, this_land, this_species, this_gsad, d_kernel, v_kernel, imm_rate, save_steps, report, ID=j)
				})

				# Save results
				save(results, this_species, this_land, this_metacomm, this_gsad, file=file.path(save_dir, paste0(simID, '_run', j, '.RData')))
 			
				gc()
			}
		}
		sim_complete=T
	}

	# Read individual runs back into a list
	if(return_results & sim_complete){
		sim_results = lapply(1:nruns, function(j){
			this_run = file.path(save_dir, paste0(simID, '_run', j, '.RData'))
			load(this_run)
			if(is.null(save_sim))  file.remove(this_run)
			
			results
		})
	

		# Save results
		sim_results = list(results = sim_results, species = species_N, lands = lands_N, gsads = gsad_N)

		if(!is.null(save_sim)){
			save(sim_results, file=file.path(save_dir, paste0(simID, '_results.RData')))
		}
		
		# Return results
		sim_results
	}
}
