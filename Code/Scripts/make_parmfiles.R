## This script makes sets of parameter values for CAMM to be run on the cluster

options(stringsAsFactors=F)

# Load CTSim
library(CTSim)

library(pracma) #erf, erfinv

# Set directories
parm_dir = 'Code/Parameters/'
setwd(parm_dir)

# Read in baseline parameter file
source('baseline_parameter_file.txt')


# Run code for desired Experiment below (not all)



#-------------------------------------------------------------------------------------------
### Experiment 1 ###

# ID for this set of parameters
expID = 'EXP1'
expID = 'EXP1-turn' # Used for set of runs that records species turnover 9/5/2016

# For EXP1-turn, save colonization rates
calc_rates=T

# Make directory
dir.create(file.path(parm_dir, expID))

# Define number of runs, for EXP1-turn
nruns=8

## Define dispersal parameters

# Find mean of halfnormal distribution that corresponds to P(x <= D) = 0.99 for D in {1,2,4,8,16,32}
find_d = function(x, p) x/(sqrt(pi)*erfinv(p)) # function that calculates mean of halfnormal under given quantile
D = find_d(2^(0:5), 0.99)

# Adjacent cell dispersal probabilities
adj_P = c(0,.5, 1)

# Set of dispersal parameters
gaus = data.frame(kern='gaussian', d=D, id=paste0('g',2^(0:5)))
adj = data.frame(kern='adjacent', d=adj_P, id=paste0('a', adj_P)); rownames(adj) = adj$id
unif = data.frame(kern='uniform', d=32*sqrt(2), id=paste0('u',round(32*sqrt(2),0)))
d_parms = rbind(gaus, adj,unif)
rownames(d_parms) = d_parms$id

## Define spatial parameters
dcorr=2^(0:4)

# Make parameter files
for(id in d_parms$id){

	# Set dispersal parameters
	d = d_parms[id, 'd']
	dist_d = list(mu=d, var=0)
	d_kernel = list(type=d_parms[id, 'kern'])
	if(d_parms[id, 'kern']=='adjacent') d_kernel = c(d_kernel, moves=1)

	# Make new parameter directory for each set of dispersal parameters
	this_dir = paste('d',id, sep='-')
	dir.create(file.path(expID,this_dir))
	
	# Go through each set of spatial parameters
	for(s in dcorr){
		
		# Set spatial paramters
		vgm_dcorr = s

		# Set movement parameters
		if(id %in% c('g1','g4','g16','u45')){

			# Restricted adult movement
			for(v_id in adj$id){

				# Set movement parameters
				dist_v = list(mu=c(0, d_parms[v_id,'d']), var=c(0,0))
				v_kernel = list(type=d_parms[v_id, 'kern'])
				if(d_parms[v_id, 'kern'] =='adjacent') v_kernel = c(v_kernel, moves=1)
				
				# Set simID
				simID = paste0('d-', id, '_v-', v_id, '_dcorr-', s)

				# Write parameter file
				parmlist = make_parmlist()
				CTSim:::write_parms(parmlist, file.path(expID,this_dir, paste0('p_', simID, '.txt')))
			}
		}	
			
		# Set adult movement = propagule dispersal
		v_id = id
		dist_v = list(mu=c(0, d_parms[v_id, 'd']), var=c(0,0))
		v_kernel = list(type = d_parms[v_id, 'kern'])
		if(d_parms[v_id, 'kern']=='adjacent') v_kernel = c(v_kernel, moves=1)

		# Set simID
		simID = paste0('d-', id, '_v-', v_id, '_dcorr-', s)

		# Write parameter file
		parmlist = make_parmlist()
		CTSim:::write_parms(parmlist, file.path(expID, this_dir, paste0('p_', simID, '.txt')))			
	}
}

# Make directories
for(id in c('g1','g4','g16','u45')){
for(v_id in adj$id){
	dir.create(file.path(expID, paste('d',id, sep='-'), paste('v',v_id, sep='-')))
}}


#-------------------------------------------------------------------------------

# Experiment 2: varying % of the landscape that is habitat A

# ID for this set of parameters
expID = 'EXP2'

# save colonization rates
calc_rates = TRUE

# Define number of runs, for EXP1-turn
nruns=100

hp_parms = seq(0.5, 1, by = 0.1)

# Need to create simID for every parameter combination and then 
# make_parmlist() and write_parms()

# Will need nested loops as multiple variables are explored simultaneously
# Multiple folders needed if you wanted to submit multiple runs on cluster
# (see EXP 1 example)

# Make parameter files
for(id in hp_parms){
  # Set proportion of landscape as habitat A
  habA_prop = id
  
  # Set simID
  simID = paste0('hp-', id)
  
  # Write parameter file
  parmlist = make_parmlist()
  CTSim:::write_parms(parmlist, file.path(expID, paste0('p_', simID, '.txt')))

}



#---------------------------------------------------------------------------------

# ID for this set of parameters
expID = 'EXP3'

# save colonization rates
calc_rates = TRUE

# Define number of runs, for EXP1-turn
nruns=50

# Make directory
dir.create(file.path(parm_dir, expID))

## Define dispersal parameters

# Find mean of halfnormal distribution that corresponds to P(x <= D) = 0.99 for D in {1,2,4,8,16}
find_d = function(x, p) x/(sqrt(pi)*erfinv(p)) # function that calculates mean of halfnormal under given quantile
D = find_d(2^(0:4), 0.99)

# Set of dispersal parameters
gaus = data.frame(kern='gaussian', d=D, id=paste0('g',2^(0:4)))
d_parms = gaus
rownames(d_parms) = d_parms$id

# Set spatial correlation structure parameter
vgm_dcorr = 8

# Two immigration rate scenarios (both less than 0.01 of EXP 1 and 2)
imm_rates = c(0, 0.001)

# Make parameter files
for(id in d_parms$id){

  # Make new parameter directory for each set of dispersal parameters
  this_dir = paste('d',id, sep='-')
  dir.create(file.path(expID,this_dir))
  
  for (imm in imm_rates) {
    # Set dispersal parameters
    d = d_parms[id, 'd']
    dist_d = list(mu=d, var=0)
    d_kernel = list(type=d_parms[id, 'kern'])
    
    dist_v = list(mu=c(0, d_parms[id,'d']), var=c(0,0))
    v_kernel = list(type=d_parms[id, 'kern'])
    
    # Set immigration
    imm_rate = imm
    
    # Set simID
    simID = paste0('d-', id, '_imm-', imm)

    # Write parameter file
    parmlist = make_parmlist()
    CTSim:::write_parms(parmlist, file.path(expID, this_dir, paste0('p_', simID, '.txt')))

  }  
}

######### EXP 4 TEST ############
# goal is vary dispersal and landscape similarity

# ID for this set of parameters
expID = 'EXP4'

# save colonization rates
calc_rates = TRUE

# Define number of runs, for EXP1-turn
nruns=50

# Make directory
dir.create(file.path(parm_dir, expID))

## Define dispersal parameters

# Find mean of halfnormal distribution that corresponds to P(x <= D) = 0.99 for D in {1,2,4,8,16}
find_d = function(x, p) x/(sqrt(pi)*erfinv(p)) # function that calculates mean of halfnormal under given quantile
D = find_d(2^(0:4), 0.99)

# Set of dispersal parameters
gaus = data.frame(kern='gaussian', d=D, id=paste0('g',2^(0:4)))
d_parms = gaus
rownames(d_parms) = d_parms$id

# Set spatial correlation structure parameter
vgm_dcorr = 8

# Two immigration rate scenarios (both less than 0.01 of EXP 1 and 2)
imm_rates = 0.01


hp_parms = seq(0.5, 1, by = 0.1)

# Make parameter files
for(did in d_parms$id){
  
  # Make new parameter directory for each set of dispersal parameters
  this_dir = paste('d',did, sep='-')
  dir.create(file.path(expID,this_dir))
  
    # Need to create simID for every parameter combination and then 
    # make_parmlist() and write_parms()
    
    # Will need nested loops as multiple variables are explored simultaneously
    # Multiple folders needed if you wanted to submit multiple runs on cluster
    # (see EXP 1 example)
    
    # Make parameter files
    for(id in hp_parms){
      # Set proportion of landscape as habitat A
      habA_prop = id
      
      # Set simID
      simID = paste0('hp-', id)
      
      # Set dispersal parameters
    d = d_parms[did, 'd']
    dist_d = list(mu=d, var=0)
    d_kernel = list(type=d_parms[id, 'kern'])
    
    dist_v = list(mu=c(0, d_parms[id,'d']), var=c(0,0))
    v_kernel = list(type=d_parms[id, 'kern'])
    
    # Set simID
    simID = paste0('d-', did, '_hp-', id)
    print(simID)
       
    # Write parameter file
    parmlist = make_parmlist()
    CTSim:::write_parms(parmlist, file.path(expID, this_dir, paste0('d_', simID, '.txt')))
    }
}  





