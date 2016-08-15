## This script is used for building and maintaining the CTSim package
#devtools::install_github("hadley/staticdocs")

# Load packages needed to build CTSim
library(devtools)
library(roxygen2)
library(staticdocs)

setwd('C:/Users/jrcoyle/Documents/Research/CT-Sim/GitHub/Code/')

current_code = as.package('CTSim')

# Load functions
load_all(current_code)

# Update documentation
document(current_code)

# Add Imports and Suggests to DESCRIPTION
setwd('./CTSim')
use_package('abind')
use_package('fdrtool')
use_package('gstat')
use_package('poweRlaw')
use_package('raster')
use_package('reshape2')
use_package('sads')
use_package('sp')
use_package('doParallel','Suggests')
use_package('foreach','Suggests')

# Check the package
setwd('../')
check('CTSim')

# Build the package
build('CTSim')
#build_win('CTSim')


# Check install
install.packages('CTSim_0.1.3.zip', repos=NULL)
library(CTSim)

# Make static html documentation
build_site('CTSim', 'HTML')

###################################################
### Misc code for testing functions


# Testing Functions
myland = make_landscape(c(5,5))
mysp = make_species(3, 4, dist_b=list(type='lognormal', maxN=5, P_maxN=0.01), 
	m=c(.5, .5), r=c(1, .5), dist_d=list(mu=1, var=0), dist_v=list(mu=c(0,.8), var=c(0,.02), type='adjacent')
)
mygsad = sapply(1:7, function(i) mysp[i, get_sptype(mysp)[i],'b'])
mymeta = populate_landscape(myland, mysp, mygsad, K=10, distribution='uniform', p=.8)

mymeta_t1 = run_timestep(mymeta, myland, mysp, mygsad, d_kernel=list(type='gaussian'),
	v_kernel=list(type='adjacent', moves=2), imm_rate=.2
)

mymeta_t10 = run_sim(10, mymeta, myland, mysp, mygsad, d_kernel=list(type='gaussian'),
	v_kernel=list(type='adjacent', moves=2), imm_rate=.2, save_steps=seq(2, 10, 2),
	report=2, ID='testrun'
)

run_sim_P(2, report=2)

mylocs = aggregate_cells(X=c(5,5), dX=c(2,2), form='partition')

calc_abun(mymeta_t1, N_S=20, only_species=T)

myabuns = calc_abun_profile(mylocs, list(start=6, stop=10), mymeta_t10, 20) 

myocc = calc_occupancy(abuns = myabuns)

calc_abun_CT(myabuns, myocc, seq(.1, .9, .1))

sapply(mylocs, function(x) average_habitat(x, myland))


myobs = sapply(c(.5, .8), function(p) sample_sim(myabuns, p), simplify='array')

mysum= summarize_sim(mymeta_t10, .5, mylocs, list(start=6, stop=10), mysp, myland, mygsad, 
	P_obs=c(.5, .8), sum_parms=list(time_sum='none'))



	abuns_global = sapply(1:2, function(i){
		apply(myobs[,,,i], 1, function(x) calc_abun(x, 20, only_species=T))
	}, simplify='array')



