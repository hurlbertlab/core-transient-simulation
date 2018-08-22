## This script is used for building and maintaining the CTSim package
#devtools::install_github("hadley/staticdocs")

# Load packages needed to build CTSim
library(devtools)
library(roxygen2)
library(gstat)

setwd('C:/git/core-transient-simulation/Code/')
# /Users/sheldontaylor/core-transient-simulation/Code/

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
install.packages('CTSim_0.1.7.zip')
library(CTSim)

# Make static html documentation
#build_site('CTSim')

################################################### DO NOT RUN THIS, D/N WORK
### Misc code for testing functions
library(raster)

# Testing Functions
myland = make_landscape(c(5,5), prop=.5)
mysp = make_species(3, 4, dist_b=list(type='same', p1=2), 
	m=c(.5, .5), r=c(1, .5), dist_d=list(mu=1, var=0), dist_v=list(mu=c(0,.8), var=c(0,.02), type='adjacent')
)
mygsad = sapply(1:7, function(i) mysp[i, get_sptype(mysp)[i],'b'])
mymeta = populate_landscape(myland, mysp, mygsad, K=10, distribution='uniform', p=.8)

mymeta_t1 = run_timestep(mymeta, myland, mysp, mygsad, d_kernel=list(type='gaussian'),
	v_kernel=list(type='adjacent', moves=2), imm_rate=.2
)

mymeta_t20 = run_sim(20, mymeta, myland, mysp, mygsad, d_kernel=list(type='gaussian'),
	v_kernel=list(type='adjacent', moves=2), imm_rate=.2, save_steps=seq(2, 20, 2),
	report=2, ID='testrun', calc_rates=T
)
mymeta_t8 = run_sim(8, mymeta, myland, mysp, mygsad, d_kernel=list(type='gaussian'),
	v_kernel=list(type='adjacent', moves=2), imm_rate=.2, save_steps=seq(2, 8, 2),
	report=2, ID='testrun'
)

# Testing summary functions
setwd('/Users/sheldontaylor/core-transient-simulation/scratch')
run_sim_P(2, report=2)

run_sum = summarize_sim_P('d-a0_v-a0_dcorr-1')
load('Summaries/test_summary.RData') # d/n work


# Check summary manually
load('Results/d-a0_v-a0_dcorr-1/d-a0_v-a0_dcorr-1_run1.RData')



mylocs = aggregate_cells(X=c(5,5), dX=c(2,2), form='partition')
mylocs = aggregate_cells(X=c(5,5), dX=c(1,1), form='partition')

calc_abun(mymeta_t1, N_S=20, only_species=T)

myabuns = calc_abun_profile(mylocs, list(start=1, stop=10), mymeta_t10, 20) 

myocc = calc_occupancy(abuns = myabuns)

calc_abun_CT(myabuns, myocc, seq(.1, .9, .1))

hab = sapply(mylocs, function(x) average_habitat(x, myland))

trates = calc_species_turnover(mymeta_t20$turnover, mylocs, 20, which_species=c(1,1,1,2,2,2,2))

tapply(trates[1,,'gain',1], hab, mean)
tapply(trates[1,,'gain',2], hab, mean)
tapply(trates[1,,'loss',1], hab, mean)
tapply(trates[1,,'loss',2], hab, mean)

data.frame(hab, trates[1,,,1], trates[1,,,2])


myobs = sapply(c(.5, .8), function(p) sample_sim(myabuns, p), simplify='array')

mysumA= summarize_sim(mymeta_t20, .5, mylocs, list(start=10, stop=20), mysp, myland, mygsad, 
	P_obs=1, sum_parms=list(time_sum='none', hab='A'),sum_turn=T)

mysumA= summarize_sim(mymeta_t20, .5, mylocs, list(start=10, stop=20), mysp, myland, mygsad, 
	P_obs=1, sum_parms=list(time_sum='none', hab='A'))

mysumB= summarize_sim(mymeta_t20, .5, mylocs, list(start=10, stop=20), mysp, myland, mygsad, 
	agg_times=3, P_obs=1, sum_parms=list(time_sum='mean', hab='B'))

mysum =  summarize_sim(mymeta_t8, .5, mylocs, list(start=1, stop=8), mysp, myland, mygsad, 
	P_obs=1, sum_parms=list(time_sum='none', hab='AB'))

hab = aggregate_hab_type(myland, mylocs)
mylocs


mymeta_t10$turnover

