## This script is used for building and maintaining the CTSim package

library(devtools)
library(roxygen2)

setwd('C:/Users/jrcoyle/Documents/Research/CT-Sim/GitHub/Code/')

current_code = as.package('CTSim')

# Load functions
load_all(current_code)

# Update documentation
document(current_code)





# Testing Functions
myland = make_landscape(c(5,5))
mysp = make_species(10, 10, dist_b=list(type='lognormal', maxN=5, P_maxN=0.01), 
	m=c(.5, .5), r=c(1, .5), dist_d=list(mu=1, var=0), dist_v=list(mu=c(0,.8), var=c(0,0))
)
mygsad = sapply(1:20, function(i) mysp[i, get_sptype(mysp)[i],'b'])
mymeta = populate_landscape(myland, mysp, mygsad, K=10, distribution='uniform', p=.8)

mymeta_t1 = run_timestep(mymeta, myland, mysp, mygsad, d_kernel=list(type='gaussian'),
	v_kernel=list(type='adjacent', moves=2), imm_rate=.2
)

