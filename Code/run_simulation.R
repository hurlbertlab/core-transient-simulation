## This script is used to run a Core-Transient Simulation


sim_dir = 'C:/Users/jrcoyle/Documents/Research/CT-Sim/GitHub/Code/'
source(paste0(sim_dir, 'simulation_functions.R'))



# TESTING

nruns = 3
nparallel = 2
simID='test'
sim_dir = './Research/CT-Sim/GitHub/Code/'
save_sim = './Research/CT-Sim/testing'

parms = list(
	dimX = 10,
	dimY = 10,
	vgm_dcorr = 1,
	habA_prop = 0.5,
	S_A = 10,
	S_B = 10,
	S_AB = 5,
	dist_b = list(type='lognormal', maxN=5, P_maxN=0.001),
	m_rates = c(.1, .2, .1),
	r_rates = c(.9, .5, .9),
	dist_d = list(mu=1.5, var=0.5),
	dist_gsad = list(type='same'),
	K = 20,
	prop_full = .5,
	init_distribute = 'designated',
	cells_distribute = data.frame(x=5, y=5, n=15),
	nsteps = 20,
	d_kernel = list(type='gaussian'),
	imm_rate = 0.05,
	save_steps = seq(0, 20, 2)
)

sim_results = run_sim_N(nruns, parms, nparallel, simID, sim_dir=sim_dir)


## Examine species richness trajectories through time
use_locs = as.matrix(expand.grid(x=1:10, y=1:10))
abun_profs = calc_abun_profile(use_locs, t_window=list(start=1, stop=20), sim_results$results[[1]], 25)
occupancy = calc_occupancy(abuns=abun_profs, do_freq=F)
tot_rich = calc_rich(abuns=abun_profs)
ct_rich = calc_rich_CT(abun_profs, occupancy, breaks=c(0.25,.5,.75), agg_times=list(1:5, 2:6, 3:7, 4:8, 5:9, 6:10, 7:11))



# TESTING
myland = make_landscape(10,10)
mysp = make_species(S_A=10, S_B=15, S_AB=3, dist_b = list(maxN=5, P_maxN=0.01, type='poisson'), m=c(0.1, 0.5), r =c(.9,.7), dist_d=list(mu=3, var=0.5))
mycomm = populate_landscape(myland, mysp, K=20, p=.5, distribution='designated', which_cells = data.frame(x=1,y=1,n=15))
mygsad = make_sad(dim(mysp)[1], distribution=list(type='same'))

mycomm = populate_landscape(myland, mysp, K=100, p=.5, distribution='designated', which_cells=data.frame(x=1,y=1,n=80))

myresults = run_sim(20, mycomm, myland, mysp, mygsad, imm_rate=.02, save_steps=seq(0,20,2))
myresults = run_timestep(mycomm, myland, mysp, mygsad, imm_rate=0.02)

which_cells = data.frame(x=1,y=1,n=80)

# TESTING

parm_list = list(
	d_kernel=list(type='gaussian'),
	imm_rate = 0.2
)
test_run = run_sim(10, mycomm, myland, mysp, mygsad, parm_list, save_steps=seq(2,10,2))
test_run = run_sim(2, mycomm, myland, mysp, mygsad)


# TESTING
mylocs = expand.grid(x=seq(2,30,3), y=seq(3,30,3))
abun_profiles = calc_abun_profile(mylocs, list(start=0, stop=10), test_run, length(mygsad))
calc_occupancy(mylocs, list(start=0, stop=10), test_run, length(mygsad), do_freq=T)
calc_occupancy(abuns = abun_profiles, which_species=10, agg_times=list(1:5))
calc_rich(mylocs, list(start=0, stop=10), test_run, length(mygsad), which_species=1:10, agg_times=2)
calc_rich(mylocs, list(start=0, stop=10), test_run, length(mygsad))


for(i in 1:nrow(locs)){
	plot_abun(abun_profiles[,,i])
	Sys.sleep(2)
}

abun_profiles

abun_profiles = calc_abun_profile(locs = as.matrix(data.frame(x=5, y=5)), t_window=list(start=1, stop=20), sim = foo, N_S=dim(this_species)[1])



*dimX
*dimY
vgm_mod
vgm_dcorr
habA_prop
*S_A
*S_B
S_AB
dist_b
*m_rates
*r_rates
dist_d
dist_gsad
K
prop_full
init_distribute
cells_distribute
*nsteps
d_kernel
imm_rate


