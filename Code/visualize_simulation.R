## This script is used to visualize simulation summaries


# Set options and load libraries
options(stringsAsFactors=F)
library(abind)



# Define working directory
sum_dir = 'C:/Users/jrcoyle/Documents/Research/CT-Sim/GitHub/Results'
sim_dir = 'C:/Users/jrcoyle/Documents/Research/CT-Sim/GitHub/Code'

setwd(sum_dir)

# Load simulation functions
source(file.path(sim_dir, 'simulation_functions.R'))

#### 32x32 grid, dispersal: p=0.99 of traveling 32 units or p=0.80 of traveling 16 units ####

sumID = 'converge32-d9'


# Read in summaries across scales
load(file.path(sumID, 'L1_summary.RData'))
bio = sim_sum$bio
occ = sim_sum$occ
xclass = sim_sum$xclass

bio_ind = sim_sum_ind$bio
occ_ind = sim_sum_ind$occ
xclass_ind = sim_sum_ind$xclass

for(sp in 2:5){

	this_file = file.path(sumID, paste0('L',sp,'_summary.RData'))
	load(this_file)

	bio = abind(bio, sim_sum$bio, along=ifelse(sp==2, 0, 1))
	occ = abind(occ, sim_sum$occ, along=ifelse(sp==2, 0, 1))
	xclass = abind(xclass, sim_sum$xclass, along=ifelse(sp==2, 0, 1))	

	bio_ind = abind(bio_ind, sim_sum_ind$bio, along=ifelse(sp==2, 0, 1))
	occ_ind = abind(occ_ind, sim_sum_ind$occ, along=ifelse(sp==2, 0, 1))
	xclass_ind = abind(xclass_ind, sim_sum_ind$xclass, along=ifelse(sp==2, 0, 1))	
}

# Assign dimension names
dimnames(bio)[[1]] = 2^(0:4); names(dimnames(bio)) = c('scale', names(dimnames(sim_sum$bio)))
dimnames(occ)[[1]] = 2^(0:4); names(dimnames(occ)) = c('scale', names(dimnames(sim_sum$occ)))
dimnames(xclass)[[1]] = 2^(0:4); names(dimnames(xclass)) = c('scale', names(dimnames(sim_sum$xclass)))

dimnames(bio_ind)[[1]] = 2^(0:4); names(dimnames(bio_ind)) = c('scale', 'run', 'comm_stat','cross_space','category','p_obs')
dimnames(occ_ind)[[1]] = 2^(0:4); names(dimnames(occ_ind)) = c('scale', 'run', 'comm_stat','cross_space','category','p_obs')
dimnames(xclass_ind)[[1]] = 2^(0:4); names(dimnames(xclass_ind)) = c('scale', 'run','cross_space','ab_ct','p_obs')

dimnames(bio_ind)[[2]] = 1:8
dimnames(occ_ind)[[2]] = 1:8
dimnames(xclass_ind)[[2]] = 1:8


## Plot statistics vs. detectability for each spatial grain

nruns=8
p_obs = seq(.1, 1, .1)

pdf(file.path(sum_dir, paste0(sumID, '_ind_run_summary.pdf')), height=8, width=8)
par(mfrow=c(2,2))
par(lend=1)
par(oma=c(0,0,3,0))

for(sp in dimnames(bio)[[1]]){

maxN = 100*(as.numeric(sp)^2)

## Biological
use_dat = bio_ind
maxS = max(use_dat[sp, ,'rich','97.5%',,])

# Core richness
par(mar=c(3,4,2,0))
plot.new()
plot.window(xlim=c(0, 1.05), ylim=c(0,maxS)) 
axis(1, at=c(0,p_obs))
abline(h=par('usr')[3], lwd=3)
axis(2, las=1)
abline(v=par('usr')[1], lwd=3)
mtext('Num. species', 2, 3)
mtext('Core', 3, 0)

for(i in 1:nruns){
	polygon(c(p_obs, rev(p_obs)), c(use_dat[sp, i,'rich', '2.5%', 'core',], rev(use_dat[sp, i,'rich', '97.5%', 'core',])),
		border=NA, col='#00000020')
	lines(p_obs, use_dat[sp, i,'rich', 'mean', 'core',], lwd=2)
}

# Trans richness
par(mar=c(3,3,2,1))
plot.new()
plot.window(xlim=c(0, 1.05), ylim=c(0,maxS)) 
axis(1, at=c(0,p_obs))
abline(h=par('usr')[3], lwd=3)
axis(2, las=1)
abline(v=par('usr')[1], lwd=3)
mtext('Transient', 3, 0)

for(i in 1:nruns){
	polygon(c(p_obs, rev(p_obs)), c(use_dat[sp, i,'rich', '2.5%', 'trans',], rev(use_dat[sp, i,'rich', '97.5%', 'trans',])),
		border=NA, col='#00000020')
	lines(p_obs, use_dat[sp, i,'rich', 'mean', 'trans',], lwd=2)
}

# Core abundance
par(mar=c(4,4,1,0))
plot.new()
plot.window(xlim=c(0, 1.05), ylim=c(0,maxN)) 
axis(1, at=c(0,p_obs))
abline(h=par('usr')[3], lwd=3)
axis(2, las=1)
abline(v=par('usr')[1], lwd=3)
mtext('Num. individuals', 2, 3)
mtext('Detection probability', 1, 2.5)

for(i in 1:nruns){
	polygon(c(p_obs, rev(p_obs)), c(use_dat[sp, i,'abun', '2.5%', 'core',], rev(use_dat[sp, i,'abun', '97.5%', 'core',])),
		border=NA, col='#00000020')
	lines(p_obs, use_dat[sp, i,'abun', 'mean', 'core',], lwd=2)
}

# Trans abundance
par(mar=c(4,3,1,1))
plot.new()
plot.window(xlim=c(0, 1.05), ylim=c(0,maxN)) 
axis(1, at=c(0,p_obs))
abline(h=par('usr')[3], lwd=3)
axis(2, las=1)
abline(v=par('usr')[1], lwd=3)
mtext('Detection probability', 1, 2.5)

for(i in 1:nruns){
	polygon(c(p_obs, rev(p_obs)), c(use_dat[sp, i,'abun', '2.5%', 'trans',], rev(use_dat[sp, i,'abun', '97.5%', 'trans',])),
		border=NA, col='#00000020')
	lines(p_obs, use_dat[sp, i,'abun', 'mean', 'trans',], lwd=2)
}

mtext(paste('Birth rate-based categories: spatial grain =', sp, 'x', sp), 3, 1, outer=T)

## Occupancy
use_dat = occ_ind
maxS = max(use_dat[sp, ,'rich','97.5%',,])

# Core richness
par(mar=c(3,4,2,0))
plot.new()
plot.window(xlim=c(0, 1.05), ylim=c(0,maxS)) 
axis(1, at=c(0,p_obs))
abline(h=par('usr')[3], lwd=3)
axis(2, las=1)
abline(v=par('usr')[1], lwd=3)
mtext('Num. species', 2, 3)
mtext('Core', 3, 0)

for(i in 1:nruns){
	polygon(c(p_obs, rev(p_obs)), c(use_dat[sp, i,'rich', '2.5%', dim(use_dat)[5],], rev(use_dat[sp, i,'rich', '97.5%', dim(use_dat)[5],])),
		border=NA, col='#00000020')
	lines(p_obs, use_dat[sp, i,'rich', 'mean', dim(use_dat)[5],], lwd=2)
}

# Trans richness
par(mar=c(3,3,2,1))
plot.new()
plot.window(xlim=c(0, 1.05), ylim=c(0,maxS)) 
axis(1, at=c(0,p_obs))
abline(h=par('usr')[3], lwd=3)
axis(2, las=1)
abline(v=par('usr')[1], lwd=3)
mtext('Transient', 3, 0)

for(i in 1:nruns){
	polygon(c(p_obs, rev(p_obs)), c(use_dat[sp, i,'rich', '2.5%', 1,], rev(use_dat[sp, i,'rich', '97.5%', 1,])),
		border=NA, col='#00000020')
	lines(p_obs, use_dat[sp, i,'rich', 'mean', 1,], lwd=2)
}

# Core abundance
par(mar=c(4,4,1,0))
plot.new()
plot.window(xlim=c(0, 1.05), ylim=c(0,maxN)) 
axis(1, at=c(0,p_obs))
abline(h=par('usr')[3], lwd=3)
axis(2, las=1)
abline(v=par('usr')[1], lwd=3)
mtext('Num. individuals', 2, 3)
mtext('Detection probability', 1, 2.5)

for(i in 1:nruns){
	polygon(c(p_obs, rev(p_obs)), c(use_dat[sp, i,'abun', '2.5%', dim(use_dat)[5],], rev(use_dat[sp, i,'abun', '97.5%', dim(use_dat)[5],])),
		border=NA, col='#00000020')
	lines(p_obs, use_dat[sp, i,'abun', 'mean', dim(use_dat)[5],], lwd=2)
}

# Trans abundance
par(mar=c(4,3,1,1))
plot.new()
plot.window(xlim=c(0, 1.05), ylim=c(0,maxN)) 
axis(1, at=c(0,p_obs))
abline(h=par('usr')[3], lwd=3)
axis(2, las=1)
abline(v=par('usr')[1], lwd=3)
mtext('Detection probability', 1, 2.5)

for(i in 1:nruns){
	polygon(c(p_obs, rev(p_obs)), c(use_dat[sp, i,'abun', '2.5%', 1,], rev(use_dat[sp, i,'abun', '97.5%', 1,])),
		border=NA, col='#00000020')
	lines(p_obs, use_dat[sp, i,'abun', 'mean', 1,], lwd=2)
}

mtext(paste('Temporal occupancy-based categories: spatial grain =', sp, 'x', sp), 3, 1, outer=T)

## Cross-classification
use_dat = xclass_ind
maxS = max(use_dat[sp, ,'97.5%',,])

for(j in 1:4){

	par(mar=c(ifelse(j<=2, 3, 4), ifelse(j %in% c(1,3), 4, 3), ifelse(j<=2, 2, 1), ifelse(j %in% c(1,3), 0, 1)))
	plot.new()
	plot.window(xlim=c(0, 1.05), ylim=c(0,maxS)) 
	axis(1, at=c(0,p_obs))
	abline(h=par('usr')[3], lwd=3)
	axis(2, las=1)
	abline(v=par('usr')[1], lwd=3)
	if(j>=3) mtext('Detection probability', 1, 2.5)
	if(j %in% c(1,3)) mtext('Num. species', 2, 2)
	if(j==1) text(0, maxS, 'Biologically core, classified core', pos=4)
	if(j==2) text(0, maxS, 'Biologically core, classified transient', pos=4)
	if(j==3) text(0, maxS, 'Biologically transient, classified core', pos=4)
	if(j==4) text(0, maxS, 'Biologically transient, classified transient', pos=4)

	for(i in 1:nruns){
		polygon(c(p_obs, rev(p_obs)), c(use_dat[sp, i, '2.5%', j,], rev(use_dat[sp, i, '97.5%', j,])),
			border=NA, col='#00000020')
		lines(p_obs, use_dat[sp, i, 'mean', j,], lwd=2)
	}

}

mtext(paste('Cross-classification: spatial grain =', sp, 'x', sp), 3, 1, outer=T)

}

dev.off()











