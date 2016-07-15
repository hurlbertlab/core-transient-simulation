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


##########################
## Functions


make_plot = function(xlim, ylim, xlab=NULL, ylab=NULL){	
	plot.new()
	plot.window(xlim=xlim, ylim=ylim)
	axis(1)
	abline(h=par('usr')[3], lwd=3)
	axis(2, las=1)
	abline(v=par('usr')[1], lwd=3)
	if(!is.null(xlab)) mtext(xlab, 1, 2.5)
	if(!is.null(ylab)) mtext(ylab, 2, 3)

}

##########################





#### 32x32 grid, dispersal: p=0.99 of traveling 32 units or p=0.80 of traveling 16 units ####

sumID = 'converge32_d-9'


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



#### Examine variation across time and between runs

sumID = 'converge64_d-9'

# Define times and scales
endTs = seq(25, 1000, 25)
scales = paste0('L',1:5)

# Read in summaries through time
for(i in 1:length(scales)){

	load(file.path(sumID, paste0(scales[i], '-T', endTs[1], '_summary.RData')))
	this_bio = sim_sum_ind$bio
	this_occ = sim_sum_ind$occ
	this_xclass = sim_sum_ind$xclass

	for(j in 2:length(endTs)){
	
		this_endT = endTs[j]	
	
		this_file = file.path(sumID, paste0(scales[i], '-', 'T', this_endT,'_summary.RData'))
		load(this_file)

		this_bio = abind(this_bio, sim_sum_ind$bio, along=ifelse(j==2, 0, 1))
		this_occ = abind(this_occ, sim_sum_ind$occ, along=ifelse(j==2, 0, 1))
		this_xclass = abind(this_xclass, sim_sum_ind$xclass, along=ifelse(j==2, 0, 1))	
	}
	
	if(i==1){
		bio_ind_T = this_bio
		occ_ind_T = this_occ
		xclass_ind_T = this_xclass
	} else {

		bio_ind_T = abind(bio_ind_T, this_bio, along=ifelse(i==2, 0, 1))
		occ_ind_T = abind(occ_ind_T, this_occ, along=ifelse(i==2, 0, 1))
		xclass_ind_T = abind(xclass_ind_T, this_xclass, along=ifelse(i==2, 0, 1))		
	}
}
	
# Assign dimension names
dimnames(bio_ind_T)[[1]] = 2^(0:4)
dimnames(bio_ind_T)[[2]] = endTs
dimnames(occ_ind_T)[[1]] = 2^(0:4)
dimnames(occ_ind_T)[[2]] = endTs
dimnames(xclass_ind_T)[[1]] = 2^(0:4)
dimnames(xclass_ind_T)[[2]] = endTs
names(dimnames(bio_ind_T)) = c('scale', 'endT', 'run', 'comm_stat','cross_space','category','p_obs')
names(dimnames(occ_ind_T)) = c('scale', 'endT', 'run', 'comm_stat','cross_space','category','p_obs')
names(dimnames(xclass_ind_T)) = c('scale', 'endT', 'run','cross_space','ab_ct','p_obs')

# Save objects for later comparison of dispersal method

#bio_dunif = bio_ind_T; occ_dunif = occ_ind_T; xclass_dunif = xclass_ind_T
#bio_d9 = bio_ind_T; occ_d9 = occ_ind_T; xclass_d9 = xclass_ind_T
#bio_d4 = bio_ind_T; occ_d4 = occ_ind_T; xclass_d4 = xclass_ind_T


# Scale on different pages
# Community statistic on different pages
# Detectability in different plots

nruns=8
p_obs = seq(.1, 1, .1)

pdf(file.path(sum_dir, paste0(sumID, '_ind_run_summary_through_time.pdf')), height=10, width=8)

par(mfrow=c(5, 2))
#par(mfrow=c(10, 2)) 
par(mar=c(2, 2, 1, 1))
par(oma=c(2,2,3,0))
par(lend=1)

for(p in rev(p_obs)){

	# Biological Core/Trans Richness
	for(sp in 2^(0:4)){
		maxS = max(bio_ind_T[as.character(sp), , , 'rich', '97.5%', , as.character(p)])	

		for(type in c('core','trans')){
			make_plot(c(0, 1000), c(0, maxS), xlab=ifelse(p==p_obs[1], 'Time', ''), ylab=ifelse(type=='core', 'Num. species', ''))
	
			if(sp==1&type=='core') mtext('Core', 3, 0)
			if(sp==1&type=='trans') mtext('Transient', 3, 0)
			if(type=='core') text(0, 0.95*maxS, paste('Grain =',sp,'x',sp), pos=4)

			for(i in 1:nruns){
				polygon(c(endTs, rev(endTs)), c(bio_ind_T[as.character(sp), , i,'rich', '2.5%',type,as.character(p)], rev(bio_ind_T[as.character(sp), , i, 'rich', '97.5%',type,as.character(p)])),
					border=NA, col='#00000020')
				lines(endTs, bio_ind_T[as.character(sp), , i, 'rich', 'mean', type, as.character(p)], lwd=2)
			}
		}
	}
	
	mtext(paste('Birth rate-based categories: detection prob. =', p), 3, 1, outer=T)


	# Occupancy Core/Trans Richness
	for(sp in 2^(0:4)){
		maxS = max(occ_ind_T[as.character(sp), , , 'rich', '97.5%', , as.character(p)])	

		for(type in c(dim(occ_ind_T)[6], 1)){
			make_plot(c(0, 1000), c(0, maxS), xlab=ifelse(p==p_obs[1], 'Time', ''), ylab=ifelse(type==dim(occ_ind_T)[6], 'Num. species', ''))
	
			if(sp==1&type==dim(occ_ind_T)[6]) mtext('Core', 3, 0)
			if(sp==1&type==1) mtext('Transient', 3, 0)
			if(type==dim(occ_ind_T)[6]) text(0, 0.95*maxS, paste('Grain =',sp,'x',sp), pos=4)

			for(i in 1:nruns){
				polygon(c(endTs, rev(endTs)), c(occ_ind_T[as.character(sp), , i,'rich', '2.5%',type,as.character(p)], rev(occ_ind_T[as.character(sp), , i, 'rich', '97.5%',type,as.character(p)])),
					border=NA, col='#00000020')
				lines(endTs, occ_ind_T[as.character(sp), , i, 'rich', 'mean', type, as.character(p)], lwd=2)
			}
		}
	}
	
	mtext(paste('Temporal occupancy-based categories: detection prob. =', p), 3, 1, outer=T)

	
	# Biological Core/Trans Abundance
	for(sp in 2^(0:4)){
		maxS = max(bio_ind_T[as.character(sp), , , 'abun', '97.5%', , as.character(p)])	

		for(type in c('core','trans')){
			make_plot(c(0, 1000), c(0, maxS), xlab=ifelse(p==p_obs[1], 'Time', ''), ylab=ifelse(type=='core', 'Num. individuals', ''))
	
			if(sp==1&type=='core') mtext('Core', 3, 0)
			if(sp==1&type=='trans') mtext('Transient', 3, 0)
			if(type=='core') text(0, 0.95*maxS, paste('Grain =',sp,'x',sp), pos=4)

			for(i in 1:nruns){
				polygon(c(endTs, rev(endTs)), c(bio_ind_T[as.character(sp), , i,'abun', '2.5%',type,as.character(p)], rev(bio_ind_T[as.character(sp), , i, 'abun', '97.5%',type,as.character(p)])),
					border=NA, col='#00000020')
				lines(endTs, bio_ind_T[as.character(sp), , i, 'abun', 'mean', type, as.character(p)], lwd=2)
			}
		}
	}
	
	mtext(paste('Birth rate-based categories: detection prob. =', p), 3, 1, outer=T)


	# Occupancy Core/Trans Richness
	for(sp in 2^(0:4)){
		maxS = max(occ_ind_T[as.character(sp), , , 'abun', '97.5%', , as.character(p)])	

		for(type in c(dim(occ_ind_T)[6], 1)){
			make_plot(c(0, 1000), c(0, maxS), xlab=ifelse(p==p_obs[1], 'Time', ''), ylab=ifelse(type==dim(occ_ind_T)[6], 'Num. individuals', ''))
	
			if(sp==1&type==dim(occ_ind_T)[6]) mtext('Core', 3, 0)
			if(sp==1&type==1) mtext('Transient', 3, 0)
			if(type==dim(occ_ind_T)[6]) text(0, 0.95*maxS, paste('Grain =',sp,'x',sp), pos=4)

			for(i in 1:nruns){
				polygon(c(endTs, rev(endTs)), c(occ_ind_T[as.character(sp), , i,'abun', '2.5%',type,as.character(p)], rev(occ_ind_T[as.character(sp), , i, 'abun', '97.5%',type,as.character(p)])),
					border=NA, col='#00000020')
				lines(endTs, occ_ind_T[as.character(sp), , i, 'abun', 'mean', type, as.character(p)], lwd=2)
			}
		}
	}
	
	mtext(paste('Temporal occupancy-based categories: detection prob. =', p), 3, 1, outer=T)
	

	# Cross Classification
	for(sp in 2^(0:4)){
		maxS = max(xclass_ind_T[as.character(sp), , , '97.5%', 1:2 , as.character(p)])	

		for(type in c('core','trans')){
			make_plot(c(0, 1000), c(0, maxS), xlab=ifelse(p==p_obs[1], 'Time', ''), ylab=ifelse(type=='core', 'Num. species', ''))
	
			if(sp==1&type=='core') mtext('Classified Core (by occupancy)', 3, 0)
			if(sp==1&type=='trans') mtext('Classified Transient (by occupancy)', 3, 0)
			if(type=='core') text(0, 0.95*maxS, paste('Grain =',sp,'x',sp), pos=4)

			for(i in 1:nruns){
				polygon(c(endTs, rev(endTs)), c(xclass_ind_T[as.character(sp), , i, '2.5%',paste0('core_',type),as.character(p)], rev(xclass_ind_T[as.character(sp), , i, '97.5%', paste0('core_',type),as.character(p)])),
					border=NA, col='#00000020')
				lines(endTs, xclass_ind_T[as.character(sp), , i, 'mean', paste0('core_',type), as.character(p)], lwd=2)
			}
		}
	}
	
	mtext(paste('Birth rate-based Core Species: detection prob. =', p), 3, 1, outer=T)

	for(sp in 2^(0:4)){
		maxS = max(xclass_ind_T[as.character(sp), , , '97.5%', 3:4 , as.character(p)])	

		for(type in c('core','trans')){
			make_plot(c(0, 1000), c(0, maxS), xlab=ifelse(p==p_obs[1], 'Time', ''), ylab=ifelse(type=='core', 'Num. species', ''))
	
			if(sp==1&type=='core') mtext('Classified Core (by occupancy)', 3, 0)
			if(sp==1&type=='trans') mtext('Classified Transient (by occupancy)', 3, 0)
			if(type=='core') text(0, 0.95*maxS, paste('Grain =',sp,'x',sp), pos=4)

			for(i in 1:nruns){
				polygon(c(endTs, rev(endTs)), c(xclass_ind_T[as.character(sp), , i, '2.5%',paste0('trans_',type),as.character(p)], rev(xclass_ind_T[as.character(sp), , i, '97.5%', paste0('trans_',type),as.character(p)])),
					border=NA, col='#00000020')
				lines(endTs, xclass_ind_T[as.character(sp), , i, 'mean', paste0('trans_',type), as.character(p)], lwd=2)
			}
		}
	}
	
	mtext(paste('Birth rate-based Transient Species: detection prob. =', p), 3, 1, outer=T)

}

dev.off()

### Comparing dispersal on 32 x 32 grid

bio = abind(bio_d4, bio_d9, bio_dunif, along=0, use.dnns=T); dimnames(bio)[[1]] = c('d-4','d-9','d-unif')
occ = abind(occ_d4, occ_d9, occ_dunif, along=0, use.dnns=T); dimnames(occ)[[1]] = c('d-4','d-9','d-unif')
xclass = abind(xclass_d4, xclass_d9, xclass_dunif, along=0, use.dnns=T); dimnames(xclass)[[1]] = c('d-4','d-9','d-unif')

dnames = c('Gaussian (d=4)', 'Gaussian (d=9)', 'Uniform (unlimited)'); names(dnames) = dimnames(bio)[[1]]

# Just look at dectibility = 1

# Scale on different pages
pdf(file.path(sum_dir, 'compare_dispersal_ind_run_summary_through_time.pdf'), height=8, width=8)

par(mfrow=c(3, 2))
par(mar=c(2, 2, 1, 1))
par(oma=c(2,2,3,0))
par(lend=1)

p = 1

for(sp in 2^(0:4)){


	# Biological Core/Trans Richness
	maxS = max(bio[,as.character(sp), , , 'rich', '97.5%', , as.character(p)])
	for(d in dimnames(bio)[[1]]){
		for(type in c('core','trans')){
			make_plot(c(0, 1000), c(0, maxS), xlab=ifelse(p==p_obs[1], 'Time', ''), ylab=ifelse(type=='core', 'Num. species', ''))
	
			if(d=='d-4'&type=='core') mtext('Core', 3, 0)
			if(d=='d-4'&type=='trans') mtext('Transient', 3, 0)
			if(type=='core') text(0, 0.98*maxS, dnames[d], pos=4)

			for(i in 1:nruns){
				polygon(c(endTs, rev(endTs)), c(bio[d, as.character(sp), , i,'rich', '2.5%',type,as.character(p)], rev(bio[d, as.character(sp), , i, 'rich', '97.5%',type,as.character(p)])),
					border=NA, col='#00000020')
				lines(endTs, bio[d, as.character(sp), , i, 'rich', 'mean', type, as.character(p)], lwd=2)
			}
		}
	}
	
	mtext(paste('Birth rate-based categories: spatial grain =', sp, 'x', sp), 3, 1, outer=T)

	# Occupancy Core/Trans Richness
	maxS = max(occ[, as.character(sp), , , 'rich', '97.5%', , as.character(p)])	
	for(d in dimnames(bio)[[1]]){
		for(type in c(dim(occ)[7], 1)){
			make_plot(c(0, 1000), c(0, maxS), xlab=ifelse(p==p_obs[1], 'Time', ''), ylab=ifelse(type==dim(occ)[7], 'Num. species', ''))
	
			if(d=='d-4'&type==dim(occ)[7]) mtext('Core', 3, 0)
			if(d=='d-4'&type==1) mtext('Transient', 3, 0)
			if(type==dim(occ)[7]) text(0, 0.98*maxS, dnames[d], pos=4)

			for(i in 1:nruns){
				polygon(c(endTs, rev(endTs)), c(occ[d, as.character(sp), , i,'rich', '2.5%',type,as.character(p)], rev(occ[d, as.character(sp), , i, 'rich', '97.5%',type,as.character(p)])),
					border=NA, col='#00000020')
				lines(endTs, occ[d, as.character(sp), , i, 'rich', 'mean', type, as.character(p)], lwd=2)
			}
		}
	}
	
	mtext(paste('Temporal occupancy-based categories: spatial grain =', sp, 'x', sp), 3, 1, outer=T)

	# Biological Core/Trans Abundance
	maxS = max(bio[d, as.character(sp), , , 'abun', '97.5%', , as.character(p)])	
	for(d in dimnames(bio)[[1]]){
		for(type in c('core','trans')){
			make_plot(c(0, 1000), c(0, maxS), xlab=ifelse(p==p_obs[1], 'Time', ''), ylab=ifelse(type=='core', 'Num. individuals', ''))
	
			if(d=='d-4'&type=='core') mtext('Core', 3, 0)
			if(d=='d-4'&type=='trans') mtext('Transient', 3, 0)
			if(type=='core') text(0, 0.98*maxS, dnames[d], pos=4)

			for(i in 1:nruns){
				polygon(c(endTs, rev(endTs)), c(bio[d, as.character(sp), , i,'abun', '2.5%',type,as.character(p)], rev(bio[d, as.character(sp), , i, 'abun', '97.5%',type,as.character(p)])),
					border=NA, col='#00000020')
				lines(endTs, bio[d, as.character(sp), , i, 'abun', 'mean', type, as.character(p)], lwd=2)
			}
		}
	}

	mtext(paste('Birth rate-based categories: spatial grain =', sp, 'x', sp), 3, 1, outer=T)

	# Occupancy Core/Trans Richness
	maxS = max(occ[, as.character(sp), , , 'abun', '97.5%', , as.character(p)])	
	for(d in dimnames(bio)[[1]]){
		for(type in c(dim(occ)[7], 1)){
			make_plot(c(0, 1000), c(0, maxS), xlab=ifelse(p==p_obs[1], 'Time', ''), ylab=ifelse(type==dim(occ)[7], 'Num. individuals', ''))
	
			if(d=='d-4'&type==dim(occ)[7]) mtext('Core', 3, 0)
			if(d=='d-4'&type==1) mtext('Transient', 3, 0)
			if(type==dim(occ)[7]) text(0, 0.98*maxS, dnames[d], pos=4)

			for(i in 1:nruns){
				polygon(c(endTs, rev(endTs)), c(occ[d, as.character(sp), , i,'abun', '2.5%',type,as.character(p)], rev(occ[d, as.character(sp), , i, 'abun', '97.5%',type,as.character(p)])),
					border=NA, col='#00000020')
				lines(endTs, occ[d, as.character(sp), , i, 'abun', 'mean', type, as.character(p)], lwd=2)
			}
		}
	}

	mtext(paste('Temporal occupancy-based categories: spatial grain =', sp, 'x', sp), 3, 1, outer=T)
	
	# Cross Classification
	maxS = max(xclass[, as.character(sp), , , '97.5%', 1:2 , as.character(p)])	
	for(d in dimnames(xclass)[[1]]){
		for(type in c('core','trans')){
			make_plot(c(0, 1000), c(0, maxS), xlab=ifelse(p==p_obs[1], 'Time', ''), ylab=ifelse(type=='core', 'Num. species', ''))
	
			if(d=='d-4'&type=='core') mtext('Classified Core (by occupancy)', 3, 0)
			if(d=='d-4'&type=='trans') mtext('Classified Transient (by occupancy)', 3, 0)
			if(type=='core') text(0, 0.98*maxS, dnames[d], pos=4)

			for(i in 1:nruns){
				polygon(c(endTs, rev(endTs)), c(xclass[d, as.character(sp), , i, '2.5%',paste0('core_',type),as.character(p)], rev(xclass[d, as.character(sp), , i, '97.5%', paste0('core_',type),as.character(p)])),
					border=NA, col='#00000020')
				lines(endTs, xclass[d, as.character(sp), , i, 'mean', paste0('core_',type), as.character(p)], lwd=2)
			}
		}
	}
	
	mtext(paste('Birth rate-based Core Species: spatial grain =', sp, 'x', sp), 3, 1, outer=T)

	maxS = max(xclass[, as.character(sp), , , '97.5%', 3:4, as.character(p)])	
	for(d in dimnames(xclass)[[1]]){
		for(type in c('core','trans')){
			make_plot(c(0, 1000), c(0, maxS), xlab=ifelse(p==p_obs[1], 'Time', ''), ylab=ifelse(type=='core', 'Num. species', ''))
	
			if(d=='d-4'&type=='core') mtext('Classified Core (by occupancy)', 3, 0)
			if(d=='d-4'&type=='trans') mtext('Classified Transient (by occupancy)', 3, 0)
			if(type=='core') text(0, 0.98*maxS, dnames[d], pos=4)

			for(i in 1:nruns){
				polygon(c(endTs, rev(endTs)), c(xclass[d, as.character(sp), , i, '2.5%',paste0('trans_',type),as.character(p)], rev(xclass[d, as.character(sp), , i, '97.5%', paste0('trans_',type),as.character(p)])),
					border=NA, col='#00000020')
				lines(endTs, xclass[d, as.character(sp), , i, 'mean', paste0('trans_',type), as.character(p)], lwd=2)
			}
		}
	}
	
	mtext(paste('Birth rate-based Transient Species: spatial grain =', sp, 'x', sp), 3, 1, outer=T)
}

dev.off()



##### Frequency distribution plots

sumID = 'converge32_d-9'

# Read in summaries across scales
load(file.path(sumID, 'occdist_L1_summary.RData'))
occ = sim_sum$occ
occ_ind = sim_sum_ind$occ

for(sp in 2:5){

	this_file = file.path(sumID, paste0('occdist_L',sp,'_summary.RData'))
	load(this_file)

	occ = abind(occ, sim_sum$occ, along=ifelse(sp==2, 0, 1))
	occ_ind = abind(occ_ind, sim_sum_ind$occ, along=ifelse(sp==2, 0, 1))
}

# Assign dimension names
dimnames(occ)[[1]] = 2^(0:4); names(dimnames(occ)) = c('scale', names(dimnames(sim_sum$occ)))
dimnames(occ_ind)[[1]] = 2^(0:4); names(dimnames(occ_ind)) = c('scale', 'run', 'comm_stat','cross_space','category','p_obs')
dimnames(occ_ind)[[2]] = 1:8


nruns=8
p_obs = seq(.1, 1, .1)
grains = 2^(0:4)
binwidth=0.05
bins = seq(0, 1, binwidth)

pdf(file.path(sum_dir, paste(sumID, 'occupancy_distributions.pdf')), height=11, width=8)

# Page layout
par(mfrow=c(5, 1))
par(mar=c(2,5,1,1))
par(oma=c(2,0,2,0))
par(lend=1)

# Dectability on separate pages
for(p in as.character(rev(p_obs))){

	# Spatial grain in separate plots
	for(sp in as.character(grains)){
		# Calculate maximum number of species in a category
		maxS = max(occ_ind[sp,,'rich','97.5%',,p])
	
		# Set up plot
		make_plot(c(0,1), c(0,maxS), xlab=ifelse(sp=='16', 'Temporal Occupancy',''), ylab='Num. species')

		# Add each run
		for(i in 1:nruns){
		
			# Add 95% intervals
			rect(bins[-1]-binwidth, occ_ind[sp,i,'rich','2.5%',,p], bins[-1], occ_ind[sp,i,'rich','97.5%',,p], col='#00000020', border=NA)

			# Add means
			segments(bins[-1]-binwidth, occ_ind[sp,i,'rich','mean',,p], bins[-1], occ_ind[sp,i,'rich','mean',,p], lwd=2)
		}

		# Add labels
		mtext(paste('Spatial grain =', sp, 'x', sp), 3, -2, adj=0.05)
		if(sp=='1') mtext(paste('Detection prob. =',p), 3, 0, outer=T)
	}


}

dev.off()































