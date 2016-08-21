## This script is used to visualize simulation summaries


# Set options and load libraries
options(stringsAsFactors=F)
library(CTSim)
library(abind)
library(reshape2)



# Define working directory
sum_dir = 'C:/Users/jrcoyle/Documents/Research/CT-Sim/GitHub/Results/Summary'
fig_dir = 'C:/Users/jrcoyle/Documents/Research/CT-Sim/GitHub/Results/Plots'
sim_dir = 'C:/Users/jrcoyle/Documents/Research/CT-Sim/GitHub/Code'

setwd(sum_dir)

# Load simulation functions
#source(file.path(sim_dir, 'simulation_functions.R'))


##########################
## Functions


# A function that extracts parameter values from a runID
#	runID = a string in the filename that gives the parameters
#	assoc_str = character pairing a parameter name with its value
#	sep_str = character separating different parameter-value pairs
get_parms = function(runID, assoc_str='-', sep_str='_'){
	parm_pairs = sapply(strsplit(runID, sep_str), function(x) strsplit(x, assoc_str))
	
	vals = sapply(parm_pairs, function(x) x[2])
	vals = as.data.frame(t(vals))
	names(vals) = sapply(parm_pairs, function(x) x[1])

	vals
}


make_plot = function(xlim, ylim, xlab=NULL, ylab=NULL, cex=1){	
	plot.new()
	plot.window(xlim=xlim, ylim=ylim)
	axis(1)
	abline(h=par('usr')[3], lwd=3)
	axis(2, las=1)
	abline(v=par('usr')[1], lwd=3)
	if(!is.null(xlab)) mtext(xlab, 1, 2.5, cex=cex)
	if(!is.null(ylab)) mtext(ylab, 2, 3, cex=cex)

}

##########################



##### EXPERIMENT 1 #####
## Effects of dispersal and landscape autocorrelation

## Parameter definitions

# Dispersal kernels: d is for propagule dispersal, v is for adult movement 
# a0 = adjacent cell dispersal with probability = 0 of moving
# a0.5 = adjacent cell dispersal with probability = 0.5 of moving
# a1 = adjacent cell dispersal with probability = 1 of moving
# g<D> = gaussian dispersal with probability = 0.99 that dispersal distance is <= D
# u45 = uniform dispersal within a radius of 45 (unlimited dispersal on the 32 x 32 grid)

# Landscape autocorrelation
# dcorr: number gives the range of the variogram model, the distance at which habitat values are uncorrelated


setwd(file.path(sum_dir, 'EXP1'))

# Get list of runs
runlist = list.dirs('.', full.names=F)

# Define objects to hold data
bio = data.frame() 
occ = data.frame()

# Define subset of data to focus on
use_subset = expression((cross_run %in% c('50%','2.5%','97.5%'))&(cross_space %in% c('mean','var')))

# Go through parameter sets
for(d in runlist){
	parm_vals = get_parms(d)
	
	# Go through summaries at different scales
	for(f in list.files(d)){

		# Define scale
		find_scale = regexpr('L([1-5])_', f, perl=TRUE)
		scale = 2^(as.numeric(substr(f, attr(find_scale, 'capture.start'), attr(find_scale, 'capture.start') + attr(find_scale, 'capture.length') - 1))-1)
		
		# Adds two data objects to environment: sim_sum and sim_sum_ind
		# Both are lists and sim_sum summarizes across runs, whereas sim_sum_ind contains data from each of the 100 runs
		load(file.path(d,f))

		# Extract landscape heterogeneity
		land_het = sim_sum$land['mean','het']

		# Get stats for biologically-based categories
		dat_bio = melt(sim_sum$bio)
		dat_bio = do.call('subset', list(x=dat_bio, subset=use_subset))
		dat_bio = dcast(dat_bio, cross_space + cross_run + p_obs ~ category + comm_stat)
		dat_bio$scale = scale	
		dat_bio$land_het = land_het
		bio = rbind(bio, cbind(dat_bio, parm_vals))
		
		# Get data of occupancy distribution
		dat_occ = melt(sim_sum$occ)
		dat_occ = do.call('subset', list(x=dat_occ, subset=use_subset))
		dat_occ = dcast(dat_occ, cross_space + cross_run + p_obs + comm_stat ~ category)
		dat_occ$scale = scale
		dat_occ$land_het = land_het
		occ = rbind(occ, cbind(dat_occ, parm_vals))

		# For this run it doesn't make sense to look at xclass since it only tallies
		# the species in the first and last occupancy classes (<0.1 and >0.9) as
		# transient or core
	
	}
}


# Save summary files as csv
write.csv(bio, 'EXP1_bio.csv', row.names=F)
write.csv(occ, 'EXP1_occ.csv', row.names=F)

# Load saved files
bio = read.csv(file.path(sum_dir,'EXP1','EXP1_bio.csv'), check.names=F)
occ = read.csv(file.path(sum_dir,'EXP1','EXP1_occ.csv'), check.names=F)

# Set working directory to figure directory
setwd(file.path(fig_dir, 'EXP1'))



# Define occupancy categories
xvals = seq(0.05, 1, 0.05)

# Define dispersal kernels
dkerns = c('a0','a0.5','a1','g1','g2','g4','g8','g16','g32','u45')

# Define landscape auto correlation
dcorrs = 2^(0:4)

# Define scales
scales = 2^(0:4)

# Define detection probabilities
pobs = seq(0.1, 1, .1)

# Define print symbols
use_pch = c(16,1)

## Plots focusing on complete detection
use_P = 1


## Species occupancy distributions
dat = subset(occ, p_obs==use_P & cross_space=='mean' & comm_stat=='rich')

pdf('EXP1_occupancy_distributions_P=1.pdf', height=10, width=10)
par(mfrow=c(length(scales), length(dcorrs)))
par(mar=c(3, 3, 1, 0))
par(oma = c(3, 8, 3, 1))

# Dispersal kernels on different pages
for(k in dkerns){

# Spatial grain and range on grid
for(sp in scales){
for(ac in dcorrs){


	this_dat = subset(dat, d==k & v==k & scale==sp & dcorr==ac)

	# CIs are across runs
	low95 = subset(this_dat, cross_run=='2.5%')[,5:24]
	up95 = subset(this_dat, cross_run=='97.5%')[,5:24]
	med = subset(this_dat, cross_run=='50%')[,5:24]

	make_plot(c(0,1), c(0,max(up95)))

	polygon( c(xvals, rev(xvals)), c(low95, rev(up95)), border=NA, col='#00000050')
	lines(xvals, med)

	if(sp==scales[1]) mtext(paste('Range =', ac), 3, 0.5, cex=0.8)
	if(ac==dcorrs[1]){
		par(xpd = NA)
		text(-.5, par('usr')[4]*.5, paste('Grain =', sp), cex=1.3, pos=2)
		par(xpd = F)
		mtext('Num. Species', 2, 2.5, cex=0.8)		
	}

}}

mtext(paste('Dispersal =', k), 3, 1, outer=T)
mtext('Temporal Occupancy', 1, 1, outer=T)

}
dev.off()


## Occupancy distributions when juvenile dispersal and adult movement differ

pdf('EXP1_occupancy_distributions_P=1_dcorr=4.pdf', height=10, width=8)
par(mfrow=c(length(scales), 3))
par(mar=c(3, 3, 1, 0))
par(oma = c(3, 8, 3, 1))

# Juvenile dispersal kernels on different pages
for(k in c('g1','g4','g16','u45')){

# Spatial grain and adult movement on grid
for(sp in scales){
for(mv in c('a0','a0.5','a1')){
	
	this_dat = subset(dat, d==k & v==mv & scale==sp & dcorr==4)

	# CIs are across runs
	low95 = subset(this_dat, cross_run=='2.5%')[,5:24]
	up95 = subset(this_dat, cross_run=='97.5%')[,5:24]
	med = subset(this_dat, cross_run=='50%')[,5:24]

	make_plot(c(0,1), c(0,max(up95)))

	polygon( c(xvals, rev(xvals)), c(low95, rev(up95)), border=NA, col='#00000050')
	lines(xvals, med)

	if(sp==scales[1]) mtext(paste('Movement =', mv), 3, 0.5, cex=0.8)
	if(mv=='a0'){
		par(xpd = NA)
		text(-.5, par('usr')[4]*.5, paste('Grain =', sp), cex=1.3, pos=2)
		par(xpd = F)
		mtext('Num. Species', 2, 2.5, cex=0.8)		
	}

}}

mtext(paste('Dispersal =', k), 3, 1, outer=T)
mtext('Temporal Occupancy', 1, 1, outer=T)

}
dev.off()


## Occupancy distributions under different detection probs

use_dcorr = 4
dat = subset(occ, cross_space=='mean' & comm_stat=='rich' & dcorr==use_dcorr)

dat$scale = as.numeric(dat$scale)

# Spatial grain and detection probability on grid
pdf('EXP1_occupancy_distributions_detectionP_dcorr=4.pdf', height=12, width=10)
par(mfrow=c(8, length(scales)))
par(mar=c(3, 3, 1, 0))
par(oma = c(3, 8, 3, 1))

# Dispersal kernels on separate pages
for(k in dkerns){

# Spatial grain and adult movement on grid
# FOR SOME REASON, CAN'T SUBSET OF p = 0.7 OR p = 0.3
for(p in c(1, .9, .8, .6, .5, .4, .2, .1)){
for(sp in scales){

	this_dat = subset(dat, scale==sp & d==k & v==k & p_obs==p)

	# CIs are across runs
	low95 = subset(this_dat, cross_run=='2.5%')[,5:24]
	up95 = subset(this_dat, cross_run=='97.5%')[,5:24]
	med = subset(this_dat, cross_run=='50%')[,5:24]

	make_plot(c(0,1), c(0,max(up95)))

	polygon( c(xvals, rev(xvals)), c(low95, rev(up95)), border=NA, col='#00000050')
	lines(xvals, med)

	if(p==1) mtext(paste('Grain =', sp), 3, 0.5, cex=0.8)
	if(sp==scales[1]){
		par(xpd = NA)
		text(-.5, par('usr')[4]*.5, paste('P =', p), cex=1.3, pos=2)
		par(xpd = F)
		mtext('Num. Species', 2, 2.5, cex=0.8)		
	}

}}

mtext(paste('Dispersal =', k), 3, 1, outer=T)
mtext('Temporal Occupancy', 1, 1, outer=T)

}
dev.off()




## Biological designations
dat = subset(bio, p_obs==use_P & cross_space=='mean')

jit = 0.1
xvals = dcorrs

pdf('EXP1_bio_rich-abun_P=1.pdf', height=9, width=6)
par(mfrow=c(length(scales),2))
par(mar=c(3, 5, 1, 0))
par(oma = c(3, 8, 3, 1))

for(k in dkerns){
for(sp in scales){

	this_dat = subset(dat, d==k & v==k & scale==sp)
	this_dat =  this_dat[order(this_dat$dcorr),]
	
	low95 = subset(this_dat, cross_run=='2.5%')
	up95 = subset(this_dat, cross_run=='97.5%')
	med = subset(this_dat, cross_run=='50%')

	# Richness
	make_plot(c(0,max(dcorrs)), c(0, max(up95[,grep('rich', names(up95))])), ylab='Num. Species', cex=0.8)
	
	segments(xvals+jit, low95$core_rich, xvals+jit, up95$core_rich)
	segments(xvals-jit, low95$trans_rich, xvals-jit, up95$trans_rich)

	points(xvals+jit, med$core_rich, pch=use_pch[1])
	points(xvals-jit, med$trans_rich, pch=use_pch[2])
	
	# Add x label
	if(sp==scales[length(scales)]) mtext('Range', 1, 2.5, cex=0.8)

	# Add scale label
	par(xpd = NA)
	text(-7, par('usr')[4]*.5, paste('Grain =', sp), cex=1.3, pos=2)
	par(xpd = F)

	# Abundance
	make_plot(c(0,max(dcorrs)), c(0, max(up95[,grep('abun', names(up95))])), ylab='Num. Individuals', cex=0.8)
	
	segments(xvals+jit, low95$core_abun, xvals+jit, up95$core_abun)
	segments(xvals-jit, low95$trans_abun, xvals-jit, up95$trans_abun)

	points(xvals+jit, med$core_abun, pch=use_pch[1])
	points(xvals-jit, med$trans_abun, pch=use_pch[2])

	# Add legend
	if(sp==scales[1]){
		par(xpd = NA)
		legend(par('usr')[2], 1.2*par('usr')[4], c('Core','Transient'), pch=use_pch, xjust=1, yjust=0, horiz=T)
		par(xpd = F)
	}

	# Add x label
	if(sp==scales[length(scales)]) mtext('Range', 1, 2.5, cex=0.8)
}

# Add dispersal label
mtext(paste('Dispersal =', k), 3, 1, outer=T)
}
dev.off()






###########################################################################################

##### Runs testing convergence ######

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


sumID = 'converge32_b-2'

# Define times and scales
endTs = seq(25, 1000, 25)
scales = paste0('L',1:5)

# Read in summaries through time
for(i in 1:length(scales)){

	load(file.path(sumID, paste0(scales[i], '-T', endTs[1], '_summary.RData')))
	this_bio = sim_sum_ind$bio
	this_occ = sim_sum_ind$occ
	this_xclass = sim_sum_ind$xclass
	this_abun = sim_sum_ind$abun

	for(j in 2:length(endTs)){
	
		this_endT = endTs[j]	
	
		this_file = file.path(sumID, paste0(scales[i], '-', 'T', this_endT,'_summary.RData'))
		load(this_file)

		this_bio = abind(this_bio, sim_sum_ind$bio, along=ifelse(j==2, 0, 1))
		this_occ = abind(this_occ, sim_sum_ind$occ, along=ifelse(j==2, 0, 1))
		this_xclass = abind(this_xclass, sim_sum_ind$xclass, along=ifelse(j==2, 0, 1))
		this_abun = abind(this_abun, sim_sum_ind$abun, along=ifelse(j==2, 0, 1))	
	}
	
	if(i==1){
		bio_ind_T = this_bio
		occ_ind_T = this_occ
		xclass_ind_T = this_xclass
		abun_ind_T = this_abun

	} else {

		bio_ind_T = abind(bio_ind_T, this_bio, along=ifelse(i==2, 0, 1))
		occ_ind_T = abind(occ_ind_T, this_occ, along=ifelse(i==2, 0, 1))
		xclass_ind_T = abind(xclass_ind_T, this_xclass, along=ifelse(i==2, 0, 1))	
		abun_ind_T = abind(abun_ind_T, this_abun, along=ifelse(i==2, 0, 1))	
	}
}
	
# Assign dimension names
dimnames(bio_ind_T)[[1]] = 2^(0:4)
dimnames(bio_ind_T)[[2]] = endTs
dimnames(occ_ind_T)[[1]] = 2^(0:4)
dimnames(occ_ind_T)[[2]] = endTs
dimnames(xclass_ind_T)[[1]] = 2^(0:4)
dimnames(xclass_ind_T)[[2]] = endTs
dimnames(abun_ind_T)[[1]] = 2^(0:4)
dimnames(abun_ind_T)[[2]] = endTs
names(dimnames(bio_ind_T)) = c('scale', 'endT', 'run', 'comm_stat','cross_space','category','p_obs')
names(dimnames(occ_ind_T)) = c('scale', 'endT', 'run', 'comm_stat','cross_space','category','p_obs')
names(dimnames(xclass_ind_T)) = c('scale', 'endT', 'run','cross_space','ab_ct','p_obs')
names(dimnames(abun_ind_T)) = c('scale', 'endT', 'run','sp_rank','p_obs')

# Save objects for later comparison of dispersal method

#bio_dunif = bio_ind_T; occ_dunif = occ_ind_T; xclass_dunif = xclass_ind_T
#bio_d9 = bio_ind_T; occ_d9 = occ_ind_T; xclass_d9 = xclass_ind_T
#bio_d4 = bio_ind_T; occ_d4 = occ_ind_T; xclass_d4 = xclass_ind_T

# Scale on different pages
# Community statistic on different pages
# Detectability in different plots

nruns=4
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


## Species relative (ranked) abundance distributions
nspecies = dim(abun_ind_T)[[4]]
lcol = rainbow(nspecies)

pdf(file.path(sum_dir, paste0(sumID, '_ind_run_SAD_through_time.pdf')), height=12, width=10)

layout(matrix(1:20, nrow=5, byrow=T))
par(mar=c(2, 2, 1, 1))
par(oma=c(3,3,3,0))
par(lend=1)

for(p in rev(p_obs)){
	for(sp in 2^(0:4)){
		maxS = max(abun_ind_T[as.character(sp), , , , as.character(p)])

		for(i in 1:nruns){
			make_plot(c(0, 1000), c(0, maxS), xlab=ifelse(sp==16, 'Time', ''), ylab=ifelse(i==1, 'Relative abun.',''))
			text(0, 0.99*maxS, paste('Grain =',sp,'x',sp), pos=4)
			text(1000, 0.99*maxS, paste('Run', i), pos=2)
	
			for(j in 1:nspecies){
				lines(endTs, abun_ind_T[as.character(sp), , i, j, as.character(p)], lwd=1, col=lcol[j])
			}
		}
	}

	mtext(paste('Detection prob. =', p), 3, 1, outer=T)
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

sumID = 'converge32_d-adj'

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










##### Examine species abundance distributions
setwd('C:/Users/jrcoyle/Documents/Research/CT-Sim/Runs/')

this_run = 'converge32_d-9/converge32_run1.RData'

load(this_run)

mymeta = results[,,1]


sp_rank = calc_abun(mymeta, 40, only_species=T, ranked=T)

















