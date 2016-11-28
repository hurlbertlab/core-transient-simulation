## This script is used to visualize simulation summaries from EXP 2

# in which habitat proportion (of habitat type A) was varied from 0.5 to 1

# Set options and load libraries
options(stringsAsFactors=F)
library(CTSim)
library(abind)
library(reshape2)

# Define working directory
sum_dir = 'Results/Summary/EXP2'
fig_dir = 'Results/Plots/EXP2'
sim_dir = 'Code'
data_dir = 'Z:/Lab/CTSim/Data'

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


# Function to plot core/transient dynamics for a pixel
pixdyn = function(row, col, lab = NULL, timewindow = NULL, scale = 3) {
  
  if(this_land[row, col] == 1) { hab = 'A' } else { hab = 'B' }
  
  if(!is.null(lab)) { lab = paste(lab, "; ", sep = '') }
  
  # Fraction of identical landscape over (2*scale+1)x(2*scale+1) region
  het = sum(this_land[max(row-scale, 1):min(row+scale, 32), max(col-scale, 0):min(col+scale, 32)] == this_land[row, col])/
    length(this_land[max(row-scale, 1):min(row+scale, 32), max(col-scale, 0):min(col+scale, 32)])
  
  # Species #1-20 are by definition core in Habitat B; species 21-40 in Habitat A
  if (hab == 'B') {
    core = unlist(lapply(results$sim[row, col, ], function(x) sum(unique(x) <= 20)))
    tran = unlist(lapply(results$sim[row, col, ], function(x) sum(unique(x) > 20)))
  } else {
    core = unlist(lapply(results$sim[row, col, ], function(x) sum(unique(x) > 20)))
    tran = unlist(lapply(results$sim[row, col, ], function(x) sum(unique(x) <= 20)))
  }
  plot(core, type = 'l', xlab = 'Time', ylab = 'Number of species', col = 'skyblue',
       main = paste(lab, "Landscape similarity ", round(het,2), ";\ncore (blue), transient (red)", sep = ''), 
       lwd = 2, ylim = c(0, max(c(core, tran))))
  points(tran, type = 'l', col = 'red', lwd = 2)
  
  # Optionally plot % transient species aggregated over timewindow
  if(!is.null(timewindow)) {
    pct.trans = c()
    times = timewindow*1:floor(200/timewindow) - timewindow/2
    for (t in 1:floor(200/timewindow)) {
      uniqsp = unique(unlist(results$sim[row, col, ((t-1)*timewindow + 2):(t*timewindow+1)]))
      if (hab == 'A') {
        pct.trans = c(pct.trans, sum(uniqsp <= 20)/length(uniqsp))
      } else {
        pct.trans = c(pct.trans, sum(uniqsp > 20)/length(uniqsp))
      }
    }
    par(new=T)
    plot(times, pct.trans, xlim = c(0,200), ylim = c(0,1), xlab = '', ylab = '', yaxt = 'n',
         xaxt = 'n', type = 'l')
    axis(4, at = seq(0,1, by = 0.25), tcl = .3, labels = F)
    mtext(c(1.0, 0.5, 0), 4, at = c(1, .5, 0), cex = .75)
  }
}
##########################

##### EXPERIMENT 2 #####
## Effects of habitat proportion

# Get list of runs
runlist = list.dirs(paste(sum_dir, '/.', sep =''), full.names=F)
runlist = runlist[runlist!='']

### Examine run summary files
# Define objects to hold data
bio = data.frame() 
occ = data.frame()
xclass = data.frame()

# Define subset of data to focus on
use_subset = expression((cross_run %in% c('50%','2.5%','97.5%'))&(cross_space %in% c('mean','var')))

# Go through parameter sets
for(d in runlist){
  parm_vals = get_parms(d)
  
    # Adds two data objects to environment: sim_sum and sim_sum_ind
    # Both are lists and sim_sum summarizes across runs, whereas sim_sum_ind contains data from each of the 100 runs
    load(file.path(paste(sum_dir, '/', d, sep = ''),list.files(paste(sum_dir, '/', d, sep = ''))))
    
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
    
    # Get data of occupancy distribution
    dat_xclass = melt(sim_sum$xclass)
    dat_xclass = do.call('subset', list(x=dat_xclass, subset=use_subset))
    dat_xclass = dcast(dat_xclass, cross_space + cross_run + p_obs ~ ab_ct)
    dat_xclass$land_het = land_het
    xclass = rbind(xclass, cbind(dat_xclass, parm_vals))

}

# Save summary files as csv
write.csv(bio, file.path(sum_dir, 'EXP2_bio.csv'), row.names=F)
write.csv(occ, file.path(sum_dir, 'EXP2_occ.csv'), row.names=F)
write.csv(xclass, file.path(sum_dir, 'EXP2_xclass.csv'), row.names=F)


# Load saved files
bio = read.csv(file.path(sum_dir, 'EXP2_bio.csv'), check.names=F)
occ = read.csv(file.path(sum_dir, 'EXP2_occ.csv'), check.names=F)


# Define occupancy categories
xvals = c(0.167, 0.5, 0.833)

# Define habitat proportions
hps = seq(0.5, 1, .1)

# Define detection probabilities (strangely, p = .6 not working)
pobs = c(0.2, 0.4, 0.5, 0.8, 1)

# Define print symbols
use_pch = c(16,1)

## Species occupancy distributions
dat = subset(occ, cross_space=='mean' & comm_stat=='rich')

pdf(file.path(fig_dir,'EXP2_occupancy_distributions.pdf'), height=10, width=10)
par(mfrow=c(length(pobs), length(hps)))
par(mar=c(3, 3, 1, 0))
par(oma = c(3, 8, 3, 1))

# Dispersal kernels on different pages
for(p in pobs){
  
  # Spatial grain and range on grid
  for(h in hps){
      
      this_dat = subset(dat, p_obs==p & hp==h)
      
      # CIs are across runs
      # Columns 5:7 are the occupancy classes from low to high, in this case only 3.
      # If more finely resolved 'breaks' vector specified in summary parameter file,
      # then columns will be 5:(5 + length(breaks) - 1)
      low95 = subset(this_dat, cross_run=='2.5%')[,5:7]
      up95 = subset(this_dat, cross_run=='97.5%')[,5:7]
      med = subset(this_dat, cross_run=='50%')[,5:7]
      
      make_plot(c(0,1), c(0,max(up95)))
      
      polygon( c(xvals, rev(xvals)), c(low95, rev(up95)), border=NA, col='#00000050')
      lines(xvals, med)
      
      if(p==pobs[1]) mtext(paste('HabA =', h), 3, 0.5, cex=0.8)
      if(h==hps[1]){
        par(xpd = NA)
        text(-.5, par('usr')[4]*.5, paste('P_obs =', p), cex=1.3, pos=2)
        par(xpd = F)
        mtext('Num. Species', 2, 2.5, cex=0.8)		
      }
      
    }
}
mtext('Temporal Occupancy', 1, 1, outer=T)
dev.off()

#--------------------------------------------------------------------------------
# Visualize pixel level dynamics


load(paste(data_dir, '/EXP2/hp-0.9/hp-0.9_run1.Rdata', sep = ''))

pdf('Results/Plots/EXP2/hp-0.9_run1_dynamics.pdf', height = 8, width = 8)
par(mfrow = c(3,3))
image(this_land, main = 'hp-0.9_run1')

# Pixels for investigation
sites = data.frame(id = c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'), 
                   row = c(5, 1, 5, 6, 17, 19, 19, 31), 
                   col = c(25, 9, 9, 12, 7, 8, 9, 26))

text(sites$col, 33-sites$row, sites$id, cex = .5)

sapply(1:nrow(sites), function(x) pixdyn(sites$row[x], sites$col[x], sites$id[x], timewindow = 40))

dev.off()

