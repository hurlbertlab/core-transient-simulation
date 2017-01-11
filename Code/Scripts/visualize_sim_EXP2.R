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
data_dir = 'Z:/CTSim/Data'

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
pixDyn = function(results, this_land, row, col, lab = NULL, timewindow = NULL, scale = 3) {
  
  if(this_land[row, col] == 1) { hab = 'A' } else { hab = 'B' }
  
  if(!is.null(lab)) { lab = paste(lab, "; ", sep = '') }
  
  # Fraction of identical landscape over (2*scale+1)x(2*scale+1) region
  het = sum(this_land[max(row-scale, 1):min(row+scale, 32), max(col-scale, 1):min(col+scale, 32)] == this_land[row, col])/
    length(this_land[max(row-scale, 1):min(row+scale, 32), max(col-scale, 1):min(col+scale, 32)])
  
  # Species #1-20 are by definition core in Habitat B; species 21-40 in Habitat A
  if (hab == 'B') {
    core = unlist(lapply(results[row, col, ], function(x) sum(unique(x) <= 20)))
    tran = unlist(lapply(results[row, col, ], function(x) sum(unique(x) > 20)))
  } else {
    core = unlist(lapply(results[row, col, ], function(x) sum(unique(x) > 20)))
    tran = unlist(lapply(results[row, col, ], function(x) sum(unique(x) <= 20)))
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
      uniqsp = unique(unlist(results[row, col, ((t-1)*timewindow + 2):(t*timewindow+1)]))
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

# Plot occupancy histogram for the final time window 162:201 (time 161:200)
pixOccHist = function(results, row, col, lab, timewindow = 40, binwidth = 4) {
  
  tmp = results[row, col, (202 - timewindow):201]
  unq = lapply(tmp, function(x) unique(x))
  occs = table(unlist(unq))
  
  # Define histogram breaks
  xrange = seq(0, timewindow, binwidth)
  
  # Split into two groups (species ID <= or > 20)
  ocore = occs[as.numeric(names(occs)) <= 20]
  otran = occs[as.numeric(names(occs)) > 20]
  
  # compute the counts per interval
  hv1 = hist(ocore,breaks=xrange,plot=F)$counts
  hv2 = hist(otran,breaks=xrange,plot=F)$counts
  
  # Generate a a stacked histogram
  barplot(rbind(hv1,hv2), col=c("#FF0000FF", "#FFFFBFFF"), names.arg = xrange[-1]/timewindow,
          space = 0, las = 1, ylab = "Number of species", xlab = "Occupancy",
          main = paste("Pixel", lab))
}


# Generic (fixed pixel location) dynamics for any simulation run
# data_dir : directory where raw simulation results are stored
# sim      : sim name specifying parameter combinations, e.g. hp-0.9
# run      : simulation run #
# plot_dir : directory for saving output plot
pixelSummary = function(data_dir, sim, run=1:20, plot_dir, plot.pdf = TRUE) {
    
  if (plot.pdf) {
    pdf(paste(plot_dir, '/', sim, '_run', run[1],'-',run[length(run)], '_dynamics.pdf', sep = ''), 
        height = 10, width = 8)
  }
  par(mfrow = c(4,3))
  
  for (r in run) {
    suppressWarnings(rm(list = c('results', 'res', 'this_land', 'this_metacomm', 'this_species')))
    load(paste(data_dir, '/', sim, '_run', r, '.Rdata', sep = ''))
    
    # Results from early sims (pre-turnover) are in slightly different structure
    if (class(results) == 'array') {
      res = results
    } else {
      res = results$sim
    }
    
    image(this_land, main = paste(sim, '_run', r, sep = ''))
    
    # Gridded pixels for investigation
    sites = data.frame(id = c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K'), 
                       row = c(5, 5, 5, 16, 16, 16, 16, 28, 28, 28, 28), 
                       col = c(5, 16, 27, 4, 12, 20, 28, 4, 12, 20, 28))
    
    text(sites$col, 33-sites$row, sites$id, cex = .5)
    
    sapply(1:nrow(sites), function(x) pixDyn(res, this_land, sites$row[x], sites$col[x], sites$id[x], timewindow = 40))
    
    image(this_land, main = paste(sim, '_run', r, sep = ''))
    text(sites$col, 33-sites$row, sites$id, cex = .5)
    sapply(1:nrow(sites), function(x) 
      pixOccHist(res, sites$row[x], sites$col[x], sites$id[x], timewindow = 40, binwidth = 4))
    
  }
  
  if (plot.pdf) {
    dev.off()
  }
}

# Calculate similarity of landscape to focal pixel
# (Fraction of identical landscape over (2*scale+1)x(2*scale+1) region)

#   loc.xy is a matrix of x and y coordinates on the landscape grid
#   land is the landscape grid raster
#   scale is the (square) radius over which similarity is calculated
land_similarity = function(loc.xy, land, scale) {
  het = apply(loc.xy, 1, function(i) {
    sum(land[max(i[1]-scale, 1):min(i[1]+scale, 32), max(i[2]-scale, 1):min(i[2]+scale, 32)] == land[i[1], i[2]])/
      length(land[max(i[1]-scale, 1):min(i[1]+scale, 32), max(i[2]-scale, 1):min(i[2]+scale, 32)])
  })
  return(het)
}


# Generic pixel cross-classification analysis for all (non-edge) pixels 
# data_dir : directory where raw simulation results are stored
# sim      : sim name specifying parameter combinations, e.g. hp-0.9
# run      : simulation run #
# scale    : radius in pixels over which landscape heterogeneity is calculated
pixelXclass = function(data_dir, sim, run=1, scale = 3, t_window = 186:200, 
                       ct_threshold = 1/3, return = 'percent') {
  
  timewindow = length(t_window)
  
  suppressWarnings(rm(list = c('results', 'res', 'this_land', 'this_metacomm', 'this_species')))
  load(paste(data_dir, '/', sim, '_run', run, '.Rdata', sep = ''))
  
  # Results from early sims (pre-turnover) are in slightly different structure
  if (class(results) == 'array') {
    res = results
  } else {
    res = results$sim
  }
  
  # Check that required objects exist
  if(!exists('this_species')|!exists('this_land')|!exists('this_gsad')){
    stop(paste(sim, 'may not contain this_species, this_land, or this_gsad. Summary NOT run.'))		
  } else {
    # Assign objects
    species = this_species
    land = this_land
    gsad = this_gsad
  }
  
  # Number of species
  N_S = dim(species)[1]
  
  # Define different scales of spatial aggregation (in a partition of the grid)
  # Must supply grid dimensions
  scale_locs = sapply(2^c(0:4), function(fact) aggregate_cells(X=c(32,32), dX=fact, dY=fact, form='partition'))
  
  # Locations to be aggregated and evaluated
  locs = scale_locs[[1]]
  
  # Calculate species abundance profiles at the spatial and temporal resolution given by locs and t_window
  abuns_act = calc_abun_profile(locs, t_window, res, N_S)
  
  # Range of detection probabilities to examine
  P_obs = seq(0.1, 1, by = .1)
  
  # Apply observation bias
  abuns_obs = sapply(P_obs, function(p){
    sample_sim(abuns_act, probs = p, return='abundance')
  }, simplify='array')
  dimnames(abuns_obs)[[2]] = 1:N_S # Name columns with species names
  # dims are now [timepoint, species, spatial unit, P]
  
  # Calculate species occupancy
  occ = sapply(1:length(P_obs), function(i){
    use_abun = abuns_obs[,,,i, drop=FALSE]
    use_abun = abind::adrop(use_abun, 4) 
    calc_occupancy(abuns=use_abun, agg_times=NULL, do_freq=F)
  }, simplify='array')
  # dims are now: [spatial unit, species, P]
  
  breaks = seq(ct_threshold, 0.9999, by = ct_threshold)
  
  # Get species birth rates
  b_rates = species[,,'b']
  
  # Get habitat types for spatial units
  habitats = sapply(locs, function(x) average_habitat(x, land))
  
  # Calculate classification of species based on birth rates
  cores = t(sapply(habitats, function(h) b_rates[,h]>0))
  classification = apply(cores, 1:2, function(x) ifelse(x, 'core', 'trans'))	
  
  # Calculate proportion mis-classified
  xclass = sapply(1:length(P_obs), function(i){
    use_occ = occ[,,i, drop=FALSE]
    use_occ = abind::adrop(use_occ, 3)
    tabs = cross_classify(use_occ, breaks, classification=classification, do_each=T, return='counts')
    reshape2::acast(reshape2::melt(tabs, varnames=c('bio','occ','sp_unit')), sp_unit ~ bio+occ)
  }, simplify='array')
  # dims are now [spatial unit, category, P]
  
  xclass.pct = array(dim = dim(xclass))
  for (d in 1:dim(xclass)[1]) {
    co = xclass[d,1:2,] #biologically core spp classified as either core or transient
    tr = xclass[d,3:4,] #biologically transient spp classified as either core or transient
    corepct = sapply(1:ncol(co), function(i) {
      if (colSums(co)[i] == 0) { 
        c(NA, NA) 
      } else { 
        co[,i]/matrix(colSums(co)[i], nrow=1)
      }
    })
    
    tranpct = sapply(1:ncol(tr), function(i) {
      if (colSums(tr)[i] == 0) { 
        c(NA, NA) 
      } else { 
        tr[,i]/matrix(colSums(tr)[i], nrow=1)
      }
    })
    
    xclass.pct[d,,] = rbind(corepct, tranpct)
  }

  loc.xy = matrix(unlist(locs), ncol = 2, byrow = TRUE)
  
  landsim = land_similarity(loc.xy, land, scale)
  
  landscapeS = data.frame(x = loc.xy[,1], y = loc.xy[,2], sim = landsim)
  
  if(return == 'count') {
    return(list(xclass=xclass, landSim = landscapeS))
  } else if (return == 'percent') {
    return(list(xclass = xclass.pct, landSim = landscapeS))
  }
}

run1 = r1$xclass[r1$landSim$sim > 2/3 & 
                   r1$landSim$x > scale & r1$landSim$x < max(r1$landSim) - scale &
                   r1$landSim$y > scale & r1$landSim$y < max(r1$landSim) - scale ,,]

datadir = 'z:/lab/CTSim/Data/Exp3/d-g2_imm-0.001'

xclass.out = c()
land.out = c()
for (r in 1:50) {
  temp = pixelXclass(datadir, 'd-g2_imm-0.001', run=r, scale = 3, t_window = 186:200, 
                         ct_threshold = 1/3, return = 'percent')
  
  # Eliminate pixels within a distance 'scale' from the grid edge
  temp2 = temp$xclass[#temp$landSim$sim > 2/3 & 
              temp$landSim$x > scale & temp$landSim$x < max(temp$landSim) - scale &
              temp$landSim$y > scale & temp$landSim$y < max(temp$landSim) - scale ,,]
  xclass.out = abind(xclass.out, temp2, along = 1)
  
  land2 = temp$landSim[#temp$landSim$sim > 2/3 & 
                         temp$landSim$x > scale & temp$landSim$x < max(temp$landSim) - scale &
                         temp$landSim$y > scale & temp$landSim$y < max(temp$landSim) - scale ,,]
  land.out = rbind(land.out, land2)
  
  print(paste(r, Sys.time()))
}

# Classification error as a function of detection probability
#   (only for pixels with heterogeneity/similarity > 2/3)
means = apply(xclass.out[land.out$sim > 2/3 , , ], c(2,3), function(x) mean(x, na.rm = T))
P_obs = seq(0.1, 1, by = .1)

par(mar = c(4,4,1,4), mgp = c(3, 1, 0))
plot(P_obs, 100*means[1,], type = 'l', lwd = 4, ylim = c(0,100), 
     xlab = "Detection probability", ylab = "Percent", col = 'cornflowerblue')
axis(4)
abline(h=90, lty = 'dotted')
points(P_obs, 100*means[4,], type = 'l', lwd = 4, col = 'salmon')
points(P_obs, 100*means[2,], type = 'l', lwd = 4, lty = 'dashed', col = 'cornflowerblue')
points(P_obs, 100*means[3,], type = 'l', lwd = 4, lty = 'dashed', col = 'salmon')
legend(0.8, 70, c('C-C', 'T-T', 'C-T', 'T-C'), col = rep(c('cornflowerblue', 'salmon'),2),
       lwd = 4, lty = rep(c('solid', 'dashed'), each = 2), cex = 1.5)

# Classification error as a function of landscape similarity
xc.p1 = 100*xclass.out[, 1:2, 10]

plot(land.out$sim, xc.p1[,1], pch = 16, xlab = "Landscape similarity", 
     ylab = "Percent", col = 'cornflowerblue')
points(land.out$sim, xc.p1[,2], pch = 17, col = 'salmon')







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

sapply(1:nrow(sites), function(x) pixDyn(sites$row[x], sites$col[x], sites$id[x], timewindow = 40))

image(this_land, main = 'hp-0.9_run1')
text(sites$col, 33-sites$row, sites$id, cex = .5)
sapply(1:nrow(sites), function(x) 
  pixOccHist(sites$row[x], sites$col[x], sites$id[x], timewindow = 40, binwidth = 4))

dev.off()

# pixelSummary(data_dir = 'Z:/Lab/CTSim/Data/EXP2/hp-0.5', sim = 'hp-0.5', plot_dir = 'Results/Plots/EXP2')


## Checking summary functions manually

# Define locations as single cells given above
pix_locs = aggregate_cells(32, 32, 1, 1, form='partition')
site_locs = sapply(1:nrow(sites), function(i) list(sites[i, c('row','col')]))

# Calculate abundance profiles
site_abun_profs = calc_abun_profile(site_locs, 161:200, results$sim, 40) 
apply(site_abun_profs, c(2,3), sum)

# Calculate species temporal occupancy
site_occs = calc_occupancy(abuns=site_abun_profs)
site_richCT = calc_rich_CT(site_abun_profs, site_occs, c(1/3, 2/3))

# Calculate classification of species based on birth rates
habitats = sapply(site_locs, function(x) average_habitat(x, this_land))
b_rates = this_species[,,'b']
cores = t(sapply(habitats, function(h) b_rates[,h]>0))
classification = apply(cores, 1:2, function(x) ifelse(x, 'core', 'trans'))	
occ_ab = apply(classification, 1:2, function(x) ifelse(x=='core', 1, .1))
site_richAB = calc_rich_CT(site_abun_profs, occ_ab, 0.5)

# Examine a given site:
use_site=1 # C
site_abun_profs[,,use_site] # abundance profiles of each species over last 40 timesteps
site_occs[use_site,] # temporal occupancy of each species over last 40 timesteps
sum(site_occs[use_site,]<1/3); sum(site_occs[use_site,]>2/3) # number of species classified as trans / core
site_richCT[,,use_site] # at each timestep, number of core/trans species
site_richAB[,,use_site]

# Do summary using summarize function, averaging across all pixels
site_sum = summarize_sim(sim=results$sim, breaks = c(1/3,2/3), locs=pix_locs, 
	t_window=list(start=161, stop=200), species=this_species, land=this_land, gsad=this_gsad, 
	sum_parms=list(time_sum='last'), sum_turn=F)


# Summarizing when species richness is accumulated across a wider time window 
site_sum_10 = summarize_sim(sim=results$sim, breaks = c(1/3,2/3), locs=pix_locs, 
	t_window=list(start=161, stop=200), species=this_species, land=this_land, gsad=this_gsad, 
	sum_parms=list(time_sum='last', agg_times=10), sum_turn=F)

# Compare
site_sum$occ['rich','mean',,]
site_sum_10$occ['rich','mean',,]




