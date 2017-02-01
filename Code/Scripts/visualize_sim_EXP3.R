## This script is used to visualize simulation summaries from EXP 3
#   (specifically where dispersal is set as a Gaussian (99% within 2 pixels),
#   and immigration to any pixel from the metacommunity is 0.001.)

# Set options and load libraries
options(stringsAsFactors=F)
library(CTSim)
library(dplyr)
library(tidyr)
library(abind)
library(reshape2)

# Define working directory
sum_dir = 'Results/Summary/EXP3'
fig_dir = 'Results/Plots/EXP3'
sim_dir = 'Code'
data_dir = 'z:/lab/CTSim/Data/Exp3/d-g2_imm-0.001'

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



# Generic pixel cross-classification analysis for all (non-edge) pixels 
# data_dir : directory where raw simulation results are stored
# sim      : sim name specifying parameter combinations, e.g. hp-0.9
# run      : simulation run #
# scale    : radius in pixels over which landscape heterogeneity is calculated
pixelXclassBySpecies = function(data_dir, sim, run=1, scale = 3, t_window = 186:200, 
                       ct_threshold = 1/3) {
  
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
  
  # Get total abundance of each species across the entire grid (remove 1st col describing empty cells)
  sad.grid = data.frame(sp = 1:40, Ngrid = apply(abuns_act[, -1, ], 2, sum))
  
  # Get abundance of each species in each pixel averaged over timewindow
  sad.cell = t(apply(abuns_act[, -1, ], 2:3, sum)) %>% data.frame() %>%
    mutate(pix = 1:1024) %>%
    gather("sp", "Ncell", 1:40) %>%
    mutate(sp = as.numeric(substr(sp, 2, nchar(sp)))) 
           
  
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
  
  # Flatten occupancy data
  occ2 = apply(occ, 2, I) %>% data.frame() %>%
    mutate(pix = rep(1:1024, times = 10),
           p = rep(seq(0.1, 1, 0.1), each = 1024)) %>%
    gather("sp", "occ", 1:40) %>%
    mutate(sp = as.numeric(substr(sp, 2, nchar(sp)))) 
  
  breaks = seq(ct_threshold, 0.9999, by = ct_threshold)
  
  # Get species birth rates
  b_rates = species[,,'b']
  
  # Get habitat types for spatial units
  habitats = sapply(locs, function(x) average_habitat(x, land))
  
  # Calculate classification of species based on birth rates
  cores = t(sapply(habitats, function(h) b_rates[,h]>0))
  classification = apply(cores, 1:2, function(x) ifelse(x, 'core', 'trans'))	
  
  # compare classification to occupancy-based designation
  xclass.sp = sapply(1:length(P_obs), function(i){
    use_occ = occ[,,i, drop=FALSE]
    use_occ = abind::adrop(use_occ, 3)
    xclass = classification
    xclass[classification == 'core' & use_occ > (1 - ct_threshold)] = 'c-c'
    xclass[classification == 'core' & use_occ <= ct_threshold & use_occ > 0] = 'c-t'
    xclass[classification == 'trans' & use_occ > (1 - ct_threshold)] = 't-c'
    xclass[classification == 'trans' & use_occ <= ct_threshold & use_occ > 0] = 't-t'
    xclass[xclass == 'core'] = 'c-int'  # biologically core but intermediate occupancy
    xclass[xclass =='trans'] = 't-int'  # biologically transient but intermediate occupancy
    xclass[use_occ == 0] = NA           # absent over the time period
    return(xclass)
  }, simplify='array')
  
  
  # Calculate landscape similarity to focal pixel
  loc.xy = matrix(unlist(locs), ncol = 2, byrow = TRUE)
  landsim = land_similarity(loc.xy, land, scale)
  landscapeS = data.frame(pix = 1:1024, x = loc.xy[,1], y = loc.xy[,2], sim = landsim)
  
  # xclass.sp flattened out into 2 dimensions 
  #   (rows increase in pixel.id first, then species, then P_obs)
  #   and then adding in landscape sim, abundances, etc
  xc2 = apply(xclass.sp, 2, I) %>% data.frame() %>% 
    mutate(pix = rep(1:1024, times = 10), 
           p = rep(P_obs, each = 1024)) %>% 
    gather("sp", "xc", 1:40) %>%
    mutate(sp = as.numeric(substr(sp, 2, nchar(sp)))) %>%
    full_join(landscapeS, by = "pix") %>%
    filter(x > scale & x < max(x) - scale & y > scale & y < max(y) - scale, 
           !is.na(xc)) %>%
    inner_join(sad.cell, by = c("pix", "sp")) %>%
    inner_join(sad.grid, by = "sp") %>%
    inner_join(occ2, by = c("pix", "p", "sp")) %>%
    arrange(pix, p, sp) %>%
    select(pix, sim, p, sp, Ncell, Ngrid, occ, xc)
    
  return(xclass = xc2)
}




#---------------------------------------------------------------------------------

xclass.out = c()
land.out = c()
xclass.ct = c()
scale = 3
for (r in 1:50) {
  temp = pixelXclass(datadir, 'd-g2_imm-0.001', run=r, scale = 3, t_window = 186:200, 
                        ct_threshold = 1/3, return = 'percent')
  tmp = pixelXclass(datadir, 'd-g2_imm-0.001', run=r, scale = 3, t_window = 186:200, 
                    ct_threshold = 1/3, return = 'count')
  
  # Eliminate pixels within a distance 'scale' from the grid edge
  temp2 = temp$xclass[temp$landSim$x > scale & temp$landSim$x < max(temp$landSim) - scale &
             temp$landSim$y > scale & temp$landSim$y < max(temp$landSim) - scale ,,]
  xclass.out = abind(xclass.out, temp2, along = 1)
  
  tmp2 = tmp$xclass[tmp$landSim$x > scale & tmp$landSim$x < max(tmp$landSim) - scale &
                      tmp$landSim$y > scale & tmp$landSim$y < max(tmp$landSim) - scale ,,]
  xclass.ct = abind(xclass.ct, tmp2, along = 1)
  
  
  land2 = temp$landSim[temp$landSim$x > scale & temp$landSim$x < max(temp$landSim) - scale &
                       temp$landSim$y > scale & temp$landSim$y < max(temp$landSim) - scale ,,]
  land.out = rbind(land.out, land2)
  
  print(paste(r, Sys.time()))
}


# Conduct cross-classification by species and pixel
#   (generates ~7M rows, may have memory issues)
xclass.sp = c()
for (r in 1:50) {
  tmp = pixelXclassBySpecies(data_dir, 'd-g2_imm-0.001', run=r, scale = 3, t_window = 186:200, 
                     ct_threshold = 1/3)
  xclass.sp = rbind(xclass.sp, tmp)
  rm(tmp)
  print(paste(r, Sys.time()))
}


save(xclass.out, xclass.ct, land.out, file = 'Results/Summary/EXP3/d-g2_imm-0.001/pixel_xclass_summary.Rdata')
save(xclass.sp, file = 'Results/Summary/EXP3/d-g2_imm-0.001/pixel_xclass_summary_bysp.Rdata')

# CAN START FROM HERE ONCE THE ABOVE HAS BEEN RUN ONCE:
# Load
load('Results/Summary/EXP3/d-g2_imm-0.001/pixel_xclass_summary.Rdata')


# Classification error as a function of detection probability
#   (only for pixels with heterogeneity/similarity > 2/3)
means = apply(xclass.out[land.out$sim > 2/3 , , ], c(2,3), function(x) 100*mean(x, na.rm = T))
ul95 = apply(xclass.out[land.out$sim > 2/3 , , ], c(2,3), function(x) 100*quantile(x, 0.975, na.rm = T))
ll95 = apply(xclass.out[land.out$sim > 2/3 , , ], c(2,3), function(x) 100*quantile(x, 0.025, na.rm = T))
P_obs = seq(0.1, 1, by = .1)


# Analysis based on counts (not %s) of species in each classification group

rel.error = c()
for (p in 1:10) {# loop over 10 detection probabilities
  tmp = xclass.ct[,,p]
  bio_core = tmp[,1] + tmp[,2]
  occ_core = tmp[,1] + tmp[,3]
  bio_tran = tmp[,4] + tmp[,3]
  occ_tran = tmp[,4] + tmp[,2]
  core.rel = occ_core/bio_core # estimated/actual core richness
  tran.rel = occ_tran/bio_tran # estimated/actual transient richness
  core.pct = tmp[,1]/occ_core # fraction of spp perceived to be core that actually are
  tran.pct = tmp[,4]/occ_tran # fraction of spp perceived to be transient that actually are
  core.dif = occ_core - bio_core # estimated - actual core richness
  tran.dif = occ_tran - bio_tran # estimated - actual transient richness
  rel.error = abind(rel.error, cbind(core.rel, tran.rel, core.pct, tran.pct, core.dif, tran.dif), along = 3)
}

# Mean values by detection probability for pixels with landscape similarity > 2/3
mean.rel = apply(rel.error[land.out$sim > 2/3 , , ], c(2,3), function(x) 
  mean(x[x!=Inf], na.rm = T))
ul95.rel = apply(rel.error[land.out$sim > 2/3 , , ], c(2,3), function(x) 
  quantile(x[x!=Inf], 0.975, na.rm = T))
ll95.rel = apply(rel.error[land.out$sim > 2/3 , , ], c(2,3), function(x) 
  quantile(x[x!=Inf], 0.025, na.rm = T))


# Classification error as a function of landscape similarity (when P_obs == 1)
xc2.p1 = 100*xclass.out[, , 10]
re.p1 = cbind(rel.error[, 1:2, 10], 100*rel.error[, 3:4, 10])



# Misclassification and over/underestimation as a function of detection probability

pdf('Results/Plots/EXP3/classification_errors_v_detect_prob.pdf', height = 10, width = 12)
par(mfrow = c(2,2), mar = c(5,6,3,4), mgp = c(3, 1, 0), 
    cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.25, las = 1)

# PANEL A) Classification rate based on true biological status
# Core classification
plot(P_obs, means[1,], type = 'n', lwd = 4, ylim = c(0,100), 
     xlab = "Detection probability", ylab = "Percent correctly classified",
     main = "% of spp whose occupancy matches their status")
axis(4)
polygon(c(P_obs, rev(P_obs)), c(ul95[1,], rev(ll95[1,])), col = rgb(0,0,1,0.1), border = NA)
points(P_obs, means[1,], type = 'l', lwd = 4, col = 'cornflowerblue') 

# Transient classification
polygon(c(P_obs, rev(P_obs)), c(ul95[4,], rev(ll95[4,])), col = rgb(1,0,0,0.1), border = NA)
points(P_obs, means[4,], type = 'l', lwd = 4, col = 'salmon')

abline(h=90, lty = 'dotted')
legend("bottomright", c('Core', 'Transient'), col = c('cornflowerblue', 'salmon'),
       lwd = 4, cex = 1.5)


# PANEL B) Degree of over/underestimation of true values by observed occupancy
plot(P_obs, mean.rel[1,], type = 'n', lwd = 4, ylim = c(0,3),
     xlab = "Detection probability", ylab = "Observed/True Richness",
     main = "Relative error")
axis(4)
polygon(c(P_obs, rev(P_obs)), c(ul95.rel[1,], rev(ll95.rel[1,])), col = rgb(0,0,1,0.1), border = NA)
points(P_obs, mean.rel[1,], type = 'l', lwd = 4, col = 'cornflowerblue') 

# Transient classificaiton
polygon(c(P_obs, rev(P_obs)), c(ul95.rel[2,], rev(ll95.rel[2,])), col = rgb(1,0,0,0.1), border = NA)
points(P_obs, mean.rel[2,], type = 'l', lwd = 4, col = 'salmon')

abline(h=1, lty = 'dotted')


# PANEL C) Classification rate based on occupancy class
plot(P_obs, 100*mean.rel[3,], type = 'n', lwd = 4, ylim = c(0,100),
     xlab = "Detection probability", ylab = "% correctly classified",
     main = "% of spp of a given occupancy class that are of that class")
axis(4)
polygon(c(P_obs, rev(P_obs)), c(100*ul95.rel[3,], rev(100*ll95.rel[3,])), col = rgb(0,0,1,0.1), border = NA)
points(P_obs, 100*mean.rel[3,], type = 'l', lwd = 4, col = 'cornflowerblue') 

# Transient classificaiton
polygon(c(P_obs, rev(P_obs)), c(100*ul95.rel[4,], rev(100*ll95.rel[4,])), col = rgb(1,0,0,0.1), border = NA)
points(P_obs, 100*mean.rel[4,], type = 'l', lwd = 4, col = 'salmon')


# PANEL D) Absolute degree of over/underestimation of richness in each class
par(mar = c(5,6,3,5))
plot(P_obs, mean.rel[6,], type = 'n', lwd = 4, ylim = c(-1,11),
     xlab = "Detection probability", ylab = "Overestimate of transient species",
     main = "Absolute error")
axis(4)
mtext("Underestimate of core species", 4, cex = 1.25, line = 3.5, las = 0)
polygon(c(P_obs, rev(P_obs)), c(ul95.rel[6,], rev(ll95.rel[6,])), col = rgb(0,0,0,0.1), border = NA)
points(P_obs, mean.rel[6,], type = 'l', lwd = 4, col = 'gray20') 

dev.off()


# Misclassification and over/underestimation as a function of landscape similarity

# eliminate NA and Inf values from transient vals

pdf('Results/Plots/EXP3/classification_errors_v_landscape.pdf', height = 15, width = 12)
par(mfrow = c(3, 2), mar = c(5,5,3,1), mgp = c(3, 1, 0), 
    cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.25, las = 1)
scatter.smooth(land.out$sim, xc2.p1[,1], degree = 2, xlab = 'Landscape similarity',
               ylab = 'Core % correct', col = 'cornflowerblue', pch = 16, 
               lpars = list(lwd = 3, col = 'darkblue'),
               main = "Fraction of core species that have high occupancy")
scatter.smooth(land.out$sim, xc2.p1[,4], degree = 2, xlab = 'Landscape similarity',
               ylab = 'Transient % correct', col = 'salmon', pch = 17, 
               lpars = list(lwd = 3, col = 'darkred'),
               main = "Fraction of transient species that have low occupancy")

scatter.smooth(land.out$sim[!re.p1[,3] %in% c(0, Inf)], 
               re.p1[!re.p1[,3] %in% c(0, Inf),3], degree = 2, xlab = 'Landscape similarity',
               ylab = 'Core % correct', col = 'cornflowerblue', pch = 16, 
               lpars = list(lwd = 3, col = 'darkblue'),
               main = "Fraction of high occupancy species that are actually core")
scatter.smooth(land.out$sim[!re.p1[,4] %in% c(0, Inf)], 
               re.p1[!re.p1[,4] %in% c(0, Inf),4], degree = 2, xlab = 'Landscape similarity',
               ylab = 'Transient % correct', col = 'salmon', pch = 17, 
               lpars = list(lwd = 3, col = 'darkred'),
               main = "Fraction of low occupancy species that are actually transient")

scatter.smooth(land.out$sim[!rel.error[,1,10] %in% c(0, Inf)], 
               rel.error[!rel.error[,1,10] %in% c(0, Inf),1,10], ylim = c(0,4),
               degree = 2, xlab = 'Landscape similarity',
               ylab = 'Obs/True Core Richness', col = 'cornflowerblue', pch = 16, 
               lpars = list(lwd = 3, col = 'darkblue'))
abline(h = 1, lty = 'dotted')
text(0.5, 4, "Repeated immigration by transients\n makes them appear to be core")

scatter.smooth(land.out$sim[!rel.error[,2,10] %in% c(Inf, 0) & !is.na(rel.error[,2,10])], 
               rel.error[!rel.error[,2,10] %in% c(Inf, 0) & !is.na(rel.error[,2,10]),2,10], 
               degree = 2, xlab = 'Landscape similarity', ylim = c(0,4),
               ylab = 'Obs/True Transient Richness', col = 'salmon', pch = 17, 
               lpars = list(lwd = 3, col = 'darkred'))
abline(h = 1, lty = 'dotted')
text(0.15, 1.5, "High # of true transients\nevery year due\n to dispersal")
text(0.6, 3.6, "Rare core spp show up\non connected landscapes and\nappear to be transient")

dev.off()


#----------------------------------------------------------------------------------
# Analysis at pixel/species level

load('Results/Summary/EXP3/d-g2_imm-0.001/pixel_xclass_summary_bysp.Rdata')

xc.sp.p1 = xclass.sp[xclass.sp$sim < 2/3 & xclass.sp$p == 1,]
xc.sp.p.5 = xclass.sp[xclass.sp$sim < 2/3 & xclass.sp$p == 0.5,]

# Number of species occurrences in each classification class
cat.counts.p1 = xc.sp.p1 %>% count(xc)
cat.counts.p.5 = xc.sp.p.5 %>% count(xc)

# Cross-classification %s
pct.bio_core1 = round(100*cat.counts.p1$n[1:3]/sum(cat.counts.p1$n[1:3]), 0)
pct.bio_tran1 = round(100*cat.counts.p1$n[4:6]/sum(cat.counts.p1$n[4:6]), 0)
pct.occ_core1 = round(100*cat.counts.p1$n[c(1,4)]/(cat.counts.p1$n[1] + cat.counts.p1$n[4]), 0)
pct.occ_tran1 = round(100*cat.counts.p1$n[c(3,6)]/(cat.counts.p1$n[3] + cat.counts.p1$n[6]), 0)

pct.bio_core.5 = round(100*cat.counts.p.5$n[1:3]/sum(cat.counts.p.5$n[1:3]), 0)
pct.bio_tran.5 = round(100*cat.counts.p.5$n[4:6]/sum(cat.counts.p.5$n[4:6]), 0)
pct.occ_core.5 = round(100*cat.counts.p.5$n[c(1,4)]/(cat.counts.p.5$n[1] + cat.counts.p.5$n[4]), 0)
pct.occ_tran.5 = round(100*cat.counts.p.5$n[c(3,6)]/(cat.counts.p.5$n[3] + cat.counts.p.5$n[6]), 0)

pdf('Results/Plots/EXP3/classification_errors_v_abundance.pdf', height = 8, width = 10)
par(mfcol = c(2,2), mar = c(4, 4, 2, 1), mgp = c(2.5,.5,0))

# 2 panels for p = 1
boxplot(log10(xc.sp.p1$Ngrid) ~ xc.sp.p1$xc, col = c(rep('cornflowerblue', 3), rep('salmon', 3)),
        ylab = "log10 Grid Abundance", xlab = "[biological classification]-[occupancy classification]",
        main = "Pixels >2/3 landscape similarity, P_obs = 1")
abline(h = c(4,4.5,5), lty = 'dotted')

boxplot(log10(xc.sp.p1$Ncell) ~ xc.sp.p1$xc, col = c(rep('cornflowerblue', 3), rep('salmon', 3)),
        ylab = "log10 Pixel Abundance", xlab = "[biological classification]-[occupancy classification]",
        ylim = c(0, 4))
abline(h = c(1,2), lty = 'dotted')
text(1:3, rep(3.3, 3), pct.bio_core1, col = 'cornflowerblue', cex = 1.5)
text(4:6, rep(3.3, 3), pct.bio_tran1, col = 'salmon', cex = 1.5)
text(c(1,4,3,6), rep(3.8, 4), c(pct.occ_core1, pct.occ_tran1), cex = 1.5)

# 2 panels for p = .5
boxplot(log10(xc.sp.p.5$Ngrid) ~ xc.sp.p.5$xc, col = c(rep('cornflowerblue', 3), rep('salmon', 3)),
        ylab = "log10 Grid Abundance", xlab = "[biological classification]-[occupancy classification]",
        main = "Pixels >2/3 landscape similarity, P_obs = .5")
abline(h = c(4,4.5,5), lty = 'dotted')

boxplot(log10(xc.sp.p.5$Ncell) ~ xc.sp.p.5$xc, col = c(rep('cornflowerblue', 3), rep('salmon', 3)),
        ylab = "log10 Pixel Abundance", xlab = "[biological classification]-[occupancy classification]",
        ylim = c(0, 4))
abline(h = c(1,2), lty = 'dotted')
text(1:3, rep(3.3, 3), pct.bio_core.5, col = 'cornflowerblue', cex = 1.5)
text(4:6, rep(3.3, 3), pct.bio_tran.5, col = 'salmon', cex = 1.5)
text(c(1,4,3,6), rep(3.8, 4), c(pct.occ_core.5, pct.occ_tran.5), cex = 1.5)
dev.off()

# Linear models for perfect detection
spmod.p1 = lm(occ ~ sim + Ncell, data = xc.sp.p1)
round(cor(xc.sp.p1[,c('occ', 'sim', 'Ncell', 'Ngrid')]), 2)

