#' Classify occupancy vs residency
#'
#' Cross-classify species based on their temporal occupancy
#' and true status as core residents or transients
#'
#' This function classifies species as core or transient at a set of locations
#' based on temporal occupancy and also as as biologically core or transient
#' based on birth rates, \code{b_rates}, or a pre-defined 
#' \code{classification}. It then counts the number of species occuring 
#' in each of the four categories: biologically core - occupancy core, 
#' biologically core - occupancy transient, biologically transient - occupancy 
#' core, and biologically transient - occupancy transient. 
#' Occupancy-based core/transient categories are determined by breakpoints 
#' given in \code{breaks}. This can be either a named list with explicit
#' intervals or a vector of one or two numbers on the interval
#' \code{[0,1]}. In this case, transient species are those with occupancies 
#' less than or equal to \code{breaks[1]} while core species are those with 
#' occupancies greater than \code{max(breaks)}. The status of species as 
#' biologically core or transient is specified either by a matrix classifying  
#' each species as \code{'core'} or \code{'trans'} at each spatial location
#' (\code{classification}) or by supplying both \code{b_rates} and a vector of
#' the habitat types at each location (\code{habitats}). This should be
#' calculated for the same set of locations in \code{occupancy} using 
#' \code{\link{average_habitat}}.
#' 
#' The function can return either a count of the number of species in each category
#' or a list of the species that are in each category (specify with 
#' \code{return}). It can also return a separate table for each location if
#' \code{do_each=TRUE}.
#'
#' @param occupancy (required) matrix of species temporal occupancies at
#' 	a set of locations (as returned by \code{\link{calc_occupancy}})
#' @param breaks (required) either a numeric vector of length 1 or 2 giving 
#' 	occupancy breakpoints for transient vs. core status or a named list with 
#' 	specific intervals: \code{list(trans=c(), core=c())}
#' @param b_rates matrix of birth rates for each species in each habitat type.
#' 	Dimensions are \code{[species, habitat type]}.
#' @param habitats vector of habitat types (\code{'A'} or \code{'B'}) 
#' 	in each location 
#' @param classification matrix classifying each species at each location as 
#' 	\code{'core'} or \code{'trans'}. Dimensions are 
#' 	\code{[locations, species]}.
#' @param do_each logical indicating whether one table should be returned 
#' 	summarizing across all locations (\code{FALSE}) or an array with one table 
#' 	for each spatial unit \code{TRUE}. Defaults to \code{FALSE}.
#' @param return character string indicating whether a table of \code{'counts'}
#' 	or lists of species IDs (\code{'ids'}) should be returned. 
#' 	Default is \code{'counts'}.
#' @return table or array, depending on \code{do_each}

cross_classify = function(occupancy, breaks, b_rates=NULL, habitats=NULL, classification=NULL, do_each=F, return='counts'){

	# Convert breaks to list of intervals defining core and transient occupancy levels if just specified as numeric
	if(!is.list(breaks)){
		breaks = list(trans=c(0,breaks[1]), core=c(breaks[length(breaks)], 1))
	}

	# Remove species that were never observed at a site
	occupancy[occupancy==0] = NA

	# Classify species based on occupancy matrix
	
	if(breaks$trans[2]==breaks$core[1]){
		
		# Case when no occupancy levels excluded
		classes = matrix(cut(occupancy, breaks=c(breaks$trans, breaks$core[2]), include.lowest=T, labels=c('trans','core')), nrow=nrow(occupancy))
	
	} else {
		
		# Case when intermediate occupancy levels excluded
		classes = matrix(cut(occupancy, breaks=c(breaks$trans, breaks$core), include.lowest=T, labels=c('trans',NA,'core')), nrow=nrow(occupancy))
	}
	colnames(classes) = colnames(occupancy)

	# Get predefined classification or calculate from birth rates
	if(is.null(classification)){
		
		# Catch error where not enough information specified
		if(is.null(habitats)&is.null(b_rates)) stop('Must specify either a classification of each species in each spatial unit or habitat types and birth rates.')

		cores = t(sapply(habitats, function(h) b_rates[,h]>0))
		classification = apply(cores, 1:2, function(x) ifelse(x, 'core', 'trans'))		
	}

	# Check whether classification matrix matches occupancy matrix
	if(sum(dim(classification)==dim(classes))<2) stop('Incorrect dimensions for classification matrix. Should be sites x species.')

	# Calculate contingency tables
	if(do_each){
		# Calculate for each spatial unit
		if(return=='counts'){
			cross_tab = sapply(1:nrow(classes), function(i){
				table(classification=factor(classification[i,], levels=c('core','trans')), occupancy=factor(classes[i,], levels=c('core','trans')))
			}, simplify='array')
		}
		if(return=='ids'){
			cross_tab = sapply(1:nrow(classes), function(i){
				this_tab = array(list(), dim=c(2,2), dimnames=list(classification=c('core','trans'), occupancy=c('core','trans')))
				this_tab['core','core'] = list(as.character(which(classes[i,]=='core' & classification[i,]=='core')))
				this_tab['core','trans'] = list(as.character(which(classes[i,]=='trans' & classification[i,]=='core')))
				this_tab['trans','core'] = list(as.character(which(classes[i,]=='core' & classification[i,]=='trans')))
				this_tab['trans','trans'] = list(as.character(which(classes[i,]=='trans' & classification[i,]=='trans')))
				this_tab
			}, simplify='array')
		}
	} else {
		# Sumarize across all spatial units
		if(return=='counts') cross_tab = table(classification=factor(classification, levels=c('core','trans')), occupancy=factor(classes, levels=c('core','trans')))
		if(return=='ids'){
			cross_tab = array(list(), dim=c(2,2), dimnames=list(classification=c('core','trans'), occupancy=c('core','trans')))
			cross_tab['core','core'] = list(as.character(which(classification=='core'&classes=='core')))
			cross_tab['core','trans'] = list(as.character(which(classification=='core'&classes=='trans')))
			cross_tab['trans','core'] = list(as.character(which(classification=='trans'&classes=='core')))
			cross_tab['trans','trans'] = list(as.character(which(classification=='trans'&classes=='trans')))
		}
	}

	# Return tables
	cross_tab
}

