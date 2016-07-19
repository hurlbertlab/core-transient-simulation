#' Aggregate cells
#'
#' Groups cells of a grid together and returns a list
#' of cell locations.
#'
#' This function is used for summarizing simulation output across different 
#' spatial grains. The function returns a list of matrices, where each matrix
#' specifies the locations of cells to be aggregated. This list can be passed
#' as the \code{locs} parameter in summary functions such as 
#' \code{\link{summarize_sim_N}} and the functions therein.
#' The function implements three methods for aggregating cells:
#' 	\describe{
#'		\item{'partition'}{The grid is partitioned so that groups with 
#'			dimensions \code{(dX,dY)} are non-overlapping, starting from the  
#'			first row, first column. If the number of rows and columns in 
#'			the grid is not divisible by \code{dX} or \code{dY} then the 
#'			last aggregation in each row or column will contain fewer columns 
#' 			or rows of cells than specified.}
#'		\item{'window'}{A window with dimensions \code{(dX, dY)} is slid
#'			repeatedly across the grid according the the number of steps given  
#'			in \code{slide}. If \code{slide} is less than \code{(dX, dY)}
#'			this results in overlapping groups of cells.}
#'		\item{'origin'}{Groups of cells with dimensions \code{(dX,dY)}   
#'			are aggregatedaround a focal cell or set of 
#'			cells given in \code{locs}.}
#'	}
#'
#' @param X (required) x-dimension of the grid or a vector of length 2 giving 
#' 	the dimensions of the grid \code{(x, y)}
#' @param Y y-dimension of the grid
#' @param dX (required) number of rows to aggregate or a vector of length 2 
#' 	giving the dimensions to aggregate
#' @param dY number of columns to aggregate
#' @param form character string indicating how cells should be aggregated 
#' 	(see details).
#' @param slide if \code{form = 'window'}, a vector of length 1 or 2 indicating
#' 	the number of cells to slide the window in the x and y directions.
#' 	Default is one cell in each direction.
#' @param locs if \code{form = 'origin'}, a matrix of cell locations on 
#' 	which aggregations should be centered. Default is the center cell.
#' @return a list of matrices

aggregate_cells = function(X, Y=NA, dX, dY=NA, form=NA, slide=NULL, locs=NULL){
	
	# Catch errors where not enough parameter specified
	if(length(X)==1 & is.na(Y)) stop('Must specify dimensions of landscape.')
	if(length(dX)==1 & is.na(dY)) stop('Must specify aggregation dimensions dX x dY.')

	# Define missing parameters to defaults.
	if(length(X)==2){
		Y = X[2]
		X = X[1]
	}

	if(length(dX)==2){
		dY = dX[2]
		dX = dY[1]
	}

	if(is.na(form)) form = 'partition'
	
	if(form=='window'){
		if(is.null(slide)) slide = c(1,1)
		if(length(slide)==1) slide = rep(slide, 2)
	}

	if(form=='origin' & is.null(locs)){
		locs = matrix(c(ceiling(X/2), ceiling(Y/2)), nrow=1)		
	}
	
	# Aggregate cells
	if(form=='partition'){

		# Determine starting points
		Xs = seq(1, X, dX)
		Ys = seq(1, Y, dY)
		group_origins = as.matrix(expand.grid(Xs, Ys))

		# For each origin cell, find all cells within dX, dY
		groups = sapply(1:nrow(group_origins), function(i){
			o = group_origins[i,]
			Xcoords = o[1]:min(o[1]+dX-1, X)
			Ycoords = o[2]:min(o[2]+dY-1, Y)
			cells = as.matrix(expand.grid(Xcoords, Ycoords))
			colnames(cells) = c('x','y')
			list(cells)		
		})

	}
	if(form=='window'){
	
		# Determine starting points
		Xs = seq(1, X, slide[1])
		Ys = seq(1, Y, slide[2])
		group_origins = as.matrix(expand.grid(Xs, Ys))

		# For each origin cell, find all cells within dX, dY
		groups = sapply(1:nrow(group_origins), function(i){
			o = group_origins[i,]
			Xcoords = o[1]:min(o[1]+dX-1, X)
			Ycoords = o[2]:min(o[2]+dY-1, Y)
			cells = as.matrix(expand.grid(Xcoords, Ycoords))
			colnames(cells) = c('x','y')
			list(cells)		
		})
	}
	if(form=='origin'){

		# For each given origin cell, find all cells within dX, dY centered on that cell
		groups = sapply(1:nrow(locs), function(i){
			o = locs[i,]

			if(dX%%2==1){
				Xcoords =(o[1]-floor(dX/2)):(o[1] + floor(dX/2))
			} else {
				Xcoords = (o[1]-dX/2+1):(o[1]+dX/2)
			}

			if(dY%%2==1){
				Ycoords =(o[2]-floor(dY/2)):(o[2] + floor(dY/2))
			} else {
				Ycoords = (o[2]-dY/2+1):(o[2]+dY/2)
			}
			cells = as.matrix(expand.grid(Xcoords, Ycoords))
			colnames(cells) = c('x','y')
			list(cells)	
		})
	}

	# Return list of locations
	groups
}
