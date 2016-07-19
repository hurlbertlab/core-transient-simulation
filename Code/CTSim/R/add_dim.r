#' Add array dimension
#'
#' Inserts a dimension into an array at a specified location
#'
#' The function inserts a dimension into the array \code{x}
#' at each dimension given in \code{where}.
#'
#' @param x (required) an array
#' @param where vector of natural numbers specifying where the new dimension(s) 
#' 	should occur (see details). Defaults to last dimension.
#' @return an array 
#'
#' @export

add_dim = function(x, where=NULL){

	if(is.null(where)) where = length(dim(x)) + 1
	
	# Get original max dimension
	n = ifelse(is.null(dim(x)), 1, length(dim(x))) 
	
	# Determine whether original has dimnames
	if(n>1) has_names = !is.null(dimnames(x))
	if(n==1) has_names = !is.null(names(x))
	
	# Make new dimensions and names
	new_dim = as.numeric(1:(length(where)+n) %in% where)
	new_names = lapply(1:length(new_dim), function(i) ifelse(i %in% where, 1, NA))
	old_names = which(new_dim==0)
	
	if(n==1){
		if(has_names) new_names[[old_names]] = names(x) 
		new_dim[new_dim==0] = length(x)
	} else {
		if(has_names) for(i in 1:length(old_names)) new_names[[old_names[i]]] = dimnames(x)[[i]] 
		new_dim[new_dim==0] = dim(x)
	}

	if(has_names){
		new_array = array(x, dim=new_dim, dimnames=new_names)
	} else {
		new_array = array(x, dim=new_dim)
	}

	# Return new array
	new_array	
}

