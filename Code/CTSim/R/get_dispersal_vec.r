#' Calculate dispersal vector
#'
#' Stochastically generates dispersal vectors given a dispersal kernel 
#' and expected dispersal distance. 
#'
#' Only implements isotropic dispersal (no directional bias). Three types 
#' of dispersal kernels are currently implemented and should be specified 
#' with the \code{type} argument in the list which is passed as the parameter 
#' \code{form}. These are:
#' \describe{
#' 	\item{\code{'gaussian'}}{Gaussian kernel- dispersal distance sampled from
#' 		a half-normal distribution with mean \code{d}.}
#' 	\item{\code{'adjacent'}}{Adjacent cell dispersal (rook's move)- \code{d} 
#' 		is the probability that an individual will move. An
#' 		additional argument \code{moves} must be passed in \code{form} 
#' 		giving the number of steps an individual should attempt to take.}
#' 	\item{\code{'uniform'}}{Uniform dispersal- individual has equal 
#' 		probability of occurring anywhere within a circle with radius 
#' 		\code{d}. Sufficiently large \code{d} allows for unlimited dispersal.}
#' }
#'
#' @param d (required) numeric controlling dispersal distance (see details). 
#' @param form list defining the dispersal kernel (see details).
#' 	Minimally contains an argument \code{type} defining the form of the 
#' 	dispersal kernel. Defaults to Gaussian.
#' @param N number of vectors to generate
#' @return a matrix of dispersal vectors. First column gives the magnitude 
#' and the second column gives the direction in radians, counter-clockwise
#' from East.  
#'
#' @seealso \code{\link{disperse}} for dispersal of propagules with different 
#' expected dispersal distances.
#'
#' @import fdrtool

get_dispersal_vec = function(d, form=list(type='gaussian'), N=1){
	
	# Gaussian kernel
	if(form$type=='gaussian'){
	
		# Calculate theta for mean = d
		theta = 1/d

		# Generate distance
		r = fdrtool::rhalfnorm(N, theta)	

		# Generate random angle (in radians)
		phi = runif(N, 0, 2*pi)

	}

	# Adjacent cell dispersal
	# d is the probability that a propagule will leave its origin cell
	# form$moves = maximum number of moves that a propagule will complete in one timestep
	if(form$type=='adjacent'){

		# Generate random set of directional moves
		moves = sapply(1:N, function(n){
			sample(0:1, size=form$moves, prob= c(1-d, d), replace=T) * sample(1:4, size=form$moves, replace=T)
		})

		# Add a dimension if necessary (e.g. form$moves = 1)
		if(!is.matrix(moves)) moves = add_dim(moves, 1)


		# Calculate distance and angle
		dX = apply(moves, 2, function(x) -1*sum(x==3) + 1*sum(x==1))
		dY = apply(moves, 2, function(x) -1*sum(x==4) + 1*sum(x==2))
		r = sqrt(dX^2 + dY^2)
		phi = atan2(dY,dX)
	}

	# Uniform: implements no disersal limitation
	# d is the maximum dispersal distance
	if(form$type=='uniform'){
	
		# Generate distance
		r = runif(N, 0, d)

		# Generate random angle		
		phi = runif(N, 0, 2*pi)	
	}


	# Return vector
	cbind(r, phi)
}
