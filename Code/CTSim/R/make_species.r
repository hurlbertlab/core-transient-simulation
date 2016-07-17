#' Creation of species pool
#'
#' Creates an array to hold species vital rates.
#'
#' A species pool is defined by a 3-dimensional array
#' in which the first dimension holds species, the second
#' defines which habitat the rate applies to ('A' or 'B'), and the third 
#' hols the type of rate: 
#' \enumerate{
#' 	\item \code{'b'} birth = number of offpsring produced by established
#'		 individual per timestep
#'	\item \code{'m'} mortality = probability that established individual
#'		 dies during a time step
#'	\item \code{'r'} recruitment = probability than a propagule arriving
#'		 at an open space will become an established individual
#'	\item \code{'d'} dispersal = expected distance that a propagule will 
#'		disperse from its origination point 
#' 	\item \code{'v'} movement = expected distance that an established 
#' 		individual will move from its cell
#'}
#' Birth rates are positive in a species' preferred habitat and 0 elsewhere.
#' Generalists prefer both habitat types equally. Dispersal rates control how
#' newly produced propagules move away from their cell of origin. Movement rates
#' control how established individuals move from their current cell. Movement,
#' mortality, and recruitment rates can be set to differe systematically between
#' preferred and non-preferred habitats. 
#'
#' @param S_A (required) number of specialist species for habitat type A
#' @param S_B (required) number of specialist species for habitat type B
#' @param S_AB number of generalist species (defaults to 0)
#' @param dist_b named list describing the distribution from 
#' 	which birth rates should be drawn (see \code{distribution} 
#' 	parameter in \code{\link{make_sad}}). Defaults to uniform on [1,10].
#' @param m (required) vector of length 2 or 3 specifying death rates 
#' 	[preferred habitat, non-preferred habitat, generalist rate]
#' @param r (required) vector of length 2 or 3 specifying recruitment rates 
#' 	[preferred habitat, non-referred habitat, generalist rate]
#' @param dist_d named list of parameters describing distribution from which 
#' 	propagule dispersal kernals should be drawn. Defaults to 1 for all species.
#' 	Contains:
#'	\describe{
#' 		\item{\code{u}}{mean dispersal distance}
#' 		\item{\code{var}}{variance of distribution from which kernals drawn. 
#' 			Defaults to 0 for the same dispersal kernal for all species.}
#'	}
#' @param dist_v named list of parameters describing distribution from which 
#' 	movement kernals should be drawn. Defaults to 0 for all species in all 
#' 	habitats. Contains:
#' 	\describe{
#'		\item{\code{mu}}{vector of length 2 or 3 specifying mean movement 
#' 			distance [preferred habitat, non-preferred habitat, generalist rate]}
#'		\item{\code{var}}{vector of length 2 or 3 specifying variance in mean 
#' 			movement distance [preferred habitat, non-preferred habitat, 
#' 			generalist rate]. Defaults to 0 for the same movement kernal for 
#' 			all species.}
#'	}
#' @return an array with dimensions \code{[S_A+S_B+S_AB, 2, 5]}
#'

make_species = function(S_A=NA, S_B=NA, S_AB=NA, dist_b = NULL, m, r, dist_d=NULL, dist_v=NULL){

	# Catch error if no species specified
	if(is.na(S_A)|is.na(S_B)) stop('Must specify number of species.')

	# Set generalists to 0 if undefined
	if(is.na(S_AB)) S_AB = 0

	# Set distributions to defaults if undefined
	if(is.null(dist_b)) dist_b = list(maxN=10, type='uniform')
	if(is.null(dist_d)) dist_d = list(mu=1, var=0)
	if(is.null(dist_v)) dist_v = list(mu=rep(0,3), var=rep(0,3))

	# Create array to hold species vital rates (b = birth, m = death, r = recruitment, d = dispersal, v = movement)
	# Rates are per timestep
	N_S = S_A + S_B + S_AB
	species_rates = array(NA, dim=c(N_S, 2, 5), dimnames=list(species=1:N_S, habitat=c('A','B'), rate=c('b','m','r','d','v')))

	# Specialist species birth rates in their preferred habitat ranged from 1 to maxN and are 0 in the unpreferred habitat
	species_rates[1:S_A, 'A', 'b'] = make_sad(S_A, dist_b)
	species_rates[1:S_A, 'B', 'b'] = 0
	species_rates[(S_A+1):(S_A+S_B), 'B', 'b'] = make_sad(S_B, dist_b)
	species_rates[(S_A+1):(S_A+S_B), 'A', 'b'] = 0

	# Death rates are prespecified constants
	# This could be modified to be drawn from a distribution
	species_rates[1:S_A, 'A', 'm'] = m[1]
	species_rates[1:S_A, 'B', 'm'] = m[2]
	species_rates[(S_A+1):(S_A+S_B), 'B', 'm'] = m[1]
	species_rates[(S_A+1):(S_A+S_B), 'A', 'm'] = m[2]

	# Specialist species recruitment rates are prespecified constants
	# This could be modified to be drawn from a distribution
	species_rates[1:S_A, 'A', 'r'] = r[1]
	species_rates[1:S_A, 'B', 'r'] = r[2]
	species_rates[(S_A+1):(S_A+S_B), 'B', 'r'] = r[1]
	species_rates[(S_A+1):(S_A+S_B), 'A', 'r'] = r[2]

	# Mean dispersal distances for each species are drawn from a gamma distribution with mean 'mu' and variance 'var'
	if(dist_d$var==0){
		species_rates[,,'d'] = dist_d$mu
	} else {
		theta = dist_d$var / dist_d$mu
		k = dist_d$mu / theta
		species_rates[,,'d'] = rgamma(N_S, shape=k, scale=theta)
	}

	# Mean movement distances for each species are drawn from a gamma distribution with mean 'mu' and variance 'var'.
	# Preferred habitat
	if(dist_v$var[1]==0){
		species_rates[1:S_A, 'A', 'v'] = dist_v$mu[1]
		species_rates[(S_A+1):(S_A+S_B), 'B', 'v'] = dist_v$mu[1]
	} else {
		theta = dist_v$var[1] / dist_v$mu[1]
		k = dist_v$mu[1] / theta
		species_rates[1:S_A, 'A', 'v'] = rgamma(S_A, shape=k, scale=theta)
		species_rates[(S_A+1):(S_A+S_B), 'B', 'v'] = rgamma(S_B, shape=k, scale=theta)
	}
	# Unpreferred habitat
	if(dist_v$var[2]==0){
		species_rates[1:S_A, 'B', 'v'] = dist_v$mu[2]
		species_rates[(S_A+1):(S_A+S_B), 'A', 'v'] = dist_v$mu[2]
	} else {
		theta = dist_v$var[2] / dist_v$mu[2]
		k = dist_v$mu[2] / theta
		species_rates[1:S_A, 'B', 'v'] = rgamma(S_A, shape=k, scale=theta)
		species_rates[(S_A+1):(S_A+S_B), 'A', 'v'] = rgamma(S_B, shape=k, scale=theta)
	}
	

	# Assign rates for generalists, if they exist
	if(S_AB > 0){
	
		# Generalist birth rates range from 1 to maxN in both habitats, but is not necessarily the same in both habitats
		species_rates[(N_S-S_AB+1):N_S, ,'b'] = make_sad(2*S_AB, dist_b)
	
		# Generalist death rates default to the death rate for the preferred habitat, but can be specified separately
		species_rates[(N_S-S_AB+1):N_S, ,'m'] = ifelse(length(m)==3, m[3], m[1])

		# Generalist recruitment rates default to the recruitment rate for the preferred habitat, but can be specified separately
		species_rates[(N_S-S_AB+1):N_S, ,'r'] = ifelse(length(r)==3, r[3], r[1])

		# Generalist movement rates are the same in both habitat types and default to the movement rate in the preferred habitat, but can be specificed separately
		if(length(dist_v$mu)==2) dist_v$mu = c(dist_v$mu, dist_v$mu[1])
		if(length(dist_v$var)==2) dist_v$var = c(dist_v$var, dist_v$var[1])

		if(dist_v$var[3]==0){
			species_rates[(N_S-S_AB+1):N_S, , 'v'] = dist_v$mu[3]
		} else {
			theta = dist_v$var[3] / dist_v$mu[3]
			k = dist_v$mu[3] / theta
			species_rates[(N_S-S_AB+1):N_S, ,'v'] = rgamma(S_AB, shape=k, scale=theta)
		}
	}

	# Return rates
	species_rates
}
