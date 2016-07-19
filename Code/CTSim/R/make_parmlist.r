#' Make list of simulation parameters
#'
#' Creates a list of parameters to used in a simulation from objects present
#' in a given environment.
#'
#' The function searches the current environment (or one given in 
#' \code{e}) for a set of required and optional parameters used by
#' \code{\link{run_sim}} to run a simulation. Parameters are saved in 
#' a named list which can be used in the \code{parms} argument of
#' \code{\link{run_sim_N}} or written to a file using \code{\link{write_parms}}
#' (not yet implemented).
#' The following parameters are \strong{required} for simulation and \emph{must be present
#' in the environment or the function will fail}:
#' \itemize{
#'	\item \strong{\code{dimX}} : x-dimension of landscape
#'	\item \strong{\code{dimY}} : y-dimension of landscape
#'	\item \strong{\code{S_A}} : number of specialist species on habitat type A
#'	\item \strong{\code{S_B}} : number of specialist species on habitat type B
#'	\item \strong{\code{m_rates}} : vector of mortality rates on preferred and 
#'		non-preferred habitats
#'	\item \strong{\code{r_rates}} : vector of recruitment rates on preferred 
#'		and non-preferred habitats
#'	\item \strong{\code{K}} : carrying capacity of cells. 
#'	\item \strong{\code{nsteps}} : number of timesteps to simulate
#'	\item \strong{\code{nruns}} : number of independent simulations to run
#' }
#' For further details on \code{S_A}, \code{S_B}, \code{m_rates} and 
#' \code{r_rates} see \code{\link{make_species}}. For further details on 
#' \code{K} see \code{\link{populate_landscape}}.
#'
#' The following parameters are optional and more information can 
#' be found in the documentation on the functions they are passed to:
#' \itemize{
#'	\item Parameters passed to \code{\link{make_landscape}}
#'		\itemize{
#'			\item \strong{\code{vgm_dcorr}} : distance at which habitat values become
#'				uncorrelated
#'			\item \strong{\code{vgm_mod}} : variogram model controling spatial 
#'				autocorrelation of habitat values
#'			\item \strong{\code{habA_prop}} : proportion of landscape that comprised of
#'				habitat type A
#'		}
#'	\item Parameters passed to \code{\link{make_species}}
#'		\itemize{
#'			\item \strong{\code{S_AB}} : number of generalist species
#'			\item \strong{\code{dist_b}} : list defining the distribution from which 
#'				species' birth rates are sampled
#'			\item \strong{\code{dist_d}} : list defining the distribution from which
#'				species' dispersal rates are sampled. Must contain character string
#'				named \code{type}.
#'			\item \strong{\code{dist_v}} : list defining the distribution from which
#'				species' movement rates are sampled. Must contain character string
#'				named \code{type}.
#'			\item \strong{\code{dist_gsad}} : list defining distribution from which
#'				global species abundances are sampled or 'b_rates', indicating 
#'				that the gsad should match species birth rates in their preferred 
#'				habitat.
#'		}
#'	\item Parameters passed to \code{\link{populate_landscape}}
#'		\itemize{
#'			\item \strong{\code{prop_full}} : proportion of the landscape's carrying
#'				capacity that should initially contain individuals
#'			\item \strong{\code{init_distribute}} : character string indicating how
#'				individuals should be initially distributed across the landscape
#'			\item \strong{\code{cells_distribute}} : if \code{init_distribute} is 
#'				is 'designated', a matrix giving the locations of cells in which
#'				to place propagule
#'		}
#'	\item Parameters passed to \code{\link{run_sim}}
#'		\itemize{
#'			\item \strong{\code{d_kernel}} : list defining the shape of the dispersal
#'				kernel of new propagules.
#'			\item \strong{\code{v_kernel}} : list defining the shape of the movement
#'				kernel of established individuals
#'			\item \strong{\code{imm_rate}} : immigration rate- probability than an 
#'				empty space will be colonized by a migrant from outside
#'				the metacommunity
#'		}
#'	\item Parameters passed to \code{\link{run_sim_N}}
#'		\itemize{
#'			\item \strong{\code{save_steps}} : vector of timesteps to save in each 
#'				simulation
#'			\item \strong{\code{simID}} : character string that identifies simulations
#'				run on this set of parameters
#'		}
#'}
#'
#' @param e environment in which to look for parameters
#' @return a named list
#'
#' @note Developers should note that this function must be manually  
#' updated whenever new parameters are added to simulation functions.

make_parmlist = function(e=parent.frame()){

	# Define list of required parameters
	parms = list(
		dimX = e$dimX,
		dimY = e$dimY,
		S_A = e$S_A, 
		S_B = e$S_B,
		m_rates = e$m_rates,
		r_rates = e$r_rates,
		K = e$K,								
		nsteps = e$nsteps,
		nruns = e$nruns
	)

	# Add on optional parameters
	if(exists('vgm_dcorr', e)) parms = c(parms, vgm_dcorr = e$vgm_dcorr)			
	if(exists('vgm_mod', e)) parms = c(parms, vgm_mod = list(e$vgm_mod)) 
	if(exists('habA_prop', e)) parms = c(parms, habA_prop = e$habA_prop)
	if(exists('S_AB', e)) parms = c(parms, S_AB = e$S_AB)
	if(exists('dist_b', e)) parms = c(parms, dist_b = list(e$dist_b))
	if(exists('dist_d', e)) parms = c(parms, dist_d = list(e$dist_d))
	if(exists('dist_v', e)) parms = c(parms, dist_v=list(e$dist_v))
	if(exists('dist_gsad',e)){
		if(is.list(e$dist_gsad)){
			parms = c(parms, dist_gsad = list(e$dist_gsad))
		} else {
			parms = c(parms, dist_gsad = e$dist_gsad)
		}
	}
	if(exists('prop_full', e)) parms = c(parms, prop_full = e$prop_full)
	if(exists('init_distribute', e)) parms = c(parms, init_distribute = e$init_distribute)
	if(exists('cells_distribute', e)) parms = c(parms, cells_distribute = list(e$cells_distribute))
	if(exists('d_kernel', e)) parms = c(parms, d_kernel = list(e$d_kernel))
	if(exists('v_kernel', e)) parms = c(parms, v_kernel = list(e$v_kernel))
	if(exists('imm_rate', e)) parms = c(parms, imm_rate = e$imm_rate)
	if(exists('save_steps', e)) parms = c(parms, save_steps = list(e$save_steps))
	if(exists('simID', e)) parms = c(parms, simID = e$simID)
}



