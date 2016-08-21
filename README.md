# core-transient-simulation

This project is a spatially explicit simulation of metacommunity dynamics that models habitat filtering and dispersal.

## Simulation Description
The simulation occurs on a regular grid of cells (termed the 'landscape') and each cell has discrete habitat type, 'A' or 'B'. Each cell contains a 'community' of individuals, which is fixed at a pre-determined carrying capacity. Communities may contain fewer individuals than the carrying capacity, but not more. The species pool consists of three types of species- generalists, and habitat specialists on habitat type A or B. Generalists produce offspring in both habitat types while specialists only produce offspring in their preferred habitat. After initializing a simulation landscape and species pool, the simulation proceeds through repeated iteration of four processes: birth - dispersal - establishment - death. Once established, individuals vacate their place in the community only through dispersal or death- the simulation does not model competitive displacement. 

## Directory Structure
**/Code**
  + **./CTSim**: package source code and manual files
  + **./HTML**: html files for static documentation
  + **./Parameters**: parameter files for running and summarizing simulations. Subdirectories contain parameter files for different experiments described [here](/wiki/Experiments). Parameter files for running simulations begin with 'p_' whereas parameter files for summarizing simulations begin with 's_'.
  + **./Scripts**: miscellaneous scripts for working with CTSim
   + [debugging.R](/Code/Scripts/debugging.R): code for testing CTSim package
   + [make_package.R](/Code/Scripts/make_package.R): code for building CTSim package
   + [make_parmfiles.R](/Code/Scripts/make_parmfiles.R): code for generating parameter files used in experiments
   + [submit_runs.txt](/Code/Scripts/submit_runs.txt): shell code for submitting multiple jobs to cluster
   + [visualize_simulation.R](Code/Scripts/visualize_simulation.R): code for analyzing  CTSim runs and experiments
  + Current package version for [Windows](/Code/CTSim_0.1.3.tar.gz) and [Unix](/Code/CTSim_0.1.3.zip)

**/Results**
 + **./Plots**: visualizations of simulation results, organized by experiment/run
 + **./Summary**: summary RData files of simulation results, organized by experiment/run

## HTML Help Files
`/Code/HTML`

Static html help files for every function in the CTSim package are available [here](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/index.html).

## CTSim Package 
`/Code/CTSim/`

Files for installing the package are available in the [Code](/Code/) directory. The most recent package version is the one with the highest version number. 

### General Workflow

1. Create a set of simulation parameters and save in a directory (e.g. 'Parms'). See [example_parameter_file.txt](/Code/Parameters/example_parameter_file.txt) for an example file and  [make_parmfiles.R](/Code/Scripts/make_parmfiles.R) for an example of generating multiple parameter files. All parameter file names must begin with 'p_'.
2. Decide on a directory where simulation results should be saved (e.g. 'Results').
3. Move the run_simulation.R script from the 'exec/' directory of the CTSim package to your working directory.
4. Run R in batch mode on N cores with the run_simulation.R script: 

`R CMD BATCH "--args N Parms/ Results/" run_simulation.R myrun.Rout`

5. If the simulation stops before all runs are complete, restart the simulation without writing over already finished results using the 'restart_simulation.R' script (also found in the 'exec/' directory):

`R CMD BATCH "--args N Parms/ Results/" restart_simulation.R myrun_restart.Rout`

6. Create a summary parameter file (or multiple files) and save in a directory (e.g. 'Parms'). See [example_summary_parameter_file.txt](/Code/Parameters/example_summary_parameter_file.txt) for an example file. All summary parameter file names must begin with 's_'.
7. Decide on a directory where simulation result summarise should be saved (e.g. 'Summaries').
8. Move the summarize_simulation.R script from the 'exec/' directory of the CTSim package to your working directory.
9. Run R in batch mode with the summarize_simulation.R script:

`R CMD BATCH "--args Results/ Parms/ Summaries/" summarize_simulation.R myrun_summary.Rout`

10. Examine simulation summaries using your own code, or see [visualize_simulation.R](/Code/Scripts/visualize_simulation.R) for examples of how to load and examine summaries from simulations on different sets of parameters.

### Simulation Intialization
A simulation requires three objects: a landscape, a species pool, and a global species abundance distribution (hereafter, gsad). The functions used to initialize these objects are described below.

#### Initializing a landscape:
[`make_landscape(x, y, mod, d, prop, draw_plot)`](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/make_landscape.html)

A landscape is a raster layer whose values are -1 or 1, corresponding to habitat types 'A' and 'B', respectively. 
When a landscape is initialized, the user can specify its size (`x,y`), the proportion of the cells which should belong to habitat type 'B' (`prop`), a variogram model used to define the spatial autocorrelation of habitat values (`mod`), and the distance at which habitat values become uncorrelated (e.g. the range of the variogram model: `d`). Variogram models are implemented by the `vgm()` function in gstat. Use `show.vgms()` to see available models. Only the dimensions of the landscape are required. By default, the function will return a grid with 50% of habitat type 'A' and a exponential variogram model with partial sill = 1 and range = 1/3 of the grids smallest dimension.

#### Initializing a species pool
[`make_species(S_A, S_B, S_AB, dist_b, m, r, dist_d, dist_v)`](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/make_species.html)

A species pool is a 3-dimensional array of species vital rates. The first dimension specifies the species. The second dimension defines in which habitat the rate applies ('A' or 'B') and the third dimension specifies the type of rate:

1. **b**: birth rate- number of offpsring produced by established individual per timestep.
2. **m**: mortality rate- probability that established individual dies during a time step.
3. **r**: recruitment rate- probability than a propagule arriving at an open space will become an established individual
4. **d**: dispersal rate- expected distance that a propagule will disperse from its origination point
5. **v**: movement rate- expected distance that an established individual will move from its cell

Birth rates are positive in a species' preferred habitat and 0 elsewhere. Generalists prefer both habitat types equally. Dispersal rates control how newly produced propagules move away from their cell of origin. Movement rates control how established individuals move from their current cell. Movement, mortality, and recruitment rates can be set to differ systematically between preferred and non-preferred habitats.

#### Initializing a GSAD
[`make_sad(N_S, distribution)`](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/make_sad.html)

A GSAD is a vector of species relative abundances on the "mainland". These define the relative probabilities that an immigrant from outside the landscape will belong to each species. Implemented distributions are 'same', 'uniform', 'power', 'logseries', 'lognormal', and 'poisson'. See the help file on [`make_sad`](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/make_sad.html) for further details.

### Simulation Operation
[`populate_landscape(land, species, gsad=NULL, K, distribution=NA, p=NA, which_cells=NULL)`](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/populate_landscape)

[`run_sim(steps, metacomm, land, species, gsad, d_kernel=NULL, v_kernel=NULL, imm_rate=NA, save_steps=NULL, ...)`](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/run_sim.html)

Once a landscape, species pool and GSAD have been generated, a simulation can be run. First, the landscape must be populated by species using the function `populate_landscape`. This function has several options for defining the number of individuals in each cell (`K`), the proportion of this carrying capacity which should be initially filled (`p`), as well as how individuals should be  distributed across the landscape (`distribution`). The `which_cells` option allows the user to specify the exact cells where individuals should be placed and the number to place in each cell. Individuals are drawn probabilistically from the specified GSAD (`gsad`), which defaults to the same probability for all species.

The populated landscape is a metacommunity object which takes the form of a matrix of lists. The `x` and `y` positions define the location of each community on the landscape, while the list defines which individuals are present. The length of these lists are fixed for the remainder of the simulaton to the carrying capacity defined in `K`. Integers refer to species identities and a `0` indicates that no individual is present. For example, if `metacomm[1,2]` is equal to `list(3,3,2,10,0,2,1)`, then the cell in the first row and second column can hold seven individuals and currently contains one individual of species 1, two individuals of species 2 and 3, one individual of species 10, and one empty space. 

After a metacommunity has been generated, a simulation can be run for a fixed length of time (`steps`) using the `run_sim` function. In addition to the objects defining the metacommunity (`metacomm`), landscape (`land`), species pool (`species`), and GSAD (`gsad`), the function also requires the probability that an empty space will be colonized by an individual from the "mainland" (`imm_rate`) as well as information about how propagules disperse (`d_kernel`) and individuals move after they have established in a cell (`v_kernel`). Three dispersal modes are implemented: half-gaussian, adjacent cell, and uniform. In all cases, the direction of dispersal is random (isotropic) and the expected dispersal distances are determined by dispersal and movement rates for individual species defined in the species pool object. See the help file on [`get_dispersal_vec`](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/get_dispersal_vec.html) for additional information on specifying dispersal and movement parameters. 

A simulation runs by iteratively calling the function [`run_timestep`](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/run_timestep.html), which defines the operations that occur in a single timestep. This function calls the main process functions [`die`](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/die.html), [`reproduce`](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/reproduce.html), [`disperse`](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/disperse.html), [`establish`](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/establish.html) to progress a metacommunity through one simulation timestep. Each of these functions are described in detail in their help files. The order of operations is as follows:

1. **Death**: Established individuals in each community experience	probabilistic mortality according to species- and habitat-specific mortality rates provided in `species`.
2. **Birth**: Established individuals in each community produce propagules according to species- and habitat-specific birth rates provided in `species`.
3. **Dispersal**: New propagules from each community disperse across the landscape (`land`) according to species-specific dispersal rates provided in `species`. The parameter `d_kernel` specifies the dispersal kernel for new propagules.
4. **Movement**: Established individuals in each community move from their current cell with species- and habitat-specific movement rates provided in `species`. The parameter `v_kernel` specifies the dispersal kernel for previously established individuals.
5. **Establishment**: Empty spaces in each community are colonized by either a migrant from outside the community with probability `imm_rate` or by an individual selected at random from the pool of new propagules and moving individuals that arrived in the cell. External migrants are chosen probabilistically from the relative abundances given in `gsad`. 

Once the simulation has run for a fixed number of timesteps, results are returns as an array of lists with three dimensions, where the third dimension indicates the timestep and the other two form a metacommunity object. An optional paramter `save_steps` can be used to only record data from specific timesteps, but by default all time steps are returned, including the initial metacommunity. Dimension names record the actual timestep (`0` for the initial metacommunity, `1` for the state of the community after one step, etc...).

### Multiple Simulation Runs

Most users will primarily want to run simulations using the [`run_sim_N`](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/run_sim_N.html) and [`run_sim_P`](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/run_sim_P.html) functions. These functions run multiple simulations on a set of parameters or multiple sets of parameters, respectively. They are useful because the user does not need to initialize simulation objects (landscape, species pool, etc...), as this occurs internally. Results can be returned or saved to a directory.

#### Multiple simulation runs on one set of parameters
[`run_sim_N(nruns, parms, nparallel=1, simID='test', save_sim=NULL, report=0, return_results=T, restart=F, lib_loc=NULL)`](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/run_sim_N.html)

This function runs multiple independent simulations on a set of parameters given as a list (`parms`). A parameter list can be generated by the function [`make_parmlist`](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/make_parmlist.hrml) which compiles simulation parameters stored in the current (or specified) environment into a list. See the section on Parameter Files below for a list of parameter objects which are required versus optional.

Each run generates a new landscape, species pool, global species abundance distribution and initial metacommunity. Landscapes, species pools and gsads are saved as lists (object names: `lands_N`, `species_N`, `gsad_N`) in a file named `<simID>_simobjects.RData`. Simulations can be run in parallel by specifying `nparallel > 1`, which requires the `doParallel` and `foreach` packages. By default, `nparallel = 1` and the simulations proceed serially. Each run of the simulation is temporarily saved to the working directory or permanently saved to the directory specified by `save_sim`. If this directory does not exist then it is created. Runs are saved in a subdirectory named `simID` as `<simID>_run<i>.RData`. This RData file contains five objects:

1. **results**: an array of the metacommunity through time returned by `run_sim`
2. **this_land**: the landscape used for the simulation
3. **this_species**: the species pool used in the simulation
4. **this_gsad**: the global species abundance distribution used in the simulation
5. **this_metacomm**: the initial metacommunity
 
If `return_results`is `TRUE`, then after all simulations are complete, all runs are read back into memory, compiled arrays, and returned as a list by the function. The componenets of this list are `results`, `species`, `lands`, and `gsads` and the first dimension of each array is the run to which it corresponds. For example, `results[1,,,]` is the metacommunity from run 1 through time. If `return_results` is `FALSE` then temporary files from the simulation runs are NOT removed, even if `save_sim` is not specified.  Users should be careful returning results for large simulations since the results object is generally too large for the working memory. It is better to save the results and then use summary functions (described below) to access them.

`run_sim_N` can be used to restart a set of multiple simulation runs, but does not currently allow users to restart a simulation on an existing run. If `restart` is `TRUE`, the function navigates to the `save_sim` directory and searches for the first instance of `i <= nruns` where `<simID>_run<i>.RData` does not exist. It then starts simulations for all `i` that do not have saved files, using the objects saved in `<simID>_simobjects.RData`.

By default, `run_sim_N` runs in silent mode, but by specifying a number for `report` users can request a timestamp to be written to STDOUT every time a fixed number of timesteps have passed. This can be useful for gauging how long simulations are taking. Finally, if `CTSim` is installed in a directory not on the default search path, users should indicate where the package is installed using `lib_loc`.

#### Multiple simulation runs on multiple sets of parameters
[`run_sim_P(ncores=1, parm_dir='./', results_dir='./Results/', sim_dir=NULL, report=0, restart=F)`](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/run_sim_P.html)

This function is a wrapper for `run_sim_N` which sequentially calls the function on each parameter file in the directory `parm_dir`, defaulting to the working directory. Simulation runs are saved to a subdirectory of `results_dir` named `simID`, which must be defined in each parameter file. Using the same `simID` in multiple parameter files will cause results to be saved over one another. Users can generate parameter files by creating a parameter list using `make_parmlist` and then saving this object to a file using `write_parms`. For an example see [make_parmfiles.R](./Code/Scripts/make_parmfiles.R). `sim_dir` refers to the directory where `CTSim` is installed, if not on the default search path. Specifying `ncores > 1` will run the simulations in parallel requesting the defined number of cores. 

#### Running simulations in batch mode
`R CMD BATCH "--args ncores parm_dir results_dir sim_dir report" run_simulation.R outfile.Rout`

`R CMD BATCH "--args ncores parm_dir results_dir sim_dir report" restart_simulation.R outfile.Rout`

Two scripts are provided with the `CTSim` package in the `exec/` directory which can be used to run simulations in batch mode using `R CMD BATCH` (see above for usage). [`run_simulation.R`](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/run_simulation.html) calls `run_sim_P` on command line arguments and starts new simulations, whereas [`restart_simulation.R`](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/restart_simulation.html) attempts to restart existing sets of simulation runs. Command line arguments must be speficied in order.

### Summarizing Simulation Runs
[`summarize_sim(sim, breaks, locs, t_window, species = NULL, land = NULL, gsad = NULL, agg_times = NULL, P_obs = list(1), sum_parms = NULL)`](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/summarize_sim.html)

[`summarize_sim_N(sim, breaks, locs, t_window, agg_times = NULL, P_obs = list(1), sum_parms = NULL, sum_func = NULL)`](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/summarize_sim_N.html)

These two functions are used to summarize the results of a single simulation (`summarize_sim`) or multiple runs from one set of simulation parameters (`summarize_sim_N`). `summarize_sim` calculates four main types of community descriptors on a given simulation run for a set of spatial and temporal units and then summarizes these descriptors across spatial and temporal units. A spatial unit is a collection of grid cells and a temporal unit is a collection of timepoints. Thus, by specifying different sets of spatial and temporal units the user can analyize a simulation at different spatial and temporal scales. [`aggregate_cells`](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/aggregate_cells.html)  is a useful function for creating sets of spatial units at different spatial scales which can be passed to the `locs` argument.

The community descriptors that the function calculates are returned as arrays in a list:

+ **bio**: richness and abundance of biologically core and transient species (based on birth rates) 
+ **occ**: richness and abundance of species in classes based on their temporal occupancy
+ **xclass**: number of species in each of four categories: biologically core - occupancy core, biologically core - occupancy transient, biologically transient - occupancy core, and biologically transient - occupancy transient (see [cross_classify](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/cross_classify.html)
+ **abun**: relative rank abundance of each species in the metacommunity

Details on how the summary function proceeds:

1. Actual species' abundances are calculated for each spatial unit in `locs` and temporal unit in `t_window`. See documentation of [`calc_abun_profile`](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/calc_abun_profile) for more information on these parameters.
2. Observed species' abundance profiles are calculated by sampling actual abundance profiles using the detection probabilities provided in `P_obs`.
3. Observed abundance profiles are used to calculate species temporal occupancy in each of the spatial units, after aggregating timepoints according to `agg_times` (see [`calc_occupancy`](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/calc_occupancy.html)).
4. Richness and abundance of species in each of the temporal occupancy categories defined by breaks are calculated in each spatial unit for each timepoint. Timepoints may be aggregated prior to calculations using the `agg_times` argument in `sum_parms`. See (see [`calc_rich_CT`](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/calc_rich_CT.html)) for details.
5. For each spatial unit, habitat values are averaged across the cells that comprise it (see [`average_habitat`](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/average_habitat)) and this average habitat type is used to determine which species are biologically core and transient (based on whether they have positive birth rates in that habitat). Richness and abundance of biologically core and transient species is then calculated for each spatial unit at each timepoint. As in the previous step, timepoints may be aggregated prior to calculations using the `agg_times` argument in `sum_parms`.
6. For each spatial unit, species are classified as core or transient based on whether they fall in the first or last occupancy interval defined by breaks. Then, the number of species classified as biologically core - occupancy core, biologically core - occupancy transient, biologically transient - occupancy core, and biologically transient - occupancy transient are counted for each spatial unit.
7. Thus far the simulation has calculated richness and abundance for each timepoint. The element `time_sum` in `sum_parms` defines which timepoints are used for summarizing richness and abundance. If `time_sum = 'none'` then all timepoints are summarized individually. Other options are `'mean'` for averaging across all timepoints in `t_window` or `'last'` for only using the last timepoint in `t_window`.
8. After timepoints have been summarized (in step 7), summary statistics are calculated across all spatial units. Minimally, this is the mean and variance, but may also include any quantiles specified in the `quants` argument of `sum_parms`.
9. The relative abundance of each species is calculated across the entire metacommunity for each set of observed abundances (from step 2) and these abundances are returned in rank order.

#### Summarizing multiple simulation runs

`summarize_sim_N` calls `summarize_sim` for each independent simulation run. In addition, it also summarizes landscape properties across the spatial units specified in `locs` using [`summarize_land`](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/summarize_land.html) for the landscapes used in each simulation run. If `sum_func` is provided, this function is used to summarize quantities across runs. The package includes a useful default summary function (`default_sum_func`) which can be passed to `sum_func`. It returns the mean, variance and qunatiles at 2.5%, 50% and 97.5%.

Simulation runs to summarize can be specified in two ways using the parameter `sim`: by passing the list of simulation results returned by `run_sim_N` or by passing the name of the directory where multiple simulation result files are saved. Simulation run files are identified by ending in `'run<i>.RData'`, where `i` is a number. This is the naming format that is automatically generated by run_sim_N.

#### Summarizing simulations with multiple sets of summary parameters
[`summarize_sim_P(run_dir = "./", parm_dir = "./", results_dir = "./Summaries/", cross_time = F)`](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/summarize_sim_P.html)

Users can save summary parameters passed to `summarize_sim_N` in a parameter file and then use `summarize_sim_P` to summarize simulation results using one or more summary paramter files. [Here](/Code/Parameters/example_summary_parameter_file.txt) is an exmaple of a summary parameter file. The function uses each summary parameter files in `parm_dir` to summarize a set simulation runs saved in the directory `run_dir`. Parameter filenames must start with `'s_'` and simulation runs should all be from one set of simulation parameters (e.g. as saved by `run_sim_N`). Summary parameter files are read into R using `source()` and should therefore be R-readable. Each summary parameter file should also have a unique `sumID` defined within it as well define objects that can be passed as parameters to `summarize_sim_N`. Parameters requiring a value are `breaks`, `locs` and `t_window`. The function will either summarize simulation results for a set time period defined in `t_window` (default) or for multiple consecutive time windows across the entire simulation period (use `cross_time=TRUE`), in which case `t_window` defines the time interval and must be a list with named elements start and stop. If `cross_time=FALSE` then two summary objects are saved to the .RData file: 

+ `sim_sum_ind`: includes a summary for each run
+ `sim_sum`: summarizes quantities across runs using the function defined in `sum_func`. 
 
If `cross_time=TRUE` then only `sim_sum_ind` is saved and `T=` is appended to the filename to denote the timestep at which the summary ends. For example, a simulation with 100 timepoints and `t_window = list(start=91, stop=100)` will conduct summaries for time windows T1 - T10, T11 - T20 ,T21 - T30 , ... T91 - T100, and save each of these summaries to a separate file. In contrast the if the same `t_window` is specified, but with `cross_time=FALSE`, only the last time window will be analyzed and only one summary file saved. 

Results are saved in the directory `results_dir` and are not returned by the function. Filenames follow the convention `'_summary.RData'`.

Two scripts are provided with the `CTSim` package in the `exec/` directory which can be used to summarize simulations in batch mode using `R CMD BATCH` (see below for usage). [`summarize_simulation.R`](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/summarize_simulation.html) calls `summarize_sim_P` on command line arguments with `cross_time = FALSE`, whereas [`summarize_simulation_cross_time.R`](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/summarize_simulation_cross_time.html) calls `summarize_sim_P` with `cross_time = TRUE`. Command line arguments must be speficied in order.

Usage:

Place the R script in the working directory, then use:

`R CMD BATCH "--args ncores parm_dir results_dir sim_dir report" summarize_simulation.R outfile.Rout`

`R CMD BATCH "--args ncores parm_dir results_dir sim_dir report" summarize_simulation_cross_time.R outfile.Rout`


## Parameter Files for Running Simulations
`/Code/Parameters/`

See [example_parameter_file.txt](./Code/Parameters/example_parameter_file.txt) for an example of all possible parameters that can be provided for running a simulation and see [baseline_parameter_file.txt](./Code/Parameters/baseline_parameter_file.txt) for the set of parameters used as a basis for [experiments](./wiki/Experiments). Required and optional parameters are described below:

The following parameters are **required** for simulation and must be present or the simulation will fail:

 + **`dimX`**: x-dimension of landscape
 + **`dimY`**: y-dimension of landscape
 + **`S_A`**: number of specialist species on habitat type A
 + **`S_B`**: number of specialist species on habitat type B
 + **`m_rates`**: vector of mortality rates on preferred and	non-preferred habitats
 + **`r_rates`**: vector of recruitment rates on preferred	and non-preferred habitats
 + **`K`**: carrying capacity of cells
 + **`nsteps`**: number of timesteps to simulate
 + **`nruns`**: number of independent simulations to run

For further details on `S_A`, `S_B`, `m_rates` and `r_rates` see `make_species`. For further details on `K` see `populate_landscape`.

The following parameters are optional and more information can be found in the documentation on the functions they are passed to:

 + Parameters passed to [`make_landscape`](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/make_landscape.html)
  + **`vgm_dcorr`**: distance at which habitat values become uncorrelated
  + **`vgm_mod`**: variogram model controling spatial autocorrelation of habitat values
  + **`habA_prop`**: proportion of landscape that comprised of habitat type A
 + Parameters passed to [`make_species`](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/make_species.html)
  + **`S_AB`**: number of generalist species
  + **`dist_b`**: list defining the distribution from which	species' birth rates are sampled
  + **`dist_d`**: list defining the distribution from which species' dispersal rates are sampled. Must contain character string named `type`.
  + **`dist_v`**: list defining the distribution from which	species' movement rates are sampled. Must contain character string	named `type`.
  + **`dist_gsad`**: list defining distribution from which global species abundances are sampled or `'b_rates'`, indicating that the gsad should match species birth rates in their preferred habitat
 + Parameters passed to [`populate_landscape`](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/populate_landscape.html)
  + **`prop_full`**: proportion of the landscape's carrying capacity that should initially contain individuals
  + **`init_distribute`**: character string indicating how	individuals should be initially distributed across the landscape
  + **`cells_distribute`**: if `init_distribute` is `'designated'`, a matrix giving the locations of cells in which to place propagules
 + Parameters passed to [`run_sim`](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/run_sim)
  + **`d_kernel`**: list defining the shape of the dispersal	kernel of new propagules
  + **`v_kernel`**: list defining the shape of the movement	kernel of established individuals
  + **`imm_rate`**: immigration rate- probability than an empty space will be colonized by a migrant from outside	the metacommunity
 + Parameters passed to [`run_sim_N`](http://htmlpreview.github.com/?https://github.com/hurlbertlab/core-transient-simulation/blob/master/Code/HTML/run_sim_N)
  + **`save_steps`**: vector of timesteps to save in each simulation
  + **`simID`**: character string that identifies simulations	run on this set of parameters

## Useful Scripts
`/Code/Scripts/`

+ [make_package.R](/Code/Scripts/make_package.R): R code for compiling the package. Should be run to create new package versions and update documentation.
+ [make_parmfiles.R](Code/Scripts/make_parmfiles.R): R code for generating simulation parameter files used in the experiments
+ [visualize_simulations.R](Code/Scripts/visualize_simulations.R): R code for analyzing experiment results and generating visualizations
+ [submit_runs.txt](Code/Scripts/submit_runs.txt): fragments of shell code used to submit simulation batches to computing cluster
+ [debugging.R](Code/Scripts/debugging.R): fragments of R code used in debugging CTSim functions



## Simulation Results
`/Results`

This directory holds the summary data files and visualizations from Experiments performed with CTSim. Actual data files are too large to store in this repository, but are saved in an untracked 'Data' directory here. 

See [wiki page](/wiki/Experiments) for details on parameters used in each experiment and naming conventions. In general, files are named according to parameters that were varied. A '_' separates different parameters while a '-' associates a parameter with a value. For example a file including the string: 'd-g1_v-a0.5' indicates that the data refers to a simulation where the juvenile dispersal kernel was set to g1 (gaussian with distance 1) and the adult movement kernel was set to a0.5 (adjacent cell with probability 0.5). 
 + **CONV**: Initial runs testing for long-term convergence of simulation under different grid sizes (32x32 vs 64x64) and dispersal types (d) or birth rates (b). These were just used to determine that 200 timestes were adequate and shouldn't be used for further analysis.
 + **EXP1**: Experiment 1 evaluates effects of juvenile dispersal (d), adult movement (v), and habitat spatial autocorrelation (dcorr)

See [wiki page](/wiki/Experiments) for specific parameters used in summaries or see the Parameters directory that corresponds to the experiment.
  





