# core-transient-simulation

This project is a spatially explicit simulation of metacommunity dynamics that models habitat filtering and dispersal.

## Simulation Description
The simulation occurs on a regular grid of cells (termed the 'landscape') and each cell has discrete habitat type, 'A' or 'B'. Each cell contains a 'community' of individuals, which is fixed at a pre-determined carrying capacity. Communities may contain fewer individuals than the carrying capacity, but not more. The species pool consists of three types of species- generalists, and habitat specialists on habitat type A or B. Generalists produce offspring in both habitat types while specialists only produce offspring in their preferred habitat. After initializing a simulation landscape and species pool, the simulation proceeds through repeated iteration of four processes: birth - dispersal - establishment - death. Once established, individuals vacate their place in the community only through dispersal or death- the simulation does not model competitive displacement. 


## Simulation Intialization
A simulation requires three objects: a landscape, a species pool, and a global species abundance distribution (hereafter, gsad). The functions used to initialize these objects are described below.

### Initializing a landscape:
`make_landscape(x, y, mod, d, prop, draw_plot)`

A landscape is a raster layer whose values are -1 or 1, corresponding to habitat types 'A' and 'B', respectively. 
When a landscape is initialized, the user can specify its size (`x,y`), the proportion of the cells which should belong to habitat type 'B' (`prop`), a variogram model used to define the spatial autocorrelation of habitat values (`mod`), and the distance at which habitat values become uncorrelated (e.g. the range of the variogram model: `d`). Variogram models are implemented by the `vgm()` function in gstat. Use `show.vgms()` to see available models.

**Required parameters and defaults**

Only the dimensions of the landscape are required. By default, the function will return a grid with 50% of habitat type 'A' and a exponential variogram model with partial sill = 1 and range = 1/3 of the grids smallest dimension.

### Initializing a species pool
`make_species(S_A, S_B, S_AB, dist_b, m, r, dist_d, dist_v)`

A species pool is a 3-dimensional array of species vital rates. The first dimension specifies the species. The second dimension defines in which habitat the rate applies ('A' or 'B') and the third dimension specifies the type of rate:

1. **b**: birth rate- number of offpsring produced by established individual per timestep.
2. **m**: mortality rate- probability that established individual dies during a time step.
3. **r**: recruitment rate- probability than a propagule arriving at an open space will become an established individual
4. **d**: dispersal rate- expected distance that a propagule will disperse from its origination point
5. **v**: movement rate- expected distance that an established individual will move from its cell

Birth rates are positive in a species' preferred habitat and 0 elsewhere. Generalists prefer both habitat types equally. Dispersal rates control how newly produced propagules move away from their cell of origin. Movement rates control how established individuals move from their current cell. Movement, mortality, and recruitment rates can be set to differ systematically between preferred and non-preferred habitats.

### Initializing a GSAD
`make_sad(N_S, distribution)`

A GSAD is a vector of species relative abundances on the "mainland". These define the relative probabilities that an immigrant from outside the landscape will belong to each species. Implemented distributions are 'same', 'uniform', 'power', 'logseries', 'lognormal', and 'poisson'. See the help file on `make_sad` for further details.

## Simulation Operation
`populate_landscape(land, species, gsad=NULL, K, distribution=NA, p=NA, which_cells=NULL)`
`run_sim(steps, metacomm, land, species, gsad, d_kernel=NULL, v_kernel=NULL, imm_rate=NA,...)`

Once a landscape, species pool and GSAD have been generated, a simulation can be run. First, the landscape must be populated by species using the function `populate_landscape`. This function has several options for defining the number of individuals in each cell (`K`), the proportion of this carrying capacity which should be initially filled (`p`), as well as how individuals should be  distributed across the landscape (`distribution`). The `which_cells` option allows the user to specify the exact cells where individuals should be placed and the number to place in each cell. Individuals are drawn probabilistically from the specified GSAD (`gsad`), which defaults to the same probability for all species.

The populated landscape is a metacommunity object which takes the form of a matrix of lists. The `x` and `y` positions define the location of each community on the landscape, while the list defines which individuals are present. The length of these lists are fixed for the remainder of the simulaton to the carrying capacity defined in `K`. Integers refer to species identities and a `0` indicates that no individual is present. For example, if `metacomm[1,2]` is equal to `list(3,3,2,10,0,2,1)`, then the cell in the first row and second column can hold 7 individuals and currently contains 1 individual of species 1, 2 individuals of species 2 and 3, 1 individual of species 10, and one empty space. 

After a metacommunity has been generated, a simulation can be run for a fixed length of time (`steps`) using the `run_sim` function. In addition to the objects defining the metacommunity (`metacomm`), landscape (`land`), species pool (`species`), and GSAD (`gsad`), the function also requires the probability that an empty space will be colonized by an individual from the "mainland" (`imm_rate`) as well as information about how propagules disperse (`d_kernel`) and individuals move after they have established in a cell (`v_kernel`). Three dispersal modes are implemented: half-gaussian, adjacent cell, and uniform. In all cases, the direction of dispersal is random (isotropic) and the expected dispersal distances are determined by dispersal and movement rates for individual species defined in the species pool object. See the help file on `get_dispersal_vec` for additional information on specifying dispersal and movement parameters. 

A simulation runs by iteratively calling the function `run_timestep`, which defines the operations that occur in a single timestep. This function calls the main process functions `die`, `reproduce`, `disperse`, `establish` to progress a metacommunity through one simulation timestep. Each of these functions are described in detail in their help files. The order of operations is as follows:
1. **Death** : Established individuals in each community experience	probabilistic mortality according to species- and habitat-specific mortality rates provided in `species`.
2. **Birth** : Established individuals in each community produce propagules according to species- and habitat-specific birth rates provided in `species`.
3. **Dispersal** : New propagules from each community disperse across the landscape (`land`) according to species-specific dispersal rates provided in `species`. The parameter `d_kernel` specifies the dispersal kernel for new propagules.
4. **Movement** : Established individuals in each community move from their current cell with species- and habitat-specific movement rates provided in `species`. The parameter `v_kernel` specifies the dispersal kernel for previously established individuals.
5. **Establishment** : Empty spaces in each community are colonized by either a migrant from outside the community with probability `imm_rate` or by an individual selected at random from the pool of new propagules and moving individuals that arrived in the cell. External migrants are chosen probabilistically from the relative abundances given in `gsad`. 



## Multiple Simulation Runs




## Parameter Files for Running Simulations





## Summarizing Simulation Runs



