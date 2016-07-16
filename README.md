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

A species pool is a 3-dimensional array of species vital rates. The first dimension specifies the species. The second dimension specifies the type of rate:
1. **b**: birth rate- number of offpsring produced by established individual per timestep.
2. **m**: mortality rate- probability that established individual dies during a time step.
3. **r**: recruitment rate- probability than a propagule arriving at an open space will become an established individual
4. **d**: dispersal rate- expected distance that a propagule will disperse from its origination point
5. **v**: movement rate- expected distance that an established individual will move from its cell
The third dimension allows different rates to be specified for a species in each habitat type.




## Simulation Operation



## Multiple Simulation Runs




## Parameter Files for Running Simulations





## Summarizing Simulation Runs



