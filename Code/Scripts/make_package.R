## This script is used for building and maintaining the CTSim package

library(devtools)
library(roxygen2)

setwd('C:/Users/jrcoyle/Documents/Research/CT-Sim/GitHub/Code/')

current_code = as.package('CTSim')

# Load functions
load_all(current_code)

# Update documentation
document(current_code)

