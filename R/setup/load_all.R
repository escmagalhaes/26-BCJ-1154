### Load scripts for downstream analysis ###

#Load environment and packages first
source("R/setup/renv.R")    

#Load Configuration script
source("R/setup/config.R")

#Load helper functions
lapply(list.files("R/helper_functions",full.names=TRUE),source)


