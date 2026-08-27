### RUN FULL WORKFLOW ###

## REMEMBER TO OPEN CYTOSCAPE SOFTWARE BEFORE STARTING ##

#1)Load all configuration
cat("Loading all configuration...\n")
source("R/setup/load_all.R")

#2)Run analysis
cat("Running all analyses...\n")
#Load input handling scripts
lapply(list.files("R/analyses",full.names=TRUE),source)
cat("Analysis complete!\n")

#3)Export Results
cat("Compiling and Exporting results...\n")
lapply(list.files("R/workflows",full.names=TRUE),source)
cat("Results exported!\n")




