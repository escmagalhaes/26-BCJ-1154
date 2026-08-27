### Configuration for importing and exporting files and for function parameters ###

#Get the project root
project_root<-getwd()

#Get the project folder path
project_path<-file.path(project_root)

#Define input path
input_path<-file.path(project_root,"inputs")

#Define output paths 
input_path           <-file.path(project_path,"inputs")
analyses_export_path <-file.path(project_path,"output","Analyses")
data_export_path     <-file.path(project_path,"output","Exported")

#Define biocache path
biocache_path<-file.path(project_root,"bio_cache")

#Create output folders (if they don't exist)
output_directories<-c(analyses_export_path,data_export_path,biocache_path)
for (directory in output_directories) {
  if (!dir.exists(directory)) dir.create(directory,recursive=TRUE)
}

#Check that all input files are present
input_files<-c(
  "data1.csv",
  "data2.csv",
  "data3.csv",
  "data4.csv",
  "data5.csv",
  "data6.csv",
  "data7.csv",
  "data8.csv",
  "data9.csv",
  "data10.csv"
  )
files_missing<-input_files[!file.exists(file.path(input_path,input_files))]
if (length(files_missing) > 0) {
  stop("\nThese input files are missing: ",paste(files_missing,collapse=", "))
} else {
  cat("\nAll input files found.\n")
}

## STRING DB CONFIGURATION ##
string_version<-"11.5"
string_species<-9606
string_score_threshold<-400
string_neighbor_min_connections<-4
string_neighbor_score_threshold<-700

#Define bio_cache root for STRINGdb
string_cache_path<-file.path(biocache_path,"STRINGdb")
if (!dir.exists(string_cache_path)) dir.create(string_cache_path,recursive=TRUE)

## Synapse CONFIGURATION ##
.synapse_token<-""
if (.synapse_token=="") { 
  warning(
    "Synapse Token variable is empty. This will prevent Synapse workflow from 
    occuring downstream, but all other workflows will be unaffected.
    Remember to create a Synapse account and a Personal Access Token (PAT). 
    PAT should be a long string of characters/numbers, not the name you gave to it.
    For more details, check out: 
    https://docs.synapse.org/synapse-docs/managing-your-account#ManagingYourAccount-PersonalAccessTokens"
  )
} else {
  if (!requireNamespace("synapser",quietly=TRUE)) {
    tryCatch({
      if (!requireNamespace("remotes",quietly=TRUE)) install.packages("remotes")
      remotes::install_cran(
        "synapser",
        version = "2.1.5.356",
        repos   = c("http://ran.synapse.org","https://cloud.r-project.org")
      )
    }, error = function(e) {
      message("synapser installation failed: ",conditionMessage(e))
    })
  }
  synapse_cache_path<-file.path(biocache_path,"Synapse")
  if (!dir.exists(synapse_cache_path)) dir.create(synapse_cache_path,recursive=TRUE)
}

## Print paths and files being used to double check that everything is in place ##
name_grid<-c(
  "Project root folder",
  "Inputs folder",
  "Analyses folder",
  "Exported data folder",
  "Biocache folder"
)

path_grid<-c(
  project_root,
  input_path,
  analyses_export_path,
  data_export_path,
  biocache_path
)

cat("\n=========================== Paths Being Used =========================================\n")
for(i in seq_along(name_grid)){ cat(sprintf("\n%-29s: %s\n",name_grid[ i ],path_grid[ i ])) }

cat("\n=========================== Input Files =============================================\n")
for( f in input_files){ cat(sprintf("\n%s\n",file.path(input_path, f ))) }
cat("\nYOU ARE GOOD TO GO!\n\n")



