### Load and restore packages for reproducibility ###

#Install renv if not already installed
if (!requireNamespace("renv",quietly=TRUE)) { install.packages("renv") }
suppressPackageStartupMessages(library(renv))

#Restore packages from lockfile
if (!file.exists("renv.lock")) {
  if (exists("first_time_setup") && isTRUE(first_time_setup)) {
    cat("\nFirst time setup detected. Skipping restore.\n")
  } else {
    stop("\nrenv.lock not found. If setting up for the first time, run setup_run_once.R first.\n")
  }
} else {
  renv::restore()
  cat("\nPackages restored from renv.lock\n")
}

#Load used packages for the workflow
source("R/setup/packages.R")

#Load packages and track missing ones
loaded<-c()
failed<-c()
for (pkg in packages$load) {
  success<-tryCatch({
    suppressPackageStartupMessages(library(pkg,character.only=TRUE))
    TRUE
  },error=function(e) FALSE)
  if (success) loaded<-c(loaded,pkg) else failed<-c(failed,pkg)
}

#Summary messages
cat("\n Packages loaded:",paste(loaded,collapse=", "),"\n")
if (length(failed) > 0) {
  warning("Failed to load: ", paste(failed, collapse=", "),
          "\nThese may be missing from renv.lock or require manual installation.")
}
cat("\nEnvironment ready!\n")


