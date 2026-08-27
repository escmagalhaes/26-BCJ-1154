### Define all packages used in the project ###

packages<-list(install=c(
  
  #Package management 
  "remotes",
  
  #First install remotes dependencies for rppx
  "ddsjoberg/gtforester",
  "https://cran.r-project.org/src/contrib/Archive/progenyClust/progenyClust_1.2.tar.gz", #progenyClust has been removed from CRAN and archived 
  "bioc::limma",
  "bioc::ComplexHeatmap",
  "bioc::MultiAssayExperiment",
  "bioc::cBioPortalData",
  "bioc::EnhancedVolcano",
  "bioc::RCy3",
  "bioc::STRINGdb",
  
  #Data analysis workflows
  "escmagalhaes/rppx@v0.1.0-beta",
  
  #Stats / modeling
  "broom",
  "caret",
  
  "rstatix",
  "survival",
  "survminer",
  "metafor",
  
  #Misc / utility
  "rstudioapi",
  "digest",
  "labelled",
  "corrr",
  "broom.helpers",
  "webshot2",
  "plyr",
  "reshape2",

  #Data wrangling
  "dplyr",
  "forcats",
  "janitor",
  "tibble",
  "tidyr",
  "readr",
  "sumExtras",  
  
  #Plotting / visualization
  "corrplot",
  "cowplot",
  "ggplot2",
  "ggpubr",
  "ggplotify",
  "ggsci",
  "ggVennDiagram",
  "gridExtra",
  "pheatmap",
  "RColorBrewer",
  "matlab",        # for jet.colors()
  
  #Tables / reporting
  "gt",
  
  "gtsummary",
  
  #Data Acquisition/Enrichment/Pathway analysis/Networks
  "enrichR",
  "HGNChelper",
  "igraph",
  #"synapser",
  
  #Parallel / reproducibility
  "doParallel",
  "doRNG",
  "foreach",
  
  #File input/output
  "openxlsx2",
  "readxl",
  "qs2"
),
load=c(
  "remotes",
  "gtforester","progenyClust","limma","ComplexHeatmap","MultiAssayExperiment","cBioPortalData","RCy3","EnhancedVolcano","STRINGdb", 
  "rppx",
  "broom", "caret",  "rstatix", "survival", "survminer","metafor",
  "rstudioapi", "digest", "labelled", "corrr", "broom.helpers", "webshot2", "plyr", "reshape2",
  "dplyr", "forcats", "janitor", "tibble", "tidyr", "readr", "sumExtras", 
  "corrplot", "cowplot",  "ggplot2", "ggpubr",
  "ggplotify", "ggsci", "ggVennDiagram", "gridExtra", "pheatmap","RColorBrewer", "matlab",
  "gt",  "gtsummary",
  "enrichR", "HGNChelper",  "igraph",#"synapser",
  "doParallel", "doRNG", "foreach",
  "openxlsx2", "readxl","qs2"
)
)
