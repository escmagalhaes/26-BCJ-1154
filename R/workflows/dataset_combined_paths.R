### Tables with pathway analysis ###

#Import inputs
corr_paths<-qs2::qs_read(file.path(analyses_export_path,"enriched_paths.qs2"))

#Create local objects
paths_tabs<-corr_paths$paths_tabs

#Create folder and path to export results
combined_paths_path<-file.path(data_export_path,'combined-paths-dataset')
if (!dir.exists(combined_paths_path)) { 
  dir.create(combined_paths_path,recursive=T) 
}

### Export Tables ###
export_tabs(
  df_list   = paths_tabs,
  filepath  = combined_paths_path,
  filename  ='combined-paths-dataset',
  tab_title = paste('Pathways Related to',names(paths_tabs))
)

