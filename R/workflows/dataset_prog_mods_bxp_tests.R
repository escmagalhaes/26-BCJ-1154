### Tables of pairwise tests related to MECOM/CEBPA prognostic models boxplots ###

#Import inputs
prog_mods<-qs2::qs_read(file.path(analyses_export_path,"prognostic_models.qs2"))

#Create folders and paths to export results
prog_mods_tests_path<-file.path(data_export_path,'prog-mods-bxp-tests-dataset')
if (!dir.exists(prog_mods_tests_path)) { dir.create(prog_mods_tests_path,recursive=T) }

#Extract test tables
bxps<-prog_mods$boxplots
bxp_prog_tab_list<-unlist(lapply(bxps, function(x)  lapply( x, function(y) { 
  #Remove uncecessary columns
  y$test_table<-
    y$test_table[,!colnames(y$test_table) %in% c("y.position","groups"),drop=FALSE]
  return(y$test_table)
})),recursive=FALSE)
names(bxp_prog_tab_list)<-gsub("\\..*", "",names(bxp_prog_tab_list))

### Export dataset ###
export_tabs(
  df_list   = bxp_prog_tab_list,
  filepath  = prog_mods_tests_path,
  filename  = 'prog-mods-bxp-tests-dataset',
  tab_title = paste("Pairwise comparisons related to",names(bxp_prog_tab_list),"Boxplots")
)


