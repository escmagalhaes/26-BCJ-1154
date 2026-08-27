### Tables of pairwise tests related to boxplots ###

#Import inputs
bxplts<-qs2::qs_read(file.path(analyses_export_path,"boxplots.qs2"))

#Create folders and paths to export results
bxplt_tests_path<-file.path(data_export_path,'boxplot-tests-dataset')
if (!dir.exists(bxplt_tests_path)) { dir.create(bxplt_tests_path,recursive=T) }

#Extract test tables
bxp_tab_list<-lapply(bxplts, function(x) {
  tab<-x$test_table[,!colnames(x$test_table) %in% c("y.position","groups"),drop=FALSE]
  if(nrow(tab)==1) tab<-NULL
  tab
})
names(bxp_tab_list)<-sub("^[^_]+_","",names(bxp_tab_list))
bxp_tab_list<-bxp_tab_list |> purrr::compact() #Remove all NULL items from the list

### Export dataset ###
export_tabs(
  df_list   = bxp_tab_list,
  filepath  = bxplt_tests_path,
  filename  ='boxplot-tests-dataset',
  tab_title = paste("Pairwise comparisons related to",names(bxp_tab_list),"Boxplots"))


