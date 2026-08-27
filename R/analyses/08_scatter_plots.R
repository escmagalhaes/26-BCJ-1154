### Scatter plots analyses ###

#Import input
corr_results<-qs2::qs_read(file.path(analyses_export_path,"correlations_analyses.qs2"))

#Create local object
data_corr<-corr_results$combined_correlations

#Create list of selected correlations 
sc_dataset_names<-c(
  "all_chr3_alt_neg"      = 'Chr3-WT',
  "all_chr3_alt_MECOMR"   = 'MECOM-R',
  "all_chr3_alt_Chr3abn"  = 'Chr3-abn',
  "high_chr3_alt_neg"     = 'Chr3-WT High',
  "low_chr3_alt_neg"      = 'Chr3-WT Low',
  'high_chr3_alt_Chr3abn' = 'Chr3-abn High',
  "low_chr3_alt_Chr3abn"  = 'Chr3-abn Low'
)
sc_datasets      <-names(sc_dataset_names)
comb_sc_datasets <-combn(sc_datasets,2,simplify=FALSE)

#Create varying parameters list 
sc_params_list<-lapply(comb_sc_datasets, function(x) {
  list(
    datasetA=x[1],
    datasetB=x[2],
    datasetA_name=sc_dataset_names[[x[1]]],
    datasetB_name=sc_dataset_names[[x[2]]]
  )
})
names(sc_params_list)<-sapply(comb_sc_datasets, function(x) {
  paste0(sc_dataset_names[x[1]]," vs ",sc_dataset_names[x[2]])
})

#Create fixed parameters list 
fixed_sc_params<-list(
  comb_corr     = data_corr,
  query_protein = 'MECOM'
)

#Create all Scatter plots with a single call
sc_res<-lapply(sc_params_list,function( p ) {
  args<-c(fixed_sc_params, p )
  do.call(plot_corr_scatter, args )
})

### Save analyses ###
qs2::qs_save(sc_res,file.path(analyses_export_path,"scatter_plots.qs2"))

