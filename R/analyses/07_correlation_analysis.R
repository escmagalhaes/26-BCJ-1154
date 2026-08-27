### Correlation analyses ###

#Import input datatables
main_data <-qs2::qs_read(file.path(analyses_export_path,"processed_main_data.qs2"))
ref_sets  <-qs2::qs_read(file.path(analyses_export_path,"processed_reference_sets.qs2"))

#Create local objects
df_full      <-main_data$main_df
df_normal    <-ref_sets$normal_data
all_proteins <-ref_sets$all_proteins

#Run correlations
corr_res<-run_corr_pipeline(
  df_input             = df_full,
  df_normal            = df_normal,
  all_proteins         = all_proteins,
  proteins_of_interest = 'MECOM',
  grouping_variable    = 'chr3_status',
  var_labs             = 'chr3_alt',
  df_subsetting        = TRUE,
  ncol_grid            = 1,
  p_adjust_method      = "none"
)

### Save analyses ###
qs2::qs_save(corr_res,file.path(analyses_export_path,"correlations_analyses.qs2"))




