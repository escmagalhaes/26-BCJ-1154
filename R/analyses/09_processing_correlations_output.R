### Processing correlations output ###

#Import inputs
corr_results<-qs2::qs_read(file.path(analyses_export_path,"correlations_analyses.qs2"))

#Create local objects
data_corr         <-corr_results$combined_correlations
corr_ptn_name_map <-corr_results$ptn_name_mapping

#Preprocess correlation data for downstream analyses in a single call
fixed_preprocess_cor_params<-list(
  cor_df     = data_corr,
  df_ptn_map = corr_ptn_name_map
)
cor_sig_datasets<-list(
  'MECOM-R'       = 'all_chr3_alt_MECOMR',
  'Chr3-abn High' = 'high_chr3_alt_Chr3abn',
  'Chr3-WT High'  = 'high_chr3_alt_neg'
  )
corr_res_proc<-lapply(cor_sig_datasets,function(x) {
  do.call(preprocess_corr_sig,c(fixed_preprocess_cor_params,list(datasets=x)))
})

### Save analyses ###
qs2::qs_save(
  list(
    data_corr = data_corr,
    corr_proc = corr_res_proc
  ),file.path(analyses_export_path,"processed_corr_output.qs2"))

  
            


