### Dataset of Correlation data ###

#Import inputs
corr_data<-qs2::qs_read(file.path(analyses_export_path,"correlations_analyses.qs2"))

#Create local objects
data_corr<-corr_data$combined_correlations

#Create folders and paths to export results
corr_dataset_path<-file.path(data_export_path,'corr-dataset')
if (!dir.exists(corr_dataset_path)) { dir.create(corr_dataset_path,recursive=T) }

#Extract Correlation table and adjust names for export
corr_tab_labels<-c(
  'Chr3-WT High'  = 'high_chr3_alt_neg',
  'Chr3-WT Low'   = 'low_chr3_alt_neg',
  'MECOM-R'       = 'all_chr3_alt_MECOMR',
  'Chr3-abn High' = 'high_chr3_alt_Chr3abn',
  'Chr3-abn Low'  = 'low_chr3_alt_Chr3abn'
)

corr_tab<-lapply(corr_tab_labels, function( lbl ){
  data_corr|> dplyr::filter(dataset_label == lbl ) |> dplyr::select(
    c(-dataset,-query_protein,-mean_expr_query_protein,
      -mean_expr_protein,-pair_id,-dataset_label))  |>
    dplyr::mutate(dataset_name = names(corr_tab_labels)[corr_tab_labels == lbl ]) |> 
    dplyr::select(c(dataset_name,everything()))
})

### Export dataset ###
export_tabs(
  df_list   = corr_tab,
  filepath  = corr_dataset_path,
  filename  = 'corr-dataset',
  tab_title = paste("Spearman's Correlations in",names(corr_tab)," patients")
)

