### Pathway Enrichment Analysis ###

#Import input
proc_corr <-qs2::qs_read(file.path(analyses_export_path,"processed_corr_output.qs2"))
ref_sets  <-qs2::qs_read(file.path(analyses_export_path,"processed_reference_sets.qs2"))

#Create local object
proc_corr_res <-proc_corr$corr_proc
pfg_ref_data <-ref_sets$pfg_list

#Generate pathway analyses in a single call
fixed_paths_params<-list(
  libs     = c('react'),
  pfg_list = pfg_ref_data,
  analysis = 'all'
)
sig_cor_ptn    <-lapply(proc_corr_res, function(x) x$ptn_class)
paths_corr_res <-lapply(sig_cor_ptn,function(x) {
  do.call(run_enriched_paths,c(fixed_paths_params,list(proteins=x)))
})

## Generate pathway heatmaps ##
paths_list<-lapply(paths_corr_res, function(x) x$all$enriched_paths_tab)
fixed_path_ht_params<-list(
  df_list          = paths_list,
  top_n_paths      = 20,
  row_fontsize     = 8,
  col_fontsize     = 8,
  legend_title     = "Combined Score",
  legend_direction = 'horizontal'
)

#PFG Heatmap
pfg_ht<-do.call(path_heatmap,c(
  fixed_path_ht_params,
  list(
    analysis   = "PFG",
    plot_title = "PFG Pathway Enrichment across datasets")
))
pfg_ht_plot<-as.ggplot(function() draw(
  pfg_ht$plot,
  heatmap_legend_side = "bottom",
  padding             = unit(c(2,15,2,2),"mm"))
)

#Reactome Heatmap
react_ht<-do.call(path_heatmap,c(
  fixed_path_ht_params,
  list(
    analysis      = "EnrichR",
    max_str_width = 45,
    libs          = 'react',
    plot_title    = "EnrichR Pathway Enrichment across datasets")
))
react_ht_plot<-as.ggplot(function() draw(
  react_ht$plot,
  heatmap_legend_side = "bottom",
  padding             = unit(c(2,20,2,2),"mm"))
)

## Compile pathway analysis tables to export##
paths_tabs<-c(paths_list,list(pfg_ht$table,react_ht$table))
names(paths_tabs)<-c(
  'MECOM-R Pathways related to Significantly Correlated Proteins'
  ,'Chr3-abn-high Pathways related to Significantly Correlated Proteins'
  ,'Chr3-WT-High Pathways related to Significantly Correlated Proteins'
  ,'Combined Scores of PFG pathway analyses ranked'
  ,'Combined Scores of Reactome analyses ranked'
)

### Save analyses ###
qs2::qs_save(list(
  pfg_ht     = pfg_ht_plot,
  react_ht   = react_ht_plot,
  paths_tabs = paths_tabs
  ),file.path(analyses_export_path,"enriched_paths.qs2"))
