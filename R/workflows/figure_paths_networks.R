### Analyses of Significantly correlated proteins ###

#Import inputs
proc_corr  <-qs2::qs_read(file.path(analyses_export_path,"processed_corr_output.qs2"))
corr_paths <-qs2::qs_read(file.path(analyses_export_path,"enriched_paths.qs2"))
ref_sets   <-qs2::qs_read(file.path(analyses_export_path,"processed_reference_sets.qs2"))
main_data  <-qs2::qs_read(file.path(analyses_export_path,"processed_main_data.qs2"))

#Create local objects
meta_tab      <-ref_sets$meta_table
df_full       <-main_data$main_df
data_corr     <-proc_corr$data_corr
proc_corr_res <-proc_corr$corr_proc

#Extract plots for export
pfg_ht_plot   <-corr_paths$pfg_ht
react_ht_plot <-corr_paths$react_ht

#Create folder and path to export results
paths_networks_path<-file.path(data_export_path,'paths-networks')
if (!dir.exists(paths_networks_path)) { 
  dir.create(paths_networks_path,recursive=T) 
}

#Extract mappings with proteins and metadata
corr_sig_ptn_map<-lapply(proc_corr_res, function(x) {
  x$mapping_df |> dplyr::arrange(desc(across(contains("abs_corr"))))
})

## Heatmap with protein expression ##
#Get list with protein names and turn into a dataframe for heatmap row anno
cor_ptn_list<-lapply(corr_sig_ptn_map, function(x) {
  x |> dplyr::pull(protein) |> unique() 
})
ovelap_ptn<-ComplexHeatmap::make_comb_mat(cor_ptn_list) |> 
  get_upset_ptns() |> stack() |> dplyr::distinct() |> 
  purrr::set_names(c("protein","overlap")) |>
  dplyr::mutate(overlap = factor(overlap, levels = unique(overlap)))
ht_ptn<-ovelap_ptn |> dplyr::filter(
  !overlap %in% c("MECOM-R","Chr3-abn High","Chr3-WT High")) |> 
  dplyr::mutate(overlap = droplevels(overlap))

#Get dataframe with proteins and columns annotation variables
ht_covars<-c('MECOM_gp_n2','chem2','vtx','aml_gp','cyto_risk')
ht_covar_labels<-as.vector(meta_tab[match(ht_covars,meta_tab$variable),]$label)
ht_covar_labels[is.na(ht_covar_labels)]<-"MECOM levels"

df_sub<-df_full |> dplyr::filter(
  !MECOM_gp_n2_chr3_status %in% c("Chr3-WT Low","Chr3-abn Low")) |>
  dplyr::select(upi,chr3_status,all_of(ht_covars),all_of(ht_ptn$protein)) 

#Plot heatmap
ht_corr<-plot_heatmap(
  df=df_sub,
  proteins=ht_ptn$protein,
  grouping_var='chr3_status',
  grouping_var_label="Chromosome 3 status",
  col_anno=TRUE,
  covars=ht_covars,
  covar_labels=ht_covar_labels,
  row_anno=TRUE,
  row_anno_df=ht_ptn,
  row_anno_protein_col='protein',  
  rown_anno_vars='overlap',
  rown_anno_vars_labels='Overlaps',
  covar_colors_override=list(
    'Cytogenetic Risk'=c(
      'favorable'='red',
      'intermediate'=pal_jco()(2)[1],
      'unfavorable'=pal_jco()(2)[2]
      ),
    'MECOM levels'=c(
      'High'=pal_jco()(2)[2],
      'Low'=pal_jco()(2)[1]
    )
    ),
  colors_grouping_var=c(pal_jco()(2),'red2'),
  colors_NA=pal_jco()(3)[3],
  colors_anno_row=c(pal_jco()(2)[2:1],'red2',brewer.pal(7,"Dark2")[c(7,1,2)]),
  plot_title="Significantly correlated proteins",
  cluster_rows=TRUE,
  cluster_cols=FALSE,
  show_rownames=TRUE,
  show_colnames=FALSE,
  annotation_names_row=FALSE
)

## Protein Networks ##
# Adjust dataset for protein networks generation #
corr_sig_ptn_map_net<-corr_sig_ptn_map
  
#Get list with final names and turn into a dataframe 
cor_ptn_names_list<-lapply(corr_sig_ptn_map_net, function(x) {
  x |> dplyr::pull(final_names) |> unique() 
})
ovelap_net_names<-stack(get_upset_ptns(ComplexHeatmap::make_comb_mat(cor_ptn_names_list)))
colnames(ovelap_net_names)<-c("protein","overlap")
ovelap_net_names$overlap<-factor(ovelap_net_names$overlap,levels=unique(ovelap_net_names$overlap))

#For Chr3-WT, limit the number of proteins to 50, but include all overlaps
chr3_overlaps<-ovelap_net_names |> dplyr::filter(
  overlap %in% c(
    "MECOM-R & Chr3-abn High & Chr3-WT High",
    "MECOM-R & Chr3-WT High",
    "Chr3-abn High & Chr3-WT High")) |> dplyr::pull(protein)

corr_sig_ptn_map_net$`Chr3-WT High`<-corr_sig_ptn_map_net$`Chr3-WT High`|>
  dplyr::filter(final_names %in% chr3_overlaps) |>
  dplyr::bind_rows(
    corr_sig_ptn_map_net$`Chr3-WT High` |>
      dplyr::filter(!final_names %in% chr3_overlaps) |>
      dplyr::slice_head(n = 50 - length(chr3_overlaps)))

#Adapt names to match mapping_df datasets
names(corr_sig_ptn_map_net)<-c(
  'all_chr3_alt_MECOMR',
  'high_chr3_alt_Chr3abn',
  'high_chr3_alt_neg'
)

## REMEMBER TO OPEN CYTOSCAPE SOFTWARE NOW ##
#Generate Protein networks sequentially
for ( nm in names(corr_sig_ptn_map_net)) {
  
  message("\nBuilding network: ", nm ,"...\n")
  
  #Build network
  ptn_net<-build_string_network(
    corr_sig_ptn_map_net[[ nm ]],
    mode                = "multi",
    wait_for_layout_adj = TRUE,
    hide_singletons     = TRUE
    )
  
  #Set contrasts
  set_contrast(
    contrast        = nm,
    class_col       = paste0("sign_", nm ),
    width_col       = paste0("abs_corr_", nm ),
    expr_col        = paste0("mean_expr_protein_", nm ),
    filename        = file.path(
      paths_networks_path,paste0('paths-networks-ptn-net-', nm )),
    wait_for_legend = TRUE) 
  
  message("Network ", nm ," exported.\n")
}
#Open each .SVG file in adobe illustrator and adjust

### Export other Main Figure Panels ###
ggsave(file.path(paths_networks_path,"paths-networks-expr-ht.pdf"),
       ht_corr,
       height=8,width=12,units='in')
ggsave(file.path(paths_networks_path,"paths-networks-pfg-ht.pdf"),
       pfg_ht_plot,
       height=9,width=6,units='in')
ggsave(file.path(paths_networks_path,"paths-networks-react-ht.pdf"),
       react_ht_plot,
       height=9,width=6,units='in')
#Compile final Figure in adobe illustrator

