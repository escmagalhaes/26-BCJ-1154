### Descriptive tables analyses ###

#Import input datatables
ref_sets  <-qs2::qs_read(file.path(analyses_export_path,"processed_reference_sets.qs2"))
subsets   <-qs2::qs_read(file.path(analyses_export_path,"processed_subsets.qs2"))

#Create local objects
df_all        <-subsets$all$Chr3_status$all
df_wt         <-subsets$all$Chr3_status$Chr3_status_WT
df_chr3_alt   <-subsets$all$Chr3_status$Chr3_status_altered
df_mecomr     <-subsets$all$Chr3_status$Chr3_status_MECOMR
df_chr3_abn   <-subsets$all$Chr3_status$Chr3_status_Chr3abn
df_high       <-subsets$all$MECOM_2lvl_expr$MECOM_2lvl_expr_High
df_low        <-subsets$all$MECOM_2lvl_expr$MECOM_2lvl_expr_Low
meta_tab      <-ref_sets$meta_table
decr_tab_vars <-ref_sets$meta_table$variable

#Create varying parameters list for descriptive tables
varying_params_descr_tab<-list(
  'all' = list(
    df                = df_all,
    grouping_variable = "MECOM_gp_n2"
    ),
  'wt' = list(
    df                = df_wt,
    grouping_variable = "MECOM_gp_n2"
    ),
  'all_chr3_alt' = list(
    df                = df_chr3_alt,
    grouping_variable = "MECOM_gp_n2"
    ),
  'mecomr' = list(
    df                = df_mecomr,
    grouping_variable = "MECOM_gp_n2"
    ),
  'chr3_abn' = list(
    df                = df_chr3_abn,
    grouping_variable = "MECOM_gp_n2"
    ),
  'high' = list(
    df                = df_high,
    grouping_variable = "chr3_status"
    ),
  'low' = list(
    df                = df_low,
    grouping_variable = "chr3_status"
    )
)

#Fixed parameters
fixed_params_descr_tab<-list(
  meta_table       = meta_tab,
  vars             = decr_tab_vars,
  descr_tab        = TRUE,
  uv_mod           = FALSE,
  mv_mod           = FALSE,
  label_vars       = list(aml_gp="Secondary AML"),
  relevel_vars     = list(aml_gp="primary"),
  dichotomous_vars = list(aml_gp="secondary"),
  verbose          = FALSE
)

#Create all tables in a single call
tab_res<-lapply( varying_params_descr_tab ,function( p ) {
  args<-c( fixed_params_descr_tab , p )
  do.call( run_tab_uv_mv_pipeline , args )
})

#Extract relevant objects 
df_all         <-tab_res$all$processed_df
descr_tab_list <-lapply(tab_res, function( x ) x$descr_table)

#Create tables with significant variables only to merge appropriately later
signif_tabs<-get_comb_sig_tables(
  table_list       = tab_res[!names(tab_res) %in% c("high","low")],
  exclude_sig      = c("race","cyto_unfav",'hgb','vtx','mut_ras'),
  label_vars       = list(aml_gp="Secondary AML"),
  dichotomous_vars = list(aml_gp="secondary")
)

### Save analyses ### 
#Descriptive table with significant only
qs2::qs_save(list(
  sig_tabs = signif_tabs,
  df_all   = list(df_all)
  ),file.path(analyses_export_path,"descr_sig_table.qs2"))

#Extended table 
qs2::qs_save(list(
  descr_tabs = descr_tab_list,
  df_all     = list(df_all)
),file.path(analyses_export_path,"descr_ext_table.qs2"))
