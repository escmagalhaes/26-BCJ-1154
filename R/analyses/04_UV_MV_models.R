### UV/MV CoxPH models ###

#Import input datatables
ref_sets  <-qs2::qs_read(file.path(analyses_export_path,"processed_reference_sets.qs2"))
subsets   <-qs2::qs_read(file.path(analyses_export_path,"processed_subsets.qs2"))

#Create local objects
df_wt     <-subsets$all$Chr3_status$Chr3_status_WT
df_wt_rem <-subsets$relapse$Chr3_status$Chr3_status_WT
meta_tab  <-ref_sets$meta_table
mod_vars  <-ref_sets$meta_table$variable

#Create variable parameters list for UV/MV analyses in CH3 WT patients
uv_mv_params_list<-list(
  'OS' = list(
    df        = df_wt,
    time_var  = "surv_time",
    event_var = "status"
  ),
  'RD' = list(
    df        = df_wt_rem,
    time_var  = "rem_time",
    event_var = "relapse"
  )
)

#Sensitivity analyses of MV model strategies to determine best model 
fixed_params_sens<-list(
  grouping_variable = "MECOM_gp_n2",
  meta_table        = meta_tab,
  vars              = mod_vars,
  group_var_label   = "MECOM levels",
  strategies        = c("strict","two_phase","cie"),
  priority_levels   = c("core","optional"),
  label_vars        = list(aml_gp="Secondary AML"),
  relevel_vars      = list(aml_gp="primary"),
  dichotomous_vars  = list(aml_gp="secondary"),
  block_vars        = TRUE,
  cie_final_cleanup = TRUE,
  mv_independent    = TRUE,
  run_bootstrap     = TRUE,
  parallel          = TRUE,
  verbose           = TRUE
  )

#Get all sensitivity analyses in a single call
sens_res<-lapply(uv_mv_params_list,function( p ) {
  args<-c(fixed_params_sens, p )
  do.call(run_mv_sensitivity, args )
})

#Extract sensitivity results
sens_os<-sens_res$OS$comparison
sens_rd<-sens_res$RD$comparison

#Compare Results to define best model for OS and RD
strategy_names<-c(
  'Forward Stepwise Selection model (FSS)',
  'FSS with Backward variable Prunning (FSSBP)',
  'Change in Estimate Model (CIE)')
sens_all_tab<-bind_rows(
  sens_os |> 
    dplyr::mutate(Model = "Overall Survival"  ,Strategy = strategy_names) |>
    dplyr::select(all_of('Model'),everything()),
  sens_rd |> 
    dplyr::mutate(Model = "Remission Duration",Strategy = strategy_names) |>
    dplyr::select(all_of('Model'),everything())
)

#Best Overall Survival and Remission Duration models: strict for both
uv_tab_os<-sens_res$OS$results$strict$uv_table_full
uv_tab_rd<-sens_res$RD$results$strict$uv_table_full
mv_tab_os<-sens_res$OS$results$strict$mv_table
mv_tab_rd<-sens_res$RD$results$strict$mv_table
uv_mv_df <-sens_res$OS$results$strict$processed_df

### Save analyses ### 
#MV tables and sensitivity analyses
qs2::qs_save(list(
  uv_mv_df          = uv_mv_df,
  sensitivity_table = sens_all_tab,
  mv_table_os       = mv_tab_os,
  mv_table_rd       = mv_tab_rd
),file.path(analyses_export_path,"MV_model_tables.qs2"))

#UV table OS
qs2::qs_save(uv_tab_os,file.path(analyses_export_path,"UV_model_OS_table.qs2"))

#UV table RD
qs2::qs_save(uv_tab_rd,file.path(analyses_export_path,"UV_model_RD_table.qs2"))
