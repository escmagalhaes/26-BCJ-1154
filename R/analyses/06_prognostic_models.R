### Prognostic models analyses ###

#Import inputs
main_data <-qs2::qs_read(file.path(analyses_export_path,"processed_main_data.qs2"))
ref_sets  <-qs2::qs_read(file.path(analyses_export_path,"processed_reference_sets.qs2"))
subsets   <-qs2::qs_read(file.path(analyses_export_path,"processed_subsets.qs2"))

#Create local objects
df_full  <-main_data$main_df
meta_tab <-ref_sets$meta_table
df_prog  <-subsets$all$therapy$therapy_IC_VH

bxp_colors<-get_mycolors(n_colors=nlevels(df_full$MECOM_gp_n2_chr3_status),mode='vector')
bxp_colors<-setNames(bxp_colors,levels(df_full$MECOM_gp_n2_chr3_status))
bxp_colors[c("Chr3-WT Low","Chr3-WT High")]<-bxp_colors[c("Chr3-WT High","Chr3-WT Low")]
bxp_colors[c("Chr3-abn Low","Chr3-abn High")]<-bxp_colors[c("Chr3-abn High","Chr3-abn Low")]
#Adjust color order to match the levels appropriately

#Create varying parameters vector for Boxplots
ptn_bxp<-c(GATA2="GATA2",CEBPA="CEBPA")

#Create fixed parameters list for Boxplots
fixed_bxp_prog_params<-list(
  df              = df_full,
  discrete_vars   = 'MECOM_gp_n2_chr3_status',
  disc_var_labels = 'Chr3 alterations',
  title_suffix    = 'and MECOM expression',
  colors          = bxp_colors,
  padj            = "BH"
)

#Create all boxplots with a single call
bxp_prog_res<-lapply(ptn_bxp, function(x) {
  p    <-list( continuous_vars = x, ylab = paste0( x ,' levels'))
  args <-c( fixed_bxp_prog_params , p )
  do.call( run_bxplots , args )
})

#Cluster patients according to MECOM and CEBPA prognostic labels (optimal clusters=3)
prog_cl<-progeny_clusters(
  df             = df_prog,
  n_clusters     = 3,
  cluster_prefix = 'C',
  plot_gap       = FALSE,
  proteins       = c('MECOM_prog_num','CEBPA_prog_num')
)

#Generate KM plots
## Fixed KM plots parameters ##
fixed_km_prog_params<-list(
  grouping_var          = "cluster",
  pval_coord            = c(3.5,1),
  legend_position       = c(0.6,0.75),
  ncol_legend           = 2,
  legend_mode           = "inside",
  xlim                  = c(0,11),
  break_x_by            = 1,
  table_by_time         = TRUE
  ,show_pval_table      = TRUE,
  risk_table            = FALSE,
  pval_table_adj_method = 'BH',
  add_overall           = TRUE,
  linetype_overall      = "solid",
  color_overall         = "black",
  ylab                  = 'Cumulative probability',
  xlab                  = "Time (years)",
  n_panels              = 2,
  line_types            = c("solid","twodash","dashed","dotted","dotdash","longdash")
)

## Varying KM plot parameters ##
km_prog_df_list<-list('OS'=prog_cl$df_full,'RD'=prog_cl$df_rem)
varying_km_tab<-tibble::tibble(
  name      = paste('prog',names(km_prog_df_list),sep="_"),
  data      = unname(km_prog_df_list),
  label     = 'IC vs VH (N=469)',
  time_var  = c('surv_time','rem_time'),
  event_var = c('status','relapse'),
  title     = c('Overall Survival','Remission Duration'),
) |> dplyr::mutate(title=trimws(paste(title,label,sep=' '))) |>
  dplyr::select(-label)

## One call to rule it all (i.e. make all plots) ##
km_prog_res<-purrr::pmap(
  varying_km_tab |> dplyr::select(-name),
  function(...) {
    args<-list(...)
    args<-purrr::map( args, function(x) {
      if (is.list(x) && length(x)==1) {
        x<-x[[1]]  #unwrap single-element list
      }
      if (length(x)==1 && is.na(x)) NULL else x
    })
    do.call( plot_km, c( args, fixed_km_prog_params ))
  }) |> purrr::set_names(varying_km_tab$name)

#Sensitivity analyses of MV model strategies for OS to determine best model 
sens_sv_prog<-run_mv_sensitivity(
  df                = prog_cl$df_full,
  grouping_variable = "cluster",
  meta_table        = meta_tab,
  vars              = meta_tab$variable,
  time_var          = "surv_time",
  event_var         = "status",
  group_var_label   = "Prognostic Groups",
  strategies        = c("strict","two_phase","cie"),
  priority_levels   = c("core","optional"),
  label_vars        = list(aml_gp="Secondary AML"),
  relevel_vars      = list(aml_gp="primary"),
  dichotomous_vars  = list(aml_gp="secondary"),
  all_levels_sig    = TRUE,
  block_vars        = TRUE,
  cie_final_cleanup = TRUE,
  mv_independent    = TRUE,
  run_bootstrap     = TRUE,
  parallel          = TRUE,
  verbose           = TRUE
  )

### Save analyses ### 
qs2::qs_save(list(
  dataset_clusters  = prog_cl$df_full,
  meta_table        = meta_tab,
  boxplots          = bxp_prog_res, 
  km_plots          = km_prog_res,
  mv_table          = sens_sv_prog$results$strict$mv_table,
  sensitivity_table = sens_sv_prog$comparison
),file.path(analyses_export_path,"prognostic_models.qs2"))

