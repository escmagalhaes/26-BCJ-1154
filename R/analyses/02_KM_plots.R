### KM plots analyses ###

#Import input datatables
main_data <-qs2::qs_read(file.path(analyses_export_path,"processed_main_data.qs2"))
ref_sets  <-qs2::qs_read(file.path(analyses_export_path,"processed_reference_sets.qs2"))
subsets   <-qs2::qs_read(file.path(analyses_export_path,"processed_subsets.qs2"))

#Create local objects
df_full     <-main_data$main_df
df_normal   <-ref_sets$normal_data
df_full_rem <-main_data$df_rem

df_wt_mecomr<-subsets$all$Chr3_status$Chr3_status_WT_MECOMR
df_wt_chr3_abn<-subsets$all$Chr3_status$Chr3_status_WT_Chr3abn

df_below<-subsets$all$MECOM_4lvl_expr$MECOM_4lvl_expr_Below_Normal
df_lower<-subsets$all$MECOM_4lvl_expr$MECOM_4lvl_expr_Lower_Normal
df_upper<-subsets$all$MECOM_4lvl_expr$MECOM_4lvl_expr_Upper_Normal
df_above<-subsets$all$MECOM_4lvl_expr$MECOM_4lvl_expr_Above_Normal

df_below_rem<-subsets$relapse$MECOM_4lvl_expr$MECOM_4lvl_expr_Below_Normal
df_lower_rem<-subsets$relapse$MECOM_4lvl_expr$MECOM_4lvl_expr_Lower_Normal
df_upper_rem<-subsets$relapse$MECOM_4lvl_expr$MECOM_4lvl_expr_Upper_Normal
df_above_rem<-subsets$relapse$MECOM_4lvl_expr$MECOM_4lvl_expr_Above_Normal

## Fixed KM plots parameters ##
fixed_km_params<-list(
  pval_coord=c(3.5,1),legend_position=c(0.6,0.75),ncol_legend=2
  ,legend_mode="inside",xlim=c(0,11),break_x_by=1,table_by_time=TRUE
  ,show_pval_table=TRUE,risk_table=FALSE,pval_table_adj_method='BH'
  ,add_overall=TRUE,linetype_overall="solid",color_overall="black"
  ,ylab='Cumulative probability',xlab="Time (years)"
  ,line_types=c("solid","twodash","dashed","dotted","dotdash","longdash")
)

## Varying KM plot parameters ##
#Grouping var and covar plotting ##
df_labs<-c("","and Chr3 Altered","and Chr3-WT vs MECOM-R","and Chr3-WT vs Chr3-abn")
km_df<-list('os'=df_full
            ,'os_all_chr3_alt'=df_full
            ,'os_mecomr'=df_wt_mecomr
            ,'os_chr3_abn'=df_wt_chr3_abn
)
km_df_rem<-list('rd'=df_full_rem)
covar_list<-list(NULL,"chr3_status_simple","chr3_status","chr3_status")

#KM by 3 groups
km_3gp<-get_varying_km_params(
  df_list=km_df,labs=df_labs,prefix='by_3groups'
  ,title='Overall Survival by MECOM levels'
  ,time_var="surv_time",event_var="status"
  ,single_color=FALSE,grouping_var="MECOM_gp_n3"
  ,color_mode='list',covar_list=covar_list,n_panels=12)

#KM by 6 groups
km_6gp<-get_varying_km_params(
  df_list=km_df,labs=df_labs,prefix='by_6groups'
  ,title='Overall Survival by MECOM levels'
  ,time_var="surv_time",event_var="status"
  ,single_color=FALSE,grouping_var="MECOM_gp6"
  ,color_mode='list',covar_list=covar_list,n_panels=12)

#KM by 2 groups
km_2gp<-get_varying_km_params(
  df_list=km_df,labs=df_labs,prefix='by_2groups'
  ,title='Overall Survival by MECOM levels'
  ,time_var="surv_time",event_var="status"
  ,single_color=FALSE,grouping_var="MECOM_gp6_n2"
  ,color_mode='list',covar_list=covar_list,n_panels=12)

#OS by 4 groups
os_4gp<-get_varying_km_params(
  df_list=km_df['os'],labs="",prefix='by_4groups'
  ,title='Overall Survival by MECOM levels'
  ,time_var="surv_time",event_var="status"
  ,single_color=FALSE,grouping_var="MECOM_gp_n"
  ,color_mode='list',covar_list=NULL,n_panels=6)

#RD by 4 groups
rd_4gp<-get_varying_km_params(
  df_list=km_df_rem,labs="",prefix='by_4groups'
  ,title='Remission Duration by MECOM levels'
  ,time_var="rem_time",event_var="relapse"
  ,single_color=FALSE,grouping_var="MECOM_gp_n"
  ,color_mode='list',covar_list=NULL,n_panels=6)

## Covar only plots ##
covar_labs<-paste("Chr3-ALT types in",c("Below","Lower","Upper","Above"),"Normal MECOM")
covar_df<-list('below_normal_chr3_type'=df_below
               ,'lower_normal_chr3_type'=df_lower
               ,'upper_normal_chr3_type'=df_upper
               ,'above_normal_chr3_type'=df_above)
covar_df_rem<-list('below_normal_chr3_type'=df_below_rem
                   ,'lower_normal_chr3_type'=df_lower_rem
                   ,'upper_normal_chr3_type'=df_upper_rem
                   ,'above_normal_chr3_type'=df_above_rem)

#OS covar varying parameters
vary_covar_os<-get_varying_km_params(
  df_list=covar_df,labs=covar_labs,prefix='os'
  ,title="Overall Survival by"
  ,time_var="surv_time",event_var="status"
  ,single_color=TRUE,grouping_var="MECOM_gp_n"
  ,color_mode='list',covar="chr3_status",n_panels=6)

#RD covar varying parameters
vary_covar_rd<-get_varying_km_params(
  df_list=covar_df_rem,labs=covar_labs,prefix='rd'
  ,title="Remission Duration by"
  ,time_var="rem_time",event_var="relapse"
  ,single_color=TRUE,grouping_var="MECOM_gp_n"
  ,color_mode='list',covar="chr3_status",n_panels=6)

## Combine all varying parameters in a single table ##
all_varying_params<-dplyr::bind_rows(
  os_4gp,
  rd_4gp,
  km_3gp,
  km_6gp,
  km_2gp,
  vary_covar_os,
  vary_covar_rd
)

## One single call to rule it all (i.e. make all plots) ##
km_res<-purrr::pmap(
  all_varying_params %>% dplyr::select(-name),
  function(...) {
    args<-list(...)
    args<-purrr::map( args, function(x) {
      if (is.list(x) && length(x)==1) {
        x<-x[[1]]  #unwrap single-element list
      }
      if (length(x)==1 && is.na(x)) NULL else x
    })
    do.call( plot_km, c( args, fixed_km_params ))
  }) %>% purrr::set_names(all_varying_params$name)

#Create lists of KM plots
#OS only
km_res_os<-km_res[
  c("by_4groups_os",
    "os_below_normal_chr3_type", "os_lower_normal_chr3_type",
    "os_upper_normal_chr3_type", "os_above_normal_chr3_type")]

#RD only
km_res_rd<-km_res[
  c("by_4groups_rd",
    "rd_below_normal_chr3_type", "rd_lower_normal_chr3_type",
    "rd_upper_normal_chr3_type", "rd_above_normal_chr3_type")]

#By quantile
km_res_qt<-km_res[
  c("by_3groups_os",        "by_3groups_os_all_chr3_alt",
    "by_3groups_os_mecomr", "by_3groups_os_chr3_abn",
    "by_6groups_os",        "by_6groups_os_all_chr3_alt",
    "by_6groups_os_mecomr", "by_6groups_os_chr3_abn",
    "by_2groups_os",        "by_2groups_os_all_chr3_alt",
    "by_2groups_os_mecomr", "by_2groups_os_chr3_abn")]

### Save analyses ###
qs2::qs_save(km_res_os,file.path(analyses_export_path,"km_results_os.qs2"))
qs2::qs_save(km_res_rd,file.path(analyses_export_path,"km_results_rd.qs2"))
qs2::qs_save(km_res_qt,file.path(analyses_export_path,"km_results_qt.qs2"))



