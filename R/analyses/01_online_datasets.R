### Analyses of data from public repositories ###

#Import input datatables
cohort_data<-qs2::qs_read(file.path(analyses_export_path,"processed_cross_study_data.qs2"))

#Create local object
multi_data<-cohort_data |> dplyr::arrange(cohort,type) |> tidyr::unite(
  "label",cohort,type,sep="; ",remove = FALSE) |> dplyr::mutate(
    label = factor(label)) |> dplyr::relocate(label,.after=everything())

#Create list with datasets
label_levels<-levels(multi_data$label)
dataset_list<-multi_data |> 
  dplyr::group_split(label,.keep=TRUE) |>  
  purrr::set_names(label_levels)

#Customize dataset_list order
nms<-names(dataset_list)
nms_order<-c(
  nms[grepl("UTMDACC",nms,ignore.case=TRUE)  & grepl("RPPA",nms,ignore.case=TRUE)],
  nms[grepl("UTMDACC",nms,ignore.case=TRUE)  & grepl("mRNA",nms,ignore.case=TRUE)],
  nms[grepl("BeatAML",nms,ignore.case=TRUE)  & grepl("mRNA",nms,ignore.case=TRUE)],
  nms[grepl("TCGA",nms,ignore.case=TRUE)     & grepl("mRNA",nms,ignore.case=TRUE)],
  nms[grepl("Kramer",nms,ignore.case=TRUE)   & grepl("mRNA",nms,ignore.case=TRUE)],
  nms[grepl("Kramer",nms,ignore.case=TRUE)   & grepl("MS",nms,ignore.case=TRUE)],
  nms[grepl("Jayavelu",nms,ignore.case=TRUE) & grepl("MS",nms,ignore.case=TRUE)],
  nms[grepl("BeatAML",nms,ignore.case=TRUE)  & grepl("MS",nms,ignore.case=TRUE)],
  nms[grepl("MDACC",nms,ignore.case=TRUE)    & grepl("MS",nms,ignore.case=TRUE)]
)
dataset_list_reordered<-dataset_list[nms_order]

#Remove datasets with all NA for MECOM_zscore and merge back for export later
clean_list <-dataset_list_reordered |> purrr::discard(~all(is.na(.x$MECOM_zscore)))
clean_df   <-clean_list |> purrr::map(as_tibble) |> dplyr::bind_rows()

#Make a tibble for varying parameters
uv_tab<-tibble::tibble(
  data   = clean_list,
  labels = as.list(names(clean_list))
)

#Get uv_models
uv_tab_list<-purrr::pmap(uv_tab,function(data,labels) {
  get_uv_table(
    data         = data, 
    features     = "MECOM_zscore",
    labels       = labels,
    count_vars   = c("MECOM_R","chr3_abn"),
    count_labels = c("MECOM-R","Chr3-abn"),
    sum_value    = "yes",
    time_var     = "surv_time",
    event_var    = "status",
    adjust_CI    = TRUE,
    max_CI       = 2.2,
    min_CI       = 0.55
  )
})
names(uv_tab_list)<-unlist(uv_tab$labels)

#Extract relevant models for meta-analysis
nms_use<-c(
  nms[grepl("UTMDACC"  ,nms,ignore.case=TRUE)  & grepl("RPPA",nms,ignore.case=TRUE)],
  nms[grepl("BeatAML",nms,ignore.case=TRUE)  & grepl("mRNA",nms,ignore.case=TRUE)],
  nms[grepl("TCGA"   ,nms,ignore.case=TRUE)  & grepl("mRNA",nms,ignore.case=TRUE)]
)

meta_dat<-tibble(
  cohort     = c(nms_use),
  method     = c("RPPA","mRNA","mRNA"),
  HR         = c(1.22, 1.14, 1.12),
  lowerCI    = c(1.13, 1.04, 0.93),
  upperCI    = c(1.33, 1.26, 1.34),
  Total      = c(691, 591, 169),
  Missing    = c( 0  , 0 ,  0 ),
  `MECOM-R`  = c( 20 , 10,  2 ),
  `Chr3-abn` = c( 54 , 26,  6 ),
  Model      = c(691, 591, 169),
  Event      = c(528, 355, 108),
  `p-value`  = c("<0.001", "0.008", "0.25")
)

#Calculate the SE for 95% CI and Sampling Variance
meta_dat<-meta_dat |> dplyr::mutate(
  yi  = log(HR),
  sei = (log(upperCI) - log(lowerCI)) / (2 * 1.96),
  vi  = sei^2
)

#Run fixed-effects meta-regression model
meta_res<-metafor::rma(
  yi     = yi,
  sei    = sei,
  data   = meta_dat,
  method = "REML"
)

##KM plots##
#Define datasets and grouping variable#
beataml_dataset <-clean_list[["Bottomly, D. et al. (BeatAML); mRNA"]]
tcga_dataset    <-clean_list[["TCGA et al. (LAML); mRNA"]]
km_datasets<-list(beataml_dataset,tcga_dataset)
names(km_datasets)<-c("Bottomly, D. et al. (BeatAML); mRNA","TCGA et al. (LAML); mRNA")

#Determine best RNA quantile split
beataml_rna_qt<-determine_optimal_quantile(
  data      = beataml_dataset,
  feature   = "MECOM_zscore",
  time_var  = "surv_time", 
  event_var = "status",
  q_range   = 2:6 #preserve at least 15% minimal proportion for the smaller group
)

tcga_rna_qt<-determine_optimal_quantile(
  data      = tcga_dataset,
  feature   = "MECOM_zscore",
  time_var  = "surv_time", 
  event_var = "status",
  q_range   = 2:6 #preserve at least 15% minimal proportion for the smaller group
)

#KM plots
beataml_rna_km_res<-plot_km(
  data                  = beataml_rna_qt,
  time_var              = "surv_time",
  event_var             = "status",
  grouping_var          = 'MECOM_zscore_group',
  title                 = paste0(
    "Overall Survival by MECOM mRNA levels\n(BeatAML cohort; N=", nrow(beataml_rna_qt),")"),
  n_panels              = 6,
  colors                = c('blue3','red2'),
  pval_coord            = c(3.5,1),
  legend_position       = c(0.6,0.75),
  ncol_legend           = 2,
  legend_mode           = "inside",
  xlim                  = c(0,9),
  break_x_by            = 1,
  table_by_time         = TRUE,
  show_pval_table       = TRUE,
  risk_table            = FALSE,
  pval_table_adj_method = 'BH',
  add_overall           = TRUE,
  linetype_overall      = "solid",
  color_overall         = "black",
  ylab                  = 'Cumulative probability',
  xlab                  = "Time (years)",
  line_types            = c("solid","twodash","dashed","dotted","dotdash","longdash")
  
)

beataml_rna_km_res_covar<-plot_km(
  data                  = beataml_rna_qt,
  time_var              = "surv_time",
  event_var             = "status",
  grouping_var          = 'MECOM_zscore_group',
  covar                 = "chr3_status",
  title                 = paste0(
    "Overall Survival by MECOM mRNA levels and Chr3 Status\n(BeatAML cohort; N=", nrow(beataml_rna_qt),")"),
  n_panels              = 6,
  colors                = c('blue3','red2'),
  pval_coord            = c(3.5,1),
  legend_position       = c(0.6,0.75),
  ncol_legend           = 2,
  legend_mode           = "inside",
  xlim                  = c(0,9),
  break_x_by            = 1,
  table_by_time         = TRUE,
  show_pval_table       = TRUE,
  risk_table            = FALSE,
  pval_table_adj_method = 'BH',
  add_overall           = TRUE,
  linetype_overall      = "solid",
  color_overall         = "black",
  ylab                  = 'Cumulative probability',
  xlab                  = "Time (years)",
  line_types            = c("solid","twodash","dashed","dotted","dotdash","longdash")
  
)

tcga_rna_km_res<-plot_km(
  data                  = tcga_rna_qt,
  time_var              = "surv_time",
  event_var             = "status",
  grouping_var          = 'MECOM_zscore_group',
  title                 = paste0(
    "Overall Survival by MECOM mRNA levels\n(TCGA AML cohort; N=",nrow(tcga_rna_qt),")"),
  n_panels              = 6,
  colors                = c('blue3','red2'),
  pval_coord            = c(3.5,1),
  legend_position       = c(0.6,0.75),
  ncol_legend           = 2,
  legend_mode           = "inside",
  xlim                  = c(0,9),
  break_x_by            = 1,
  table_by_time         = TRUE,
  show_pval_table       = TRUE,
  risk_table            = FALSE,
  pval_table_adj_method = 'BH',
  add_overall           = TRUE,
  linetype_overall      = "solid",
  color_overall         = "black",
  ylab                  = 'Cumulative probability',
  xlab                  = "Time (years)",
  line_types            = c("solid","twodash","dashed","dotted","dotdash","longdash")
  
)

tcga_rna_km_res_covar<-plot_km(
  data                  = tcga_rna_qt,
  time_var              = "surv_time",
  event_var             = "status",
  grouping_var          = 'MECOM_zscore_group',
  covar                 = "chr3_status",
  title                 = paste0(
    "Overall Survival by MECOM mRNA levels and Chr3 Status\n(TCGA AML cohort; N=",nrow(tcga_rna_qt),")"),
  n_panels              = 6,
  colors                = c('blue3','red2'),
  pval_coord            = c(3.5,1),
  legend_position       = c(0.6,0.75),
  ncol_legend           = 2,
  legend_mode           = "inside",
  xlim                  = c(0,9),
  break_x_by            = 1,
  table_by_time         = TRUE,
  show_pval_table       = TRUE,
  risk_table            = FALSE,
  pval_table_adj_method = 'BH',
  add_overall           = TRUE,
  linetype_overall      = "solid",
  color_overall         = "black",
  ylab                  = 'Cumulative probability',
  xlab                  = "Time (years)",
  line_types            = c("solid","twodash","dashed","dotted","dotdash","longdash")
  
)

#Create plot list
kmplt_list<-list(
  "beataml_os_2lvl"      = beataml_rna_km_res,
  "tcga_os_2lvl"         = tcga_rna_km_res,
  "beataml_os_2lvl_chr3" = beataml_rna_km_res_covar,
  "tcga_os_2lvl_chr3"    = tcga_rna_km_res_covar
)

### Save analyses ###
qs2::qs_save(list(
  meta_df  = meta_dat,
  meta_mod = meta_res
),file.path(analyses_export_path,"meta_analysis_cohorts.qs2"))

qs2::qs_save(kmplt_list,file.path(analyses_export_path,"km_plots_cohorts.qs2"))
qs2::qs_save(uv_tab_list,file.path(analyses_export_path,"UV_models_cohorts.qs2"))
qs2::qs_save(clean_df,file.path(analyses_export_path,"multi_cohort_dataset.qs2"))

