### Export multi-dataset ###

#Import input data
multi_data<-qs2::qs_read(file.path(analyses_export_path,"multi_cohort_dataset.qs2"))

#Create folder and path to export results
multi_cohort_dataset_path<-file.path(data_export_path,'multi-cohort-dataset')
if (!dir.exists(multi_cohort_dataset_path)) { 
  dir.create(multi_cohort_dataset_path,recursive=T) 
}

#Adjust dataset to export
multi_data_export<-multi_data |> dplyr::select(
  UPI=upi,Cohort=cohort,`MECOM expression`=MECOM_zscore,`Chr3 Status`=chr3_status,
  `Chr3 alteration`=chr3_status_cat,`MECOM-R variant`=MECOM_R_variant,
  Cytogenetics=karyotype) |> distinct(UPI,.keep_all = TRUE)

## Export as XLSX ##
export_tabs(
  df_list   = multi_data_export,
  filepath  = multi_cohort_dataset_path,
  filename  = 'multi-cohort-dataset',
  tab_title = "MECOM expression, Chr3 status and Cytogenetics by cohort")
