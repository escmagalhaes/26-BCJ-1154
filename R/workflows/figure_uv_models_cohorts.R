### CoxPH UV models of OS by MECOM expression levels from several cohorts ###

#Import input data
uv_mods_tab_list<-qs2::qs_read(file.path(analyses_export_path,"UV_models_cohorts.qs2"))

#Create folder and path to export results
uv_mods_cohorts_path<-file.path(data_export_path,'uv-mods-cohorts')
if (!dir.exists(uv_mods_cohorts_path)) { 
  dir.create(uv_mods_cohorts_path,recursive=T) 
  }

#Stack tables into one
uv_nms<-names(uv_mods_tab_list)
uv_nms_ignore<-c(
  uv_nms[grepl("UTMDACC"  ,uv_nms,ignore.case=TRUE)  & grepl("RPPA",uv_nms,ignore.case=TRUE)],
  uv_nms[grepl("BeatAML",uv_nms,ignore.case=TRUE)  & grepl("mRNA",uv_nms,ignore.case=TRUE)],
  uv_nms[grepl("TCGA"   ,uv_nms,ignore.case=TRUE)  & grepl("mRNA",uv_nms,ignore.case=TRUE)]
)
uv_mods_tab_list_cleaned<-uv_mods_tab_list[!uv_nms %in% uv_nms_ignore]
uv_tab_merge<-tbl_stack(uv_mods_tab_list_cleaned) |> 
  modify_header(label="**MECOM levels by cohort**") |> 
  modify_caption(
    c("**Univariate Cox Proportional Hazards models of Overall Survival (OS) 
    across different cohorts and MECOM expression layers**")) |>
  add_forest() |>
  text_transform(locations=cells_body(
    columns=c(conf.low),rows=var_label==nms[grepl("Jayavelu",nms,ignore.case=TRUE)]),
    fn=function(x) { return("0.00, Inf") })
#Adjust the CI that was changed in get_uv_table() to make the forest plot plottable

### Export Table as .HTML ###
gtsave(uv_tab_merge,file.path(uv_mods_cohorts_path,"uv-mods-cohorts.html"),inline_css=TRUE)
#Export as PDF and adjust in illustrator whatever is needed
#EXPORT AS PNG TOO!!!
