### CoxPH UV models of OS and RD by MECOM expression levels (RPPA) ###

#Import input data
uv_mod_os    <-qs2::qs_read(file.path(analyses_export_path,"UV_model_OS_table.qs2"))
uv_mod_rd    <-qs2::qs_read(file.path(analyses_export_path,"UV_model_RD_table.qs2"))
mv_mods_tabs <-qs2::qs_read(file.path(analyses_export_path,"MV_model_tables.qs2"))

#Create local objects
uv_mv_df  <-mv_mods_tabs$uv_mv_df
color_map <-get_mycolors(n_colors=nlevels(uv_mv_df$group),mode='list')
names(color_map)<-levels(uv_mv_df$group)

#Create folder and path to export results
table_uv_models_path<-file.path(data_export_path,'uv-models')
if (!dir.exists(table_uv_models_path)) { 
  dir.create(table_uv_models_path,recursive=T) 
}

#Helper function to identify rows to change design
vseq<-function(n_ini,n_var){ seq(n_ini,((n_var*3)+n_ini),3) } 
#n_ini=initial row; n_var=number of variables that are intercalated between first and last row
rows_sv  <-c(2,6,vseq(13,3),vseq(26,3),vseq(38,28))
rows_rem <-c(2,6,vseq(13,3),vseq(26,4),vseq(47,3)
            ,vseq(62,1),vseq(77,2),92,vseq(98,1),107,vseq(113,1))

# Combine OS and RD UV tables and adjust design #
uv_comb_tabs<-tbl_merge(
  list(uv_mod_os,uv_mod_rd)
  ,quiet=TRUE,tab_spanner=c("**Univariate OS**","**Univariate RD**")) |>
  bold_labels()|>
  modify_header(label~"**Variable**")|>
  modify_caption("**Univariate Cox Proportional Hazards models of 
  Overall Survival (OS) and Remission Duration (RD) in Chr3 WT patients**") |>
  as_gt() |>
  gt::sub_missing(columns=c(conf.low_1),rows=rows_sv,missing_text="_________") |> 
  gt::sub_missing(columns=c(p.value_1),rows=rows_sv,missing_text="reference") |> 
  gt::tab_style(style=cell_text(v_align="top")
                ,locations=cells_body(columns=c(conf.low_1),rows=rows_sv)) |>
  gt::sub_missing(columns=c(conf.low_2),rows=rows_rem,missing_text="_________") |> 
  gt::sub_missing(columns=c(p.value_2),rows=rows_rem,missing_text="reference") |> 
  gt::tab_style(style=cell_text(v_align="top")
                ,locations=cells_body(columns=c(conf.low_2),rows=rows_rem)) |>
  gt::tab_style(style=list(cell_text(color=color_map['Low'],weight="bold"))
                ,locations=cells_body(columns=label,rows=2)) |>
  gt::tab_style(style=list(cell_text(color=color_map['High'],weight="bold"))
                ,locations=cells_body(columns=label,rows=3))

### Export Table as .XLSx ###
gtsave(uv_comb_tabs,file.path(table_uv_models_path,"uv-models.html"),inline_css=TRUE)
#Open the file with excel and mannually save as .xlsx
#EXPORT AS PNG TOO!!!


