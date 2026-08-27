### CoxPH MV models of OS and RD by MECOM expression levels (RPPA) ###

#Import input data
mv_mods_tabs<-qs2::qs_read(file.path(analyses_export_path,"MV_model_tables.qs2"))

#Create local objects
uv_mv_df  <-mv_mods_tabs$uv_mv_df
os_mv_tab <-mv_mods_tabs$mv_table_os
rd_mv_tab <-mv_mods_tabs$mv_table_rd
color_map <-get_mycolors(n_colors=nlevels(uv_mv_df$group),mode='list')
names(color_map)<-levels(uv_mv_df$group)

#Create folder and path to export results
table_mv_models_path<-file.path(data_export_path,'mv-models')
if (!dir.exists(table_mv_models_path)) { 
  dir.create(table_mv_models_path,recursive=T) 
}

# Combine OS and RD MV tables and adjust design #
mv_comb_tabs<-tbl_stack(list(
  os_mv_tab |> modify_header(label~"**Variable**") |> 
  modify_caption(
    c("**Multivariate Cox Proportional Hazards models of Overall Survival (OS) 
      and Remission Duration (RD) in Chr3 WT patients**",
      "**Overall Survival**")
  ),
  rd_mv_tab |> modify_header(label~"**Variable**") |> 
    add_variable_group_header(header = "Remission Duration",variables = group) |> 
    sumExtras::add_group_styling(format="bold") 
  )) |>  add_forest() |>
  gt::sub_missing(columns=conf.low,rows=c(2,12),missing_text="_________") |> 
  gt::sub_missing(columns=p.value,rows=c(2,12),missing_text="reference") |> 
  gt::tab_style(style=cell_text(v_align="top")
                ,locations=cells_body(columns=conf.low,rows=c(2,12))) |> 
  gt::tab_style(style=list(cell_text(color=color_map['Low'],weight="bold"))
                ,locations=cells_body(columns=label,rows=c(2,12))) |>
  gt::tab_style(style=list(cell_text(color=color_map['High'],weight="bold"))
                ,locations=cells_body(columns=label,rows=c(3,13))) 

### Export Table as .HTML ###
gtsave(mv_comb_tabs,file.path(table_mv_models_path,"mv-models.pdf"),vwidth=1754,vheight=1240,expand=0)
#Export as PDF and adjust in illustrator whatever is needed
#EXPORT AS PNG TOO!!!

