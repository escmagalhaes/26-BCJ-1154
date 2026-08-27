### Extended Descriptive Tables ###

#Import input data
descr_tabs_full<-qs2::qs_read(file.path(analyses_export_path,"descr_ext_table.qs2"))

#Create local objects
descr_tabs <-descr_tabs_full$descr_tabs
full_df    <-as.data.frame(descr_tabs_full$df_all)
color_map  <-get_mycolors(n_colors=nlevels(full_df$group),mode='list')
names(color_map)<-levels(full_df$group)

#Create folder and path to export results
table_descr_ext_path<-file.path(data_export_path,'descr-tab-ext')
if (!dir.exists(table_descr_ext_path)) { 
  dir.create(table_descr_ext_path,recursive=T) 
}

#Extract relevant tables
descr_tab_list<-descr_tabs[!names(descr_tabs) %in% c("high","low")]

#Descriptive table with all variables of CH3 status by MECOM High vs Low levels#
comb_descr_tabs<-tbl_merge(descr_tab_list,quiet=TRUE) |>
  bold_labels()|> modify_header(label~'**Variable**')|>
  modify_column_hide(columns=paste0('stat_0_',seq_along(descr_tab_list))) |>
  modify_caption(
  "**Expanded Demographic, Clinical and Molecular Characteristics
  of AML patients stratified by MECOM levels and grouped by Chr3 alteration status.**") |>
  modify_spanning_header(
    paste0(c('stat_1','stat_2','p.value'),'_1')~'**All cases**',
    paste0(c('stat_1','stat_2','p.value'),'_2')~'**Chr3-WT**',
    paste0(c('stat_1','stat_2','p.value'),'_3')~'**All Chr3-ALT.**',
    paste0(c('stat_1','stat_2','p.value'),'_4')~'**MECOM-R**',
    paste0(c('stat_1','stat_2','p.value'),'_5')~'**Chr3-abn**'
  )|> as_gt()|> 
  gt::tab_style(
    style=list(cell_text(color=color_map['Low'],weight ="bold")),
    locations=cells_column_labels(paste0('stat_1_',seq_along(descr_tab_list))))|>
  gt::tab_style(
    style=list(cell_text(color=color_map['High'],weight ="bold")),
    locations=cells_column_labels(paste0('stat_2_',seq_along(descr_tab_list))))

### Export Table as .XLSx ###
gtsave(
  comb_descr_tabs,
  file.path(table_descr_ext_path,"descr-tab-ext.html"),
  inline_css=TRUE
  )
#Open the file with excel and mannually save as .xlsx
#EXPORT AS PNG TOO!!!

