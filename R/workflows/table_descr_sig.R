### Significant Descriptive Tables ###

#Import input data
sig_tabs_full<-qs2::qs_read(file.path(analyses_export_path,"descr_sig_table.qs2"))

#Create local objects
sig_tabs  <-sig_tabs_full$sig_tabs
full_df   <-as.data.frame(sig_tabs_full$df_all)
color_map <-get_mycolors(n_colors=nlevels(full_df$group),mode='list')
names(color_map)<-levels(full_df$group)

#Create folder and path to export results
table_descr_sig_path<-file.path(data_export_path,'descr-tab-sig')
if (!dir.exists(table_descr_sig_path)) { 
  dir.create(table_descr_sig_path,recursive=T) 
}

#Extract sig tables
sig_tab_list<-sig_tabs[grepl("^sig_table_",names(sig_tabs))] 

#Descriptive table with significant variables only #
comb_sig_tabs<-tbl_merge(sig_tab_list,quiet=TRUE) |>
  bold_labels()|>
  modify_header(label~'**Variable**') |>
  modify_caption(
    "**Significant Demographic, Clinical and Molecular Characteristics
  of AML patients stratified by MECOM levels and grouped by Chr3 alteration status**") |>
  modify_spanning_header(
    paste0(c('stat_1','stat_2','p.value'),'_1')~'**All cases**',
    paste0(c('stat_1','stat_2','p.value'),'_2')~'**Chr3-WT**',
    paste0(c('stat_1','stat_2','p.value'),'_3')~'**All Chr3-ALT.**',
    paste0(c('stat_1','stat_2','p.value'),'_4')~'**MECOM-R**',
    paste0(c('stat_1','stat_2','p.value'),'_5')~'**Chr3-abn**'  
  )|> as_gt()|> 
  gt::tab_style(
    style=list(cell_text(color=color_map['Low'],weight="bold")),
    locations=cells_column_labels(paste0('stat_1_',seq_along(sig_tab_list)))) |>
  gt::tab_style(
    style=list(cell_text(color=color_map['High'],weight="bold")),
    locations=cells_column_labels(paste0('stat_2_',seq_along(sig_tab_list))))

### Export Table as .XLSx ###
gtsave(comb_sig_tabs,file.path(table_descr_sig_path,"descr-tab-sig.html"),inline_css=TRUE)
#Open the file with excel and mannually save as .xlsx
#EXPORT AS PNG TOO!!!

