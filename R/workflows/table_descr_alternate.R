### Extended Descriptive Tables ###

#Import input data
descr_tabs_full<-qs2::qs_read(file.path(analyses_export_path,"descr_ext_table.qs2"))

#Create local objects
descr_tabs <-descr_tabs_full$descr_tabs
full_df    <-as.data.frame(descr_tabs_full$df_all)
color_map  <-get_mycolors(n_colors=nlevels(full_df$group),mode='list')
names(color_map)<-levels(full_df$group)

#Create folder and path to export results
table_descr_alternate_path<-file.path(data_export_path,'descr-tab-alternate')
if (!dir.exists(table_descr_alternate_path)) { 
  dir.create(table_descr_alternate_path,recursive=T) 
}

#Extract relevant tables
descr_tab_list_alt<-descr_tabs[names(descr_tabs) %in% c("high","low")]

# Descriptive table with all variables of MECOM High vs Low levels by CH3 status #
comb_descr_tabs_alt<-tbl_merge(descr_tab_list_alt,quiet=TRUE) %>%
  bold_labels()%>%
  modify_header(label~'**Variable**')%>%
  modify_column_hide(columns=paste0('stat_0_',seq_along(descr_tab_list_alt))) %>%
  modify_caption(
    "Expanded Demographic, Clinical and Molecular Characteristics 
   of AML patients stratified by Chr3 alteration status and grouped by MECOM levels.") %>%
  modify_spanning_header(
    'p.value_1'~'.','p.value_2'~'.'
    ,c(paste0('stat_',c(1:nlevels(full_df$chr3_status)),'_1'))~'**High MECOM**'
    ,c(paste0('stat_',c(1:nlevels(full_df$chr3_status)),'_2'))~'**Low MECOM**') %>% 
  as_gt() %>%
  gt::tab_style(
    style=cell_text(color=color_map['Low'],weight="bold"),
    location=cells_column_spanners(spanners=contains("stat_1_2"))) %>%
  gt::tab_style(
    style=cell_text(color=color_map['High'],weight="bold"),
    location=cells_column_spanners(spanners=contains("stat_1_1"))) 

### Export Table as .XLSx ###
gtsave(
  comb_descr_tabs_alt,
  file.path(table_descr_alternate_path,"descr-tab-alternate.html"),
  inline_css=TRUE
  )
#Open the file with excel and mannually save as .xlsx
#EXPORT AS PNG TOO!!!

