### Figures related to MECOM/CEBPA prognostic models ###

#Import inputs
prog_mods<-qs2::qs_read(file.path(analyses_export_path,"prognostic_models.qs2"))

#Create local objects as inputs
bxps           <-prog_mods$boxplots
km_plts        <-prog_mods$km_plots
mv_mod_tab     <-prog_mods$mv_table
df_prog_cl     <-prog_mods$dataset_clusters
meta_tab       <-prog_mods$meta_table
cluster_colors <-get_mycolors(n_colors=nlevels(df_prog_cl$cluster),mode='vector')
cluster_colors <-setNames(cluster_colors,levels(df_prog_cl$cluster))

#Create folder and path to export results
prog_mods_path<-file.path(data_export_path,'prog-mods')
if (!dir.exists(prog_mods_path)) { dir.create(prog_mods_path,recursive=T) }

## Boxplots ##
bxp_prog_list<-unlist(lapply(bxps, function(x) 
  lapply(x, function(y) y$bxplot)),recursive=FALSE)
names(bxp_prog_list)<-gsub("\\..*", "",names(bxp_prog_list))

#Combine Boxplots in a single Figure
comb_bxplot_prog_list<-ggarrange(
  plotlist   = bxp_prog_list,
  ncol       = 1,
  nrow       = 2,
  labels     = "AUTO",
  font.label = list(size=20,face="bold")
)

## Heatmap ##
ht_covars=c(
  'chem2','vtx','aml_gp','cyto_risk','complex_kar','chr3_status','mut_asxl1',
  'mut_cebpa','mut_dnmt3a','mut_ezh2','mut_flt3','mut_flt3_itd','mut_flt3_d835',
  'mut_idh1','mut_idh2','mut_jak2','mut_kit','mut_kmt2a','mut_npm1','mut_ras',
  'mut_kras','mut_nras','mut_phf6','mut_ptpn11','mut_rad21','mut_runx1',
  'mut_smc1a','mut_smc3','mut_stag2','mut_wt1'
)
ht_covar_labels<-as.vector(meta_tab[match(ht_covars,meta_tab$variable),]$label)
prog_ht<-plot_heatmap(
  df=df_prog_cl
  ,proteins=c('MECOM_prog_num','CEBPA_prog_num')
  ,protein_labels=c('MECOM','CEBPA')
  ,grouping_var='cluster'
  ,grouping_var_label='Cluster'
  ,col_anno=TRUE
  ,covars=ht_covars
  ,covar_labels=ht_covar_labels
  ,covar_colors_override=list(
    'Cytogenetic Risk'=c(
      'favorable'='red',
      'intermediate'=pal_jco()(2)[1],
      'unfavorable'=pal_jco()(2)[2]))
  ,colors_ht=c('royalblue','yellow')
  ,colors_grouping_var=cluster_colors
  ,colors_NA=pal_jco()(3)[3]
  ,plot_title="IC vs VH AML patients (N=469)"
  ,cluster_rows=TRUE,
  cluster_cols=FALSE,
  show_rownames=TRUE
  ,legend_breaks=c(0,1),
  legend_labels=c('Favorable','Unfavorable')
  ,show_colnames=FALSE,
  annotation_names_row=FALSE
  )

## KM plots ##
km_plot_os<-ggarrange(
  km_plts$prog_OS$plot,
  km_plts$prog_OS$pval_table,
  ncol=1,
  nrow=2,
  heights=c(0.85,0.15)
  )
km_plot_rd<-ggarrange(
  km_plts$prog_RD$plot,
  km_plts$prog_RD$pval_table,
  ncol=1,
  nrow=2,
  heights=c(0.85,0.15)
  )

## MV model table ##
mv_tab_mod<-mv_mod_tab |> modify_header(label~"**Variable**") |>
  modify_caption("**Overall Survival**") |> add_forest() |>
  gt::sub_missing(columns=conf.low,rows=2,missing_text="_________") |> 
  gt::sub_missing(columns=p.value,rows=2,missing_text="reference") |> 
  gt::tab_style(style=cell_text(v_align="top")
                ,locations=cells_body(columns=conf.low,rows=2)) |> 
  gt::tab_style(style=list(cell_text(color=cluster_colors['C1'],weight="bold"))
                ,locations=cells_body(columns=label,rows=2)) |>
  gt::tab_style(style=list(cell_text(color=cluster_colors['C2'],weight="bold"))
                ,locations=cells_body(columns=label,rows=3)) |>
  gt::tab_style(style=list(cell_text(color=cluster_colors['C3'],weight="bold"))
                ,locations=cells_body(columns=label,rows=4))

### Export Main Figure Panels ###
ggsave(file.path(prog_mods_path,'prog-mods-boxplots.pdf'),comb_bxplot_prog_list,width=8,height=6,units='in')
ggsave(file.path(prog_mods_path,'prog-mods-heatmap.pdf'),prog_ht,width=11.69,height=8.27,units='in')
ggsave(file.path(prog_mods_path,'prog-mods-km-plot-os.pdf'),km_plot_os,width=6,height=8,units='in')
ggsave(file.path(prog_mods_path,'prog-mods-km-plot-rd.pdf'),km_plot_rd,width=6,height=8,units='in')
gtsave(mv_tab_mod,file.path(prog_mods_path,"prog-mods-mv-table.pdf"),vwidth=1754,vheight=1240,expand=0)
#Open files in illustrator and adjust panels mannually
#Create a final PDF with all panels named 'prog-mods.pdf'
#Remember to also export 'prog-mods.pdf' as PNG!!!


