### Extra Scatter plots of Correlation networks ###

#Import inputs
scatter_plots<-qs2::qs_read(file.path(analyses_export_path,"scatter_plots.qs2"))

#Create folders and paths to export results
sct_plots_path<-file.path(data_export_path,'scatter-plots')
if (!dir.exists(sct_plots_path)) { dir.create(sct_plots_path,recursive=T) }

#Extract Scatter plots
sc_plots<-lapply(scatter_plots, function(x) x$scplot)

#Create folders and paths to export results
sct_plots_extra_path<-file.path(data_export_path,'scatter-plots-extra')
if (!dir.exists(sct_plots_extra_path)) { dir.create(sct_plots_extra_path,recursive=T) }

sc_res_plots_extra<-sc_plots[c(
  'Chr3-WT vs MECOM-R',
  'Chr3-WT vs Chr3-abn',
  'MECOM-R vs Chr3-abn',
  'Chr3-WT High vs Chr3-WT Low'
)]

#Generate multi-panel figure
comb_sc_plots_extra<-ggarrange(
  plotlist   = sc_res_plots_extra,
  ncol       = 2,
  nrow       = 2,
  labels     = "AUTO",
  font.label = list(size=20,face="bold")
)

## Export as PDF ##
ggsave(
  file.path(sct_plots_extra_path,"scatter-plots-extra.pdf"),
  comb_sc_plots_extra
  ,height=11.69,width=8.27,units='in'
)

## Export as PNG ##
ggsave(
  file.path(sct_plots_extra_path,"scatter-plots-extra.png"),
  comb_sc_plots_extra
  ,height=11.69,width=8.27,units='in',dpi=300
)




