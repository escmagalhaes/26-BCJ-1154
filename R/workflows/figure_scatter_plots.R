### Scatter plots of Correlation networks ###

#Import inputs
scatter_plots<-qs2::qs_read(file.path(analyses_export_path,"scatter_plots.qs2"))

#Create folders and paths to export results
sct_plots_path<-file.path(data_export_path,'scatter-plots')
if (!dir.exists(sct_plots_path)) { dir.create(sct_plots_path,recursive=T) }

#Extract Scatter plots
sc_plots<-lapply(scatter_plots, function(x) x$scplot)

#Combine selected Scatter plots into Figures
sc_res_plots<-sc_plots[c(
  'MECOM-R vs Chr3-WT High',
  'MECOM-R vs Chr3-abn High',
  'Chr3-WT High vs Chr3-abn High',
  'MECOM-R vs Chr3-WT Low',
  'Chr3-WT Low vs Chr3-abn Low',
  'Chr3-WT Low vs Chr3-abn High'
)]

#Generate multi-panel figure
comb_sc_plots<-ggarrange(
  plotlist   = sc_res_plots,
  ncol       = 2,
  nrow       = 3,
  labels     = "AUTO",
  font.label = list(size=20,face="bold")
  )

## Export as PDF ##
ggsave(
  file.path(sct_plots_path,"scatter-plots.pdf"),
  comb_sc_plots
  ,height=11.69,width=8.27,units='in'
)

## Export as PNG ##
ggsave(
  file.path(sct_plots_path,"scatter-plots.png"),
  comb_sc_plots
  ,height=11.69,width=8.27,units='in',dpi=300
)




