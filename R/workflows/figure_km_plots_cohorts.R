### External cohorts KM plots of OS by quantiles of RNA MECOM expression and Chr3 status ###

#Import inputs
km_cohorts<-qs2::qs_read(file.path(analyses_export_path,"km_plots_cohorts.qs2"))

#Create folder and path to export results
km_plot_cohorts_path<-file.path(data_export_path,'km-plot-cohorts')
if (!dir.exists(km_plot_cohorts_path)) { 
  dir.create(km_plot_cohorts_path,recursive=T) 
  }

#Create plot list
km_cohorts_list<-list(
  km_cohorts$beataml_os_2lvl$plot,
  km_cohorts$tcga_os_2lvl$plot,
  km_cohorts$beataml_os_2lvl_chr3$plot,
  km_cohorts$tcga_os_2lvl_chr3$plot
)
names(km_cohorts_list)<-names(km_cohorts)

#Generate multi-panel figure
comb_kmplts<-ggarrange(
  plotlist=km_cohorts_list,
  ncol=2,nrow=2,
  labels="AUTO",font.label=list(size=20,face="bold"))

## Export as PDF ##
ggsave(
  file.path(km_plot_cohorts_path,'km-plot-cohorts.pdf'),
  comb_kmplts,
  width=8.27,height=11.69,units='in'
  )

## Export as PNG ##
ggsave(
  file.path(km_plot_cohorts_path,'km-plot-cohorts.png'),
  comb_kmplts,
  width=8.27,height=11.69,units='in',dpi=300
  )
