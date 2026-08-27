### KM plots of OS by quantiles of MECOM expression level and Chr3 status ###

#Import inputs
km_qt_res<-qs2::qs_read(file.path(analyses_export_path,"km_results_qt.qs2"))

#Create folder and path to export results
km_plot_qt_path<-file.path(data_export_path,'km-plot-qt')
if (!dir.exists(km_plot_qt_path)) { dir.create(km_plot_qt_path,recursive=T) }

#Create a list of KM plots
km_plot_list_qt<-list(
  km_qt_res$by_3groups_os$plot,
  km_qt_res$by_3groups_os_all_chr3_alt$plot,
  km_qt_res$by_3groups_os_mecomr$plot,
  km_qt_res$by_3groups_os_chr3_abn$plot,
  km_qt_res$by_6groups_os$plot,
  km_qt_res$by_6groups_os_all_chr3_alt$plot,
  km_qt_res$by_6groups_os_mecomr$plot,
  km_qt_res$by_6groups_os_chr3_abn$plot,
  km_qt_res$by_2groups_os$plot,
  km_qt_res$by_2groups_os_all_chr3_alt$plot,
  km_qt_res$by_2groups_os_mecomr$plot,
  km_qt_res$by_2groups_os_chr3_abn$plot
)

#Generate multi-panel figure
comb_km_plots_qt<-ggarrange(
  plotlist=km_plot_list_qt,
  ncol=4,nrow=3
  #,labels="AUTO"
  ,font.label=list(size=20,face="bold")
  ,vjust=1.2,hjust=-0.25,heights=c(1,1,1)
)

## Export as PDF ##
ggsave(
  file.path(km_plot_qt_path,'km-plot-qt.pdf'),
  comb_km_plots_qt,
  width=11.69,height=8.27,units='in'
  )

## Export as PNG ##
ggsave(
  file.path(km_plot_qt_path,'km-plot-qt.png'),
  comb_km_plots_qt,
  width=11.69,height=8.27,units='in',dpi=300
)

