### KM plots of RD by MECOM expression level and Chr3 status ###

#Import inputs
km_rd_res<-qs2::qs_read(file.path(analyses_export_path,"km_results_rd.qs2"))

#Create folder and path to export results
km_plot_rd_path<-file.path(data_export_path,'km-plot-rd')
if (!dir.exists(km_plot_rd_path)) { dir.create(km_plot_rd_path,recursive=T) }

#Create a list of KM plots
km_plot_list_rd<-list(
  km_rd_res$by_4groups_rd$plot,
  km_rd_res$rd_below_normal_chr3_type$plot,
  km_rd_res$rd_lower_normal_chr3_type$plot,
  km_rd_res$rd_upper_normal_chr3_type$plot,
  km_rd_res$rd_above_normal_chr3_type$plot
)

#Generate multi-panel figure
comb_km_plots_rd<-ggarrange(
  plotlist=km_plot_list_rd,
  ncol=2,nrow=3,
  labels="AUTO",font.label=list(size=20,face="bold"),
  vjust=1.2,hjust=-0.25,heights=c(1,1,1)
)

## Export as PDF ##
ggsave(
  file.path(km_plot_rd_path,'km-plot-rd.pdf'),
  comb_km_plots_rd,
  width=8.27,height=11.69,units='in'
  )

## Export as PNG ##
ggsave(
  file.path(km_plot_rd_path,'km-plot-rd.png'),
  comb_km_plots_rd,
  width=8.27,height=11.69,units='in',dpi=300
)
