### Density plot and KM plots of OS by MECOM expression level and Chr3 status ###

#Import inputs
km_os_res <-qs2::qs_read(file.path(analyses_export_path,"km_results_os.qs2"))
main_data <-qs2::qs_read(file.path(analyses_export_path,"processed_main_data.qs2"))
ref_sets  <-qs2::qs_read(file.path(analyses_export_path,"processed_reference_sets.qs2"))

#Create local objects
main_df   <-main_data$main_df
normal_df <-ref_sets$normal_data

#Create folder and path to export results
km_plot_os_path<-file.path(data_export_path,'km-plot-os')
if (!dir.exists(km_plot_os_path)) { dir.create(km_plot_os_path,recursive=T) }

#Density plot
dens_plot<-plot_density(
  df              = main_df,
  df_normal       = normal_data,
  protein         = 'MECOM',
  grouping_var    = 'MECOM_gp_n',
  line_var_label  = 'Sample Type',
  color_fill_mode ='list',
  legend_mode     = "inside",
  n_panels        = 6,
  legend_position = c(0.75,0.7)
)

#Create lists of KM plots
plot_list_os<-list(
  dens_plot, 
  ggarrange(
    km_os_res$by_4groups_os$plot,
    km_os_res$by_4groups_os$pval_table,
    ncol=1,nrow=2,heights=c(0.75,0.25)),
  ggarrange(
    km_os_res$os_below_normal_chr3_type$plot,
    km_os_res$os_below_normal_chr3_type$pval_table,
    ncol=1,nrow=2,heights=c(0.825,0.175)),
  km_os_res$os_lower_normal_chr3_type$plot,
  km_os_res$os_upper_normal_chr3_type$plot,
  km_os_res$os_above_normal_chr3_type$plot
)

#Generate multi-panel figure
comb_dens_km_plots_os<-ggarrange(
  plotlist=plot_list_os,
  ncol=2,nrow=3,
  labels="AUTO",font.label=list(size=20,face="bold"),
  vjust=1.2,hjust=-0.25,heights=c(1,1,1)
)

## Export as PDF ##
ggsave(
  file.path(km_plot_os_path,'km-plot-os.pdf'),
  comb_dens_km_plots_os,
  width=8.27,height=11.69,units='in'
  )

## Export as PNG ##
ggsave(
  file.path(km_plot_os_path,'km-plot-os.png'),
  comb_dens_km_plots_os,
  width=8.27,height=11.69,units='in',dpi=300
  )


