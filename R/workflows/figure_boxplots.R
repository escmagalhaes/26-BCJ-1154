### Boxplots of categorical variables by MECOM expression levels ###

#Import inputs
bxplts<-qs2::qs_read(file.path(analyses_export_path,"boxplots.qs2"))

#Create folders and paths to export results
boxplots_path<-file.path(data_export_path,'boxplots')
if (!dir.exists(boxplots_path)) { dir.create(boxplots_path,recursive=T) }

#Extract boxplots 
bxp_res_list<-lapply(bxplts, function( x ) x$bxplot)

#Generate multi-panel figure
comb_bxplot_list<-ggpubr::ggarrange(
  plotlist=bxp_res_list,
  #labels="AUTO",
  ncol=2,nrow=3,
  font.label=list(size=20,face="bold")
)

## Export as PDF ##
## Remember to open and export as PNG too!!! ##
pdf(file.path(boxplots_path,"boxplots.pdf"),width=8.27,height=11.69)
print(comb_bxplot_list)
dev.off()


