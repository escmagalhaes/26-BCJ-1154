### Meta-analysis of CoxPH UV models of OS by MECOM expression levels ###

#Import input data
meta_model<-qs2::qs_read(file.path(analyses_export_path,"meta_analysis_cohorts.qs2"))

#Create local objects
meta_mod <-meta_model$meta_mod
meta_dat <-meta_model$meta_df

#Create folder and path to export results
meta_analysis_path<-file.path(data_export_path,'meta-analysis')
if (!dir.exists(meta_analysis_path)) { 
  dir.create(meta_analysis_path,recursive=T) 
}

#Create a function to produce the meta-analysis plot
meta_plot<-function(){
  metafor::forest(
  meta_mod,
  ilab = cbind(
    meta_dat$Total,
    meta_dat$Missing,
    meta_dat$`MECOM-R`,
    meta_dat$`Chr3-abn`,
    meta_dat$Model,
    meta_dat$Event,
    meta_dat$`p-value`),
  ilab.xpos = c(0.45,0.6,0.75,0.9,1.05,1.2,1.35),
  slab      = meta_dat$cohort,
  atransf   = exp,
  xlab      = "Hazard Ratio", 
  refline   = 0,
  xlim      = c(-0.7,1.7), #room for extra columns
  alim      = c(-0.1,0.3), 
  header    = FALSE,
  cex       = 0.75
)

#Add column headers mannually
text(
  x = c(-0.625,0.45,0.6,0.75,0.9,1.05,1.2,1.35,1.6),
  y = (meta_mod$k + 1.25),
  labels = c(
    "Cohort","Total","Missing","MECOM-R","Chr3-abn",
    "Model","Event","p-value","HR [95% CI]"
    ),
  font = 2,
  cex  = 0.75
  )

#Add pooled data mannually
text(
  x = c(0.45,0.6,0.75,0.9,1.05,1.2,1.35),
  y = -1,
  labels = c(
    sum(meta_dat$Total),
    sum(meta_dat$Missing),
    sum(meta_dat$`MECOM-R`),
    sum(meta_dat$`Chr3-abn`),
    sum(meta_dat$Model),
    sum(meta_dat$Event),
    format.pval(meta_mod$pval,digits=2,eps=0.001)),
  font = 2,
  cex  = 0.75
  )
#text(x=1,y= -1,labels=paste0("I² = ",round(meta_mod$I2,2),"%"),font=2)

#Add Title
title(main = expression(atop(
  bold("Meta-analysis of Univariate Cox Proportional Hazards of "),
  bold("MECOM expression and Overall Survival (OS) across cohorts"))))
}

## Export as PDF ##
pdf(file.path(meta_analysis_path,"meta-analysis.pdf"),width=11.69,height=8.27)
meta_plot()
dev.off()

## Export as PNG ##
png(
  file.path(meta_analysis_path,"meta-analysis.png"),
  width=11.69,height=8.27,units="in",res=300)
meta_plot()
dev.off() 


