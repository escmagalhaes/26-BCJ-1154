### Boxplot analyses ###

#Import inputs
subsets<-qs2::qs_read(file.path(analyses_export_path,"processed_subsets.qs2"))

#Create varying parameters list 
df_bxp<-list(
  'all cases' = subsets$all$Chr3_status$all,
  'in Chr3 WT' = subsets$all$Chr3_status$Chr3_status_WT
)

#Create fixed parameters list
fixed_bxp_params_df_bxp<-list(
  continuous_vars = 'MECOM',
  ylab            = 'MECOM levels',
  discrete_vars   = c(
    "cyto_risk","cyto_cat","fab2","mut_gata2","mut_sf3b1",
    "mut_npm1","mut_idh","mut_cebpa","mut_cebpa_db"
  ),
  disc_var_labels = c(
    'Cytogenetic Risk','Cytogenetic Category','FAB Classification',
    'GATA2 mutation','SF3B1 mutation','NPM1 mutation','IDH mutation',
    'CEBPA mutation','CEBPA double mutation'
  ),
  colors          = get_mycolors(mode = 'vector'),
  angle_x_var     = 45
)

#Create all boxplots with a single call
bxp_res<-purrr::imap(df_bxp, function( dataframe , name) {
  p<-list(df = dataframe , title_suffix = name )
  args<-c( fixed_bxp_params_df_bxp , p )
  do.call( run_bxplots , args )
})

#Interleave results
inter_bxp_res<-interleave_results(
  list_of_lists = bxp_res,
  group_names   = paste0('MECOM_',c(
    'Cytogenetic Risk','Cytogenetic Category','FAB Classification',
    'GATA2 mutation','SF3B1 mutation','NPM1 mutation','IDH mutation'
    ,'CEBPA mutation','CEBPA double mutation')
  ),
  suffixes      = c(' all cases',' in Chr3 WT')
)

### Save analyses ###
qs2::qs_save(inter_bxp_res,file.path(analyses_export_path,"boxplots.qs2"))


