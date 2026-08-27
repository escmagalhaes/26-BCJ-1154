### Pre-processing of inputs ###

##### Import data files ----
datatable      <-readr::read_csv(file.path(input_path,"data1.csv"),col_types=cols(.default="c"))
meta_table     <-readr::read_csv(file.path(input_path,"data2.csv"),col_types=cols(.default="c"))
normal_data    <-readr::read_csv(file.path(input_path,"data3.csv"),col_types=cols(.default="c"))
pfg_data       <-readr::read_csv(file.path(input_path,"data4.csv"),col_types=cols(.default="c"))
mdacc_rna_data <-readr::read_csv(file.path(input_path,"data5.csv"),col_types=cols(.default="c"))
beataml_data   <-readr::read_csv(file.path(input_path,"data6.csv"),col_types=cols(.default="c"))
tcga_data      <-readr::read_csv(file.path(input_path,"data7.csv"),col_types=cols(.default="c"))
kramer_data    <-readr::read_csv(file.path(input_path,"data8.csv"),col_types=cols(.default="c"))
jayavelu_data  <-readr::read_csv(file.path(input_path,"data9.csv"),col_types=cols(.default="c"))
ms_sub_set     <-readr::read_csv(file.path(input_path,"data10.csv"),col_types=cols(.default="c"))

##### Pre-processing Main datatable ----
datatable<-adjust_datatable(datatable=datatable,meta_table=meta_table)
datatable_num_cols<-setdiff(colnames(datatable),meta_table$variable)
datatable[,datatable_num_cols]<-lapply(datatable[,datatable_num_cols],as.numeric)
attr(datatable,"spec")<-NULL
attr(datatable,"problems")<-NULL

##### Pre-processing Meta table ----
meta_table_logi_cols<-c("include_table","include_uv_model","include_mv_model","forced")
meta_table[,meta_table_logi_cols]<-lapply(meta_table[,meta_table_logi_cols],as.logical)
meta_table$mv_order<-as.numeric(meta_table$mv_order)
attr(meta_table,"spec")<-NULL
attr(meta_table,"problems")<-NULL

##### Pre-processing NBM datatable ----
num_cols_normal_df<-setdiff(names(normal_data),'sample_type')
normal_data[,num_cols_normal_df]<-lapply(normal_data[,num_cols_normal_df],as.numeric)
normal_data$sample_type<-factor(normal_data$sample_type)
attr(normal_data,"spec")<-NULL
attr(normal_data,"problems")<-NULL

##### Pre-processing PFG datatable ----

#Transform pathway assignments into a list for downstream analysis
pfg_list<-membership_to_list(
  data           = pfg_data,
  membership_col = "pathway",
  yes_value      = "yes"
)

##### Pre-processing online datasets and MS data ----

#Portions of this section were adapted from Gosline, S. et al. (PMID: 35896960)
#Repository name: PNNL-CompBio/beatAMLpilotProteomics 
#Repository url:(https://github.com/PNNL-CompBio/beatAMLpilotProteomics/blob/main/beatAMLdata.R)
#Script under MIT License (See THIRD-PARTY-LICENSES.md for full license text)

if (nzchar(.synapse_token)) {
  
  #Helper function to fetch synapse data
  get_synapse_data<-function(syn_id,cache_path) {
    
    if (!dir.exists(cache_path)) dir.create(cache_path,recursive=TRUE)
    cache_file<-file.path(cache_path,paste0(syn_id,".qs2"))
    
    if (file.exists(cache_file)) {
      message("Loading ",syn_id," from cache...")
      data<-qs2::qs_read(cache_file)
      message("Done!")
      return(data)
    }
    
    message("Downloading ",syn_id," from Synapse...")
    data<-as.data.frame(synapser::synTableQuery(sprintf("SELECT * FROM %s",syn_id))) #SQL format
    message("Saving cache...")
    qs2::qs_save(data,cache_file)
    message("Done!")
    
    return(data)
  }
  
  #Login Authentication
  synapser::synLogin(authToken=.synapse_token)
  
  #Import and cache data from Synapse 
  syn22172602 <-get_synapse_data(syn_id="syn22172602",cache_path=synapse_cache_path) #Protein adn mRNA expression data
  syn22314121 <-get_synapse_data(syn_id="syn22314121",cache_path=synapse_cache_path) #Clinical data
  syn22170540 <-get_synapse_data(syn_id="syn22170540",cache_path=synapse_cache_path) #Clinical data
  cat("\nAll Done!\n")
  
  #Adjust protein and RNA data
  orig.data<-syn22172602 |> dplyr::rename(
    proteinLevels = 'LogFoldChange',
    mRNALevels    = 'transcriptCounts',
    geneMutations = 'Tumor VAF') |>
    mutate(Gene = unlist(Gene)) |> rowwise() |>
    mutate(binaryMutations = ifelse(geneMutations==0,0,1)) |> as_tibble()
  
  #Import clinical data
  pat.data<-syn22314121 |> mutate(Gene = unlist(Gene)) |>
    subset(Treatment == 'Vehicle') |> subset(`Cell number` >= 10000000) |>
    dplyr::select(Gene,LogFoldChange,`AML sample`) |> distinct() |> 
    as_tibble()
  
  #Expression data
  expr.data<-pat.data |>
    full_join(orig.data,by=c('AML sample','Gene'),relationship="many-to-many") |> 
    rowwise() |> dplyr::mutate(proteinLevels = max(proteinLevels,LogFoldChange,na.rm=T)) |>
    dplyr::select(-LogFoldChange) |> dplyr::mutate(
      mRNALevels      = tidyr::replace_na(mRNALevels,0),
      geneMutations   = tidyr::replace_na(geneMutations,0),
      binaryMutations = tidyr::replace_na(binaryMutations,0)) |>
    dplyr::select(-countMetric) |> distinct() |> as_tibble()
  
  #Determine range of mRNALevels in dataset
  #If 0 until thousands: not Log2 normalized
  #If 0 until ~20: Log2 normalized
  #If negative values: already Log2 normalized AND scaled
  range_expr.data<-range(expr.data |> select(mRNALevels))
  
  #Log2 normalize the data
  expr.data_norm<-expr.data |> dplyr::mutate(across(.cols=mRNALevels,.fns=~log2(.x + 1)))
  range_expr.data_norm<-range(expr.data_norm |> select(mRNALevels))
  
  #Get only patients with protein levels
  pats.with.prot<-expr.data_norm |>
    group_by(`AML sample`) |> summarize(hasProt=all(proteinLevels==0)) |>
    subset(hasProt==FALSE) |> dplyr::select('AML sample')
  expr.data.with.prot<-expr.data_norm |> subset(`AML sample`%in% pats.with.prot$`AML sample`)
  
  #Clinical data
  pat.drugClin<-syn22170540 |> mutate(Condition=unlist(Condition)) |> 
    dplyr::select(-ROW_ID,-ROW_VERSION,-priorMalignancyType,-Condition,-Metric,-Value) |> 
    distinct() |> dplyr::select(
      upi       = "AML sample",
      gender    = "gender",
      age       = "ageAtDiagnosis",
      status    = "vitalStatus",
      surv_time = "overallSurvival"
    ) |> dplyr::mutate(
      gender     = case_when(
        gender   == "Female" ~ "female",
        gender   == "Male"   ~ "male",
        .default = gender
      ),
      status     = case_when(
        status == "Unknown" ~ NA_character_,
        status == "Alive"   ~ "0",
        status == "Dead"    ~ "1",
        .default = status
      ),
      status     = as.integer(status),
      surv_time  = as.numeric(surv_time),
      surv_time  = round(surv_time/365,digits=2)) |> as_tibble() 
  
  #Filter for MECOM only and adjust data
  expr_filt<-expr.data_norm |> filter(grepl("MECOM",Gene,ignore.case=TRUE)) |>
    dplyr::select(
      upi          = "AML sample",
      "ptn_MECOM"  = proteinLevels,
      "mrna_MECOM" = mRNALevels) |> as_tibble()
  
  #Merge expression with clinical data and add cohort variable
  pat.data.full<-left_join(expr_filt,pat.drugClin,by='upi') |> as_tibble()
  
  #Remove NAs from data and select variables
  pat.data.clean<-pat.data.full |> filter(!is.na(status),!is.na(surv_time)) |>
    dplyr::select(upi,status,surv_time,mrna_MECOM,ptn_MECOM)
  
  #Check if only one patient per data point (If TRUE variable is unique)
  !anyDuplicated(pat.data.clean$upi)
  
  #Select only protein data and relevant variables and adjust dataset
  pat.data.final<-pat.data.clean |> dplyr::select(-mrna_MECOM,MECOM=ptn_MECOM) |> 
    dplyr::mutate(
      type  = "protein (MS)",
      type  = factor(type),
      MECOM = case_when(
        MECOM    == 0 ~ NA_real_,
        .default = MECOM),
      MECOM_zscore = as.vector(scale(MECOM))) |> dplyr::mutate(
      cohort = "Gosline, S. et al. (BeatAML)",
      cohort = factor(cohort)) |> dplyr::select(cohort,everything())
  
} else {
  pat.data.final<-NULL
  warning(
  "Synapse Login unsuccessful. Token variable is empty. 
  Double check .synapse_token assignment inside script config.R.
  Remember to create a Synapse account and a Personal Access Token (PAT). 
  PAT should be a long string of characters/numbers, not the name you gave to it.
  For more details, check out: 
  https://docs.synapse.org/synapse-docs/managing-your-account#ManagingYourAccount-PersonalAccessTokens"
  )
}

#Merge all cohorts in a single dataset
if(!is.null(pat.data.final)){
  pat.data.final.char<-pat.data.final |> mutate(across(everything(),as.character))
  bind_cohorts<-bind_rows(beataml_data,tcga_data,kramer_data,jayavelu_data,ms_sub_set,
                          mdacc_rna_data,pat.data.final.char)
} else {
  bind_cohorts<-bind_rows(beataml_data,tcga_data,kramer_data,jayavelu_data,ms_sub_set,mdacc_rna_data)
}
num_cols_bind_cohorts <-bind_cohorts |> 
  dplyr::select(status,surv_time,MECOM,MECOM_zscore) |> names()
fct_cols_bind_cohorts <-setdiff(names(bind_cohorts),c("upi","karyotype",num_cols_bind_cohorts))
bind_cohorts_adj<-bind_cohorts |> dplyr::mutate(
  dplyr::across(dplyr::all_of(num_cols_bind_cohorts),~as.numeric(.x)),
  dplyr::across(dplyr::all_of(fct_cols_bind_cohorts),~factor(.x)))
rppa_subset<-datatable |> dplyr::select(
  upi,status,surv_time,karyotype,chr3_status_cat,chr3_status,chr3_abn,MECOM_R,MECOM) |>  
  dplyr::mutate(
    MECOM_zscore = as.vector(scale(MECOM)),
    type   = "protein (RPPA)",
    type   = factor(type),
    cohort = "UTMDACC",
    cohort = factor(cohort)) |> dplyr::select(
      cohort,type,upi,status,surv_time,karyotype,chr3_status_cat,
      chr3_status,chr3_abn,MECOM_R,MECOM,MECOM_zscore)
cohort_mrg<-bind_rows(rppa_subset,bind_cohorts_adj)

##### Create patient subsets ----

#Create object with protein names
all_protein_names<-setdiff(colnames(datatable),meta_table$variable)
#easiest way to extract all protein names from main dataset

#Split data into several quantiles according to protein expression
ptn_gp_data<-get_quantile_groups(
  data        = datatable,
  features    = c('MECOM','CEBPA'),
  n_quantiles = 2:10,
  data_normal = normal_data,
  mode        = "both"
  )

#Regroup quantiles
regrouping_quantiles<-tibble::tibble(
  feature = c("MECOM_gp_n","MECOM_gp_n","MECOM_gp6","CEBPA_gp5"),
  groups  = list(
    list('Low'=c(1,2),'High'=c(3,4)),
    list('Below Normal'=1,'Lower Normal'=2,'Higher than Normal'=c(3,4)),
    list('Low 5/6ths'=c(1:5),'High 6th'=6),
    list('Low'=c(1,2),'High'=c(3:5))
  ),
  new_col = c("MECOM_gp_n2","MECOM_gp_n3","MECOM_gp6_n2","CEBPA_gp_n2")
)

for (i in seq_len(nrow(regrouping_quantiles))) {
  ptn_gp_data<-collapse_quantile_groups(
    data    = ptn_gp_data,
    feature = regrouping_quantiles$feature[ i ],
    groups  = regrouping_quantiles$groups[[ i ]],
    new_col = regrouping_quantiles$new_col[ i ]
  )
}

#Create variable according to CH3 status and MECOM levels
ptn_gp_data<-ptn_gp_data %>%
  dplyr::mutate(
    MECOM_gp_n2_chr3_status = dplyr::case_when(
      is.na(MECOM_gp_n2) | is.na(chr3_status) | is.na(chr3_status_cat) ~ NA_character_,
      chr3_status     == 'Chr3-WT'  & MECOM_gp_n2 == 'Low'  ~ "Chr3-WT Low",
      chr3_status     == 'Chr3-WT'  & MECOM_gp_n2 == 'High' ~ "Chr3-WT High",
      chr3_status     == 'Chr3-abn' & MECOM_gp_n2 == 'Low'  ~ "Chr3-abn Low",
      chr3_status     == 'Chr3-abn' & MECOM_gp_n2 == 'High' ~ "Chr3-abn High",
      chr3_status_cat == 'MECOM-R'           ~ 'MECOM-R',
      chr3_status_cat == 'atypical MECOM-R'  ~ 'atypical MECOM-R'),
    MECOM_gp_n2_chr3_status = factor( 
      MECOM_gp_n2_chr3_status,
      levels=c("Chr3-WT Low","Chr3-WT High","MECOM-R",
               "atypical MECOM-R","Chr3-abn Low","Chr3-abn High")))

#Create subset of main data without NAs for CH3 status
ptn_gp_data_filt<-ptn_gp_data[!is.na(ptn_gp_data$chr3_status),]

#Create patient subsets according to grouping variables
#For therapy subset, dataset should have all cases
subset_groups<-tibble::tibble(
  df           = list(ptn_gp_data_filt,ptn_gp_data,ptn_gp_data_filt,ptn_gp_data_filt),
  var          = c("chr3_status","chem3","MECOM_gp_n2","MECOM_gp_n"),
  var_labs     = c("Chr3_status","therapy","MECOM_2lvl_expr","MECOM_4lvl_expr"),
  ref_patterns = list(c("WT","wt","no"),c("other","Other"),c("low","Low"),c("Below Normal","below normal")),
  ref_suffix   = c("WT","other","Low",'Below_Normal'),
  alt_suffix   = c("altered","other_excluded","High",'Below_Normal_excluded')
)

subsets_list<-list()
for (i in seq_len(nrow(subset_groups))) {
  df_to_use<-subset_groups$df[[ i ]]
  grouped_subsets<-get_var_groups(
    df           = df_to_use,
    var          = subset_groups$var[ i ],
    var_labs     = subset_groups$var_labs[ i ],
    ref_patterns = subset_groups$ref_patterns[[ i ]],
    ref_suffix   = subset_groups$ref_suffix[[ i ]],
    alt_suffix   = subset_groups$alt_suffix[[ i ]]
  )

  #Create subsets and droplevels of the variable used to create subsets
  subsets<-lapply(grouped_subsets,function(mask) {
    sub_df<-df_to_use[mask,,drop=FALSE]
    labelled::var_label(sub_df)<-meta_table$label[match(names(sub_df),meta_table$variable)]
    sub_df
  })
  subsets_list[[ subset_groups$var[ i ] ]]<-Filter(function(s) nrow(s) > 0, subsets)
}
names(subsets_list)<-subset_groups$var_labs

#Create grouped CH3 subsets 
subsets_list$Chr3_status$Chr3_status_WT_altered<-dplyr::bind_rows(
  subsets_list$Chr3_status$Chr3_status_WT,
  subsets_list$Chr3_status$Chr3_status_altered
)

subsets_list$Chr3_status$Chr3_status_WT_MECOMR<-dplyr::bind_rows(
  subsets_list$Chr3_status$Chr3_status_WT,
  subsets_list$Chr3_status$Chr3_status_MECOMR
)
subsets_list$Chr3_status$Chr3_status_WT_MECOMR$chr3_status<-droplevels(
  subsets_list$Chr3_status$Chr3_status_WT_MECOMR$chr3_status)

subsets_list$Chr3_status$Chr3_status_WT_Chr3abn<-dplyr::bind_rows(
  subsets_list$Chr3_status$Chr3_status_WT,
  subsets_list$Chr3_status$Chr3_status_Chr3abn
)
subsets_list$Chr3_status$Chr3_status_WT_Chr3abn$chr3_status<-droplevels(
  subsets_list$Chr3_status$Chr3_status_WT_Chr3abn$chr3_status)

#Create subset with IC vs VH therapies only and adjust MECOM/CEBPA variables
subsets_list$therapy$therapy_IC_VH<-dplyr::bind_rows(
  subsets_list$therapy$therapy_IC_only,
  subsets_list$therapy$therapy_IC_VTX,
  subsets_list$therapy$therapy_VH
)
#Simplify therapy levels
subsets_list$therapy$therapy_IC_VH$chem3<-factor(
  ifelse(subsets_list$therapy$therapy_IC_VH$chem3=='VH','VH','IC'),levels=c('IC','VH'))

#For MECOM: Low = favorable, High = unfavorable, regardless of therapy
#For CEBPA: High IC and Low VH = favorable, Low IC and High VH = unfavorable
favorable_conditions<-tibble::tibble(
  protein = c("MECOM","MECOM","CEBPA","CEBPA"),
  therapy = c( "IC"  , "VH"  , "IC"  , "VH"  ),
  ptn_lvl = c( "Low" , "Low" , "High", "Low" )
)

for ( ptn in unique(favorable_conditions$protein)) {
  therapy_col<-"chem3"
  ptn_lvl_col<-paste0( ptn ,"_gp_n2")
  fav_conditions<-favorable_conditions[favorable_conditions$protein== ptn , ]
  
  prog_col<-paste0(ptn,"_prog")
  prog_num_col<-paste0(ptn,"_prog_num")
  
  subsets_list$therapy$therapy_IC_VH<-subsets_list$therapy$therapy_IC_VH %>%
    dplyr::mutate(
      !!prog_col := dplyr::case_when(
        is.na(.data[[ptn_lvl_col]]) | is.na(.data[[therapy_col]])~NA_character_,
        paste(.data[[therapy_col]],.data[[ptn_lvl_col]]) %in% 
          paste(fav_conditions$therapy,fav_conditions$ptn_lvl)~"favorable",
        TRUE~"unfavorable"
      ),
      !!prog_num_col := dplyr::case_when(
        is.na(.data[[ptn_lvl_col]]) | is.na(.data[[therapy_col]])~ NA_real_,
        paste(.data[[therapy_col]],.data[[ptn_lvl_col]]) %in% 
          paste(fav_conditions$therapy,fav_conditions$ptn_lvl)~0,
        TRUE~1
      ),
      !!prog_col := factor(.data[[prog_col]],levels=c("favorable","unfavorable"))
    )
}

#Create subsets without relapse=NA cases
relapse_subsets_list<-lapply(subsets_list,function(var_subsets) {
  lapply(var_subsets,function(subset) {
    sub_df<-subset[!is.na(subset$relapse),]
  })
})
ptn_gp_data_filt_rem<-ptn_gp_data_filt[!is.na(ptn_gp_data_filt$relapse),]

##### Save processed inputs ----

#Main datasets
qs2::qs_save(
  list(
    main_df = ptn_gp_data_filt,
    df_rem  = ptn_gp_data_filt_rem
    ),file.path(analyses_export_path,"processed_main_data.qs2"))

#Reference datatsets for main dataset
qs2::qs_save(
  list(
    meta_table   = meta_table,
    normal_data  = normal_data,
    pfg_list     = pfg_list,
    all_proteins = all_protein_names
    ),file.path(analyses_export_path,"processed_reference_sets.qs2"))

#Subsets of main data
qs2::qs_save(
  list(
    all     = subsets_list,
    relapse = relapse_subsets_list
    ),file.path(analyses_export_path,"processed_subsets.qs2"))

#Cross-study dataset
qs2::qs_save(
  cohort_mrg
  ,file.path(analyses_export_path,"processed_cross_study_data.qs2"))




