### Helper Functions ###

#Function to make tibble with varying parameters for KM plots
get_varying_km_params<-function(df_list,labs,title,prefix,time_var,event_var,
                                grouping_var=NULL,single_color=FALSE,
                                color_list=NULL,color_vector=NULL,
                                color_mode=c("list","vector"),
                                covar_list=NULL,n_panels=1) {
  
  mode<-match.arg(color_mode)
  
  #Safety checks
  if (isTRUE(single_color) && mode=="vector")
    stop("single_color=TRUE is only valid when color_mode='list'")
  
  if (!is.null(covar_list) && length(covar_list) > 1 && length(covar_list) != length(df_list))
    stop("covar must be length 1 or the same length as df_list.")
  
  #Normalize covar to list
  if (!is.list(covar_list)) {
    if (is.null(covar_list)) {
      covar_list<-vector("list",length(df_list))
    } else {
      covar_list<-as.list(rep(covar_list,length.out=length(df_list)))
    }
  }
  
  tab<-tibble::tibble(
    name         = paste(prefix,names(df_list),sep="_"),
    data         = unname(df_list),
    label        = labs,
    time_var     = time_var,
    event_var    = event_var,
    covar        = covar_list,
    n_panels     = n_panels
  ) %>% dplyr::mutate(title=trimws(paste(title,label,sep='\n'))) %>% dplyr::select(-label)
  
  if (mode=="list") {
    if (is.null(color_list)) color_list<-get_mycolors(mode="list")
    
    if (isTRUE(single_color)) {
      if (is.null(grouping_var))
        stop("grouping_var must be provided to compute single_color.")
      
      #Create a single color vector and add it to the table
      single_color_vec<-purrr::map2_chr(
        df_list,seq_along(df_list),
        ~get_mycolors(n_colors=nlevels(.x[[grouping_var]]),mode="list",color_list=color_list)[.y])
      
      tab<-tab %>% dplyr::mutate(grouping_var  = list(NULL)
                                 ,colors       = list(NULL)
                                 ,single_color = as.list(single_color_vec))
    } else {
      if (is.null(grouping_var))
        stop("grouping_var must be provided when single_color=FALSE")
      
      color_vec<-purrr::map(
        df_list,
        ~get_mycolors(n_colors=nlevels(.x[[grouping_var]]),mode="list",color_list=color_list))
      
      tab<-tab %>% dplyr::mutate(grouping_var = as.list(rep(grouping_var,length.out=nrow(tab))),
                                 colors       = color_vec,
                                 single_color = list(NULL))
    }
  } else {
    
    #grouping_var passed through, single_color=FALSE
    if (is.null(color_vector)) color_vector<-get_mycolors(mode="vector")
    if (is.null(grouping_var))
      stop("grouping_var must be provided when color_mode='vector'")
    
    color_vec<-purrr::map(df_list,~color_vector)
    
    tab<-tab %>% dplyr::mutate(grouping_var = as.list(rep(grouping_var,length.out=nrow(tab))),
                               colors       = color_vec,
                               single_color = list(NULL))
  }
  return(tab)
}


#Function to make specific gtsummary UV tables
get_uv_table<-function(data,features,labels=NULL,count_vars=NULL,count_labels=NULL,
                       sum_value=NULL,adjust_CI=FALSE,max_CI=50,min_CI=0.1,
                       time_var="surv_time",event_var="status"){
  
  #Coerce features to characters and create model fit
  features<-as.character(features)
  
  #Safety checks
  missing_cols<-setdiff(c(features,time_var,event_var),names(data))
  if (length(missing_cols) > 0) {
    stop("Columns not found in data: ",paste(missing_cols,collapse=", "))
  }
  if (!is.null(labels) && length(labels) != length(features)) {
    stop("`labels` must have same length as `features`.")
  }
  if (length(features)==0) {
    stop("`features` cannot be empty.")
  }
  if (anyDuplicated(features)) {
    warning("`features` contains duplicate variable names.")
  }
  if (!is.null(count_vars)) {
    missing_count_vars<-setdiff(count_vars,names(data))
    if (length(missing_count_vars) > 0) {
      stop("Columns not found in `count_vars`: ",paste(missing_count_vars,collapse=", "))
    }
  }
  
  #Helper to get summary data from variables
  make_feature_summary<-function(data,features,count_vars=NULL,
                                 count_labels=NULL,sum_value=NULL){
    
    #Summary stats
    summary_df<-tibble::tibble(
      variable = features,
      total_n  = nrow(data),
      missing  = purrr::map_int(features,~sum(is.na(data[[.x]]))))
    
    #Add optional columns if passed
    if (!is.null(count_vars)) {
      
      #Apply labels if supplied
      if (!is.null(count_labels)) {
        
        #Safety check
        if (length(count_labels) != length(count_vars)) {
          stop("`count_labels` must have the same length as `count_vars`.")
        }
        
        #Add names to variables
        names(count_vars)<-count_labels
      } else if (is.null(names(count_vars))) {
        names(count_vars)<-count_vars
      }
      
      #Count data from variables
      count_df<-purrr::imap_dfc(count_vars,function(var,label) {
        col<-data[[var]]
        count<-if (is.logical(col)) {
          sum(col,na.rm = TRUE)
        } else if (is.numeric(col)) {
          sum(col==1,na.rm=TRUE)
        } else {
          if (is.null(sum_value))
            stop(sprintf("Column '%s' requires `sum_value`.",var))
          sum(col==sum_value,na.rm=TRUE)
        }
        tibble::tibble(!!label := rep(count,length(features)))
      })
      summary_df<-dplyr::bind_cols(summary_df,count_df)
    }
    return(summary_df)
  }
  
  #Helper to adjust CI limits for models that do not converge
  clean_uv_bounds<-function(tbl,max_CI,min_CI) {
    tbl$table_body<-tbl$table_body |>
      mutate(
        #Replace upper infinity values with a practical high cap 
        conf.high=ifelse(is.infinite(conf.high) | conf.high > max_CI,max_CI,conf.high),
        # Replace 0 values with a small decimal so log-scales do not crash
        conf.low=ifelse(conf.low <= 0.1,min_CI,conf.low)
      )
    return(tbl)
  }
  
  #Call feature summary and create objects with column names and labels for gtsummary
  feature_summary<-make_feature_summary(
    data         = data,
    features     = features,
    count_vars   = count_vars,
    count_labels = count_labels,
    sum_value    = sum_value
  )
  summary_cols  <-names(feature_summary)[names(feature_summary) != "variable"]
  header_labels <-c(
    total_n = "**Total**",missing = "**Missing**", #fixed names
    setNames(paste0("**",setdiff(summary_cols,c("total_n","missing")),"**"),
             setdiff(summary_cols,c("total_n","missing"))))
  
  #Temporarily assign to global environment to avoid issue with tbl_uvregression
  #.uv_data<<-data
  .uv_surv<<-Surv(data[[time_var]],data[[event_var]])
  on.exit({ rm(list=c(".uv_surv"),envir=.GlobalEnv) },add=TRUE)
  
  #Generate UV table
  uv_tab<-tbl_uvregression(
    data         = data,
    y            = .uv_surv,
    include      = all_of(features),
    label        = if(!is.null(labels)) stats::setNames(as.list(labels),features) else NULL,
    method       = coxph,
    exponentiate = TRUE,
    hide_n = TRUE,
    add_estimate_to_reference_rows = TRUE,
    pvalue_fun = ~style_pvalue(.x,digits=2)
  ) |> bold_p() |> bold_labels() |>  add_n() |> add_nevent() |> 
    modify_table_body( 
      ~ dplyr::left_join(.x,feature_summary,by="variable") |> 
        dplyr::relocate(dplyr::all_of(summary_cols),stat_n,.before=stat_nevent)
    ) |>
    modify_header(!!!header_labels,stat_n="**Model**",stat_nevent="**Event**") |>
    modify_column_unhide(columns=dplyr::all_of(summary_cols))
  
  if(adjust_CI)  {
    uv_tab<-clean_uv_bounds(uv_tab,max_CI=max_CI,min_CI=min_CI) 
    return(uv_tab)
  } else {
    return(uv_tab)
  }
}