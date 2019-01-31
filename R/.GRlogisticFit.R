.GRlogisticFit = function(inputData, groupingVariables, force = FALSE,
                          cap = FALSE, case) {
  ## Define the biphasic dose-response function
  ## x is concentration, p is the parameter vector
  opfct_bi = function(x, p) {
    term1 = 1 + (p[1] + (1 - p[1])/(1 + (x / (10^p[2])) ^ p[3]))
    term2 = 1 + (p[4] + (1 - p[4])/(1 + (x / (10^p[5])) ^ p[6]))
    2^( 0.5*( log2(term1) + log2(term2) ) ) - 1
  }
  ## Define the sigmoidal (or logistic) dose-response function
  opfct_sig = function(x, p) {
    p[1] + (1 - p[1])/(1 + (x / (10^p[2])) ^ p[3])
  }
  ## Define the sigmoidal (or logistic) dose-response function for GR_d
  opfct_sig_GR_d = function(x, p) {
    p[1] + (0 - p[1])/(1 + (x / (10^p[2])) ^ p[3])
  }
  if(case == "static_vs_toxic") {
    ### should "time" be added to the grouping variables??
    grp = dplyr::syms(groupingVariables)
    data_grp = inputData %>% dplyr::group_by(experiment, !!!grp)
    data_grp_summ = data_grp %>% 
      dplyr::filter(concentration > 0) %>%
      dplyr::summarise(
        GR_s_mean = mean(GR_s, na.rm = TRUE),
        GR_d_mean = mean(GR_d, na.rm = TRUE),
        #rel_cell_mean = mean(rel_cell_count, na.rm = TRUE),
        #ctrl_cell_doublings = mean(ctrl_cell_doublings, na.rm = TRUE),
        #treated_cell_doublings = mean(treated_cell_doublings, na.rm = TRUE),
        concentration_points = length(unique(concentration)),
        conc = list(unique(concentration)),
        cc = list( log10(unique(concentration)) ),
        GR_s = list ( GR_s ),
        GR_d = list ( GR_d ),
        #rel_cell_count = list ( rel_cell_count ),
        concentration = list ( concentration )
      )
    ## get RSS of flat fit for GR and relative cell count curve
    data_grp_summ %<>% dplyr::mutate(
      RSS1_GR_s = sum( (unlist(GR_s) - GR_s_mean)^2, na.rm = TRUE ),
      RSS1_GR_d = sum( (unlist(GR_d) - GR_d_mean)^2, na.rm = TRUE )
      )
    ## calculate GRmax (or Emax) and GR_AOC (or AUC) directly from unfitted data points
    data_grp_conc = inputData %>% 
      dplyr::filter(concentration > 0) %>%
      dplyr::group_by(experiment, !!!grp, log10_concentration) %>%
      dplyr::summarise(
        GR_s_mean = mean(GR_s, na.rm = TRUE),
        GR_d_mean = mean(GR_s, na.rm = TRUE),
        cc = unique(log10_concentration)
      ) %>%
      ## get maximum (and second largest) concentration for each group
      dplyr::mutate(max_conc = list(sort(cc, decreasing = T)[1:2]) )
    data_grp_conc_max = data_grp_conc %>% 
      ## look at only two highest concentrations
      dplyr::filter( cc %in%  unlist(max_conc) ) %>%
      ## GRmax = min of avg. GR value at the two highest concentrations
      dplyr::summarise(
        GR_s_max = min(GR_s_mean, na.rm = T),
        GR_d_max = min(GR_d_mean, na.rm = T)
      ) %>% dplyr::ungroup()
    data_grp_conc_all = data_grp_conc %>% dplyr::arrange(cc) %>%
      ### Not sure how we want to compute AUC/AOC for GR static and toxic yet
      # dplyr::summarise(
      #   GR_s_AOC = caTools::trapz(x = cc, y = 1 - GR_s_mean),
      #   GR_d_AOC = caTools::trapz(x = cc, y = - GR_d_mean) ) %>%
      # dplyr::ungroup()
      dplyr::summarise(
        GR_s_AOC = NA,
        GR_d_AOC = NA ) %>%
      dplyr::ungroup()
    
    data_grp_conc = suppressMessages(
      dplyr::full_join(data_grp_conc_all, data_grp_conc_max, 
                       by = c("experiment", groupingVariables) )
      )
    data_grp_summ %<>% dplyr::group_by(experiment, !!!grp) %>%
      dplyr::mutate(
        p_sig_GR_s = list( tibble::tribble(
          ~parameter,                 ~lower,     ~prior,    ~upper,
          "GRinf",                       0,        0.5,         1, 
          "log10_GEC50", min(unlist(cc))-2, median(unlist(cc)), max(unlist(cc))+2,
          "h_GR",                         0.1,          2,         5
        ) ),
        p_sig_GR_d = list( tibble::tribble(
          ~parameter,                 ~lower,     ~prior,    ~upper,
          "GRinf",                       -1,        -0.5,         0, 
          "log10_GEC50", min(unlist(cc))-2, median(unlist(cc)), max(unlist(cc))+2,
          "h_GR",                         0.1,          2,         5
        ) )
      )
    sum_square_sig_GR_s = function(x, y, p) {sum((y - opfct_sig(x, p))^2)}
    sum_square_sig_GR_d = function(x, y, p) {sum((y - opfct_sig_GR_d(x, p))^2)}
    ## Fit GR_s sigmoidal curve
    data_grp_summ$sig_fit_GR_s = lapply(1:dim(data_grp_summ)[1], function(i) {
      param_df = data_grp_summ$p_sig_GR_s[[i]]
      xx = data_grp_summ$concentration[[i]]
      yy = data_grp_summ$GR_s[[i]]
      fit = suppressMessages(try(optim(par = param_df$prior, 
                      function(p, x, y) sum_square_sig_GR_s(x = xx, y = yy, p = p),
                      hessian = TRUE, method = "L-BFGS-B",
                      lower = param_df$lower, upper = param_df$upper)))
      fit$parameters = param_df$parameter
      fit$lower = param_df$lower
      fit$upper = param_df$upper
      fit$prior = param_df$prior
      return(fit)
    })
    ## Fit GR_s sigmoidal curve
    data_grp_summ$sig_fit_GR_d = lapply(1:dim(data_grp_summ)[1], function(i) {
      param_df = data_grp_summ$p_sig_GR_d[[i]]
      xx = data_grp_summ$concentration[[i]]
      yy = data_grp_summ$GR_d[[i]]
      fit = suppressMessages(try(optim(par = param_df$prior, 
                      function(p, x, y) sum_square_sig_GR_d(x = xx, y = yy, p = p),
                      hessian = TRUE, method = "L-BFGS-B",
                      lower = param_df$lower, upper = param_df$upper)))
      fit$parameters = param_df$parameter
      fit$lower = param_df$lower
      fit$upper = param_df$upper
      fit$prior = param_df$prior
      return(fit)
    })
    fit_types = list(static = NULL, toxic = NULL)
    parameters = list(GR = fit_types)
    
    for(x in c("sig_fit_GR_s", "sig_fit_GR_d")) {
      params = data_grp_summ[[x]][[1]]$parameters
      df = data_grp_summ %>% 
        dplyr::select(experiment, !!!grp, #ctrl_cell_doublings, 
          #treated_cell_doublings, 
          concentration_points) %>% 
        dplyr::ungroup() %>%
        dplyr::mutate_if(is.factor, as.character)
      ## get fitted curve parameters
      vals = sapply(data_grp_summ[[x]], function(y) { 
        if(!class(y) == "try-error") { y$par } else { rep(NA, length(params)) }
      }) %>% t() %>% as.data.frame() %>%
        magrittr::set_colnames(params) %>%
        dplyr::mutate_if(is.factor, function(z) as.numeric(as.character(z)))
      df = cbind(df, vals)
      ## Get GRmax (Emax) and GR_AOC (AUC)
      if(grepl("_GR_s", x)) { 
        df =  suppressMessages(
          dplyr::left_join(df, 
            data_grp_conc %>% dplyr::select(experiment, !!!grp, GR_s_max, GR_s_AOC),
              by = c("experiment", groupingVariables) ))
      } else if(grepl("_GR_d", x)) {
        df =  suppressMessages(
          dplyr::left_join(df, 
            data_grp_conc %>% dplyr::select(experiment, !!!grp, GR_d_max, GR_d_AOC),
              by = c("experiment", groupingVariables) ))
      }
      ## Get RSS1 
      if(grepl("_GR_s", x)) { df$RSS1 = data_grp_summ$RSS1_GR_s }
      if(grepl("_GR_d", x)) { df$RSS1 = data_grp_summ$RSS1_GR_d }
      
      ## Get RSS2
      df$RSS2 = sapply(data_grp_summ[[x]], function(y) {
        if("value" %in% names(y) && is.numeric(y$value)) { y$value } else { NA }
      })
      
      Npara_flat = 1
      Npara = length(params)
      df$df1 = Npara_flat
      df$df2 = Npara
      if(grepl("_GR_s", x)) {
        df$n =  sapply(data_grp_summ$GR_s, function(y) return(length(na.omit(y))) )
      }
      if(grepl("_GR_d", x)) {
        df$n =  sapply(data_grp_summ$GR_d, function(y) return(length(na.omit(y))) )
      }
      ## note: f_value same as before, just expressed differently
      df$f_value = with(df, ( (RSS1 - RSS2)/(df2 - df1) )/(RSS2/(n - df2) ) )
      df$f_pval = with(df, stats::pf(f_value, df1, df2, lower.tail = FALSE) )
      ## note: RSS1 = residual sum of squares of flat fit = total sum of squares
      df$R_square = with(df, 1 - RSS2/RSS1 )
      
      #pcutoff = ifelse(force == FALSE, .05 , 1)
      pcutoff = 0.05
      # Flat or sigmoid fit for GR curve
      df %<>% dplyr::mutate(fit = ifelse(f_pval > pcutoff | is.na(f_pval), "flat", "curve" ))
      
      if(grepl("_GR_s", x)) { df$flat = data_grp_summ$GR_s_mean }
      if(grepl("_GR_d", x)) { df$flat = data_grp_summ$GR_d_mean }

      if(x == "sig_fit_GR_s") { parameters$GR$static = df }
      if(x == "sig_fit_GR_d") { parameters$GR$toxic = df }
    }
    return(parameters)
  }
  grp = dplyr::syms(groupingVariables)
  data_grp = inputData %>% dplyr::group_by(experiment, !!!grp)
  data_grp_summ = data_grp %>% dplyr::summarise(
    ## mean GRvalue by group (all concentrations)
    GR_mean = mean(GRvalue, na.rm = TRUE),
    rel_cell_mean = mean(rel_cell_count, na.rm = TRUE),
    ctrl_cell_doublings = mean(ctrl_cell_doublings, na.rm = TRUE),
    treated_cell_doublings = mean(treated_cell_doublings, na.rm = TRUE),
    concentration_points = length(unique(concentration)),
    conc = list(unique(concentration)),
    cc = list( log10(unique(concentration)) ),
    GRvalue = list ( GRvalue ),
    rel_cell_count = list ( rel_cell_count ),
    concentration = list ( concentration )
  )
  ## get RSS of flat fit for GR and relative cell count curve
  data_grp_summ %<>% dplyr::mutate(
    RSS1_GR = sum( (unlist(GRvalue) - GR_mean)^2, na.rm = TRUE ),
    RSS1_rel_cell = sum( (unlist(rel_cell_count) - rel_cell_mean)^2, na.rm = TRUE )
  )
  ## calculate GRmax (or Emax) and GR_AOC (or AUC) directly from unfitted data points
  data_grp_conc = inputData %>% dplyr::group_by(experiment, !!!grp, log10_concentration) %>%
    dplyr::summarise(
      # mean GRvalue by concentration in each group
      GR_mean = mean(GRvalue, na.rm = TRUE),
      rel_cell_mean = mean(rel_cell_count, na.rm = TRUE),
      # GRvalue = list ( GRvalue ),
      # rel_cell_count = list ( rel_cell_count ),
      cc = unique(log10_concentration)
    ) %>%
    ## get maximum (and second largest) concentration for each group
    dplyr::mutate(max_conc = list(sort(cc, decreasing = T)[1:2]) )
  data_grp_conc_max = data_grp_conc %>% 
    ## look at only two highest concentrations
    dplyr::filter( cc %in%  unlist(max_conc) ) %>%
    ## GRmax (Emax) = min of avg. GR value (avg. relative cell count) at the two highest concentrations
    dplyr::summarise(
      GRmax = min(GR_mean, na.rm = T),
      Emax = min(rel_cell_mean, na.rm = T)
    ) %>% dplyr::ungroup()
  data_grp_conc_all = data_grp_conc %>% dplyr::arrange(cc) %>%
    ### GR_AOC = area OVER the GR curve
    dplyr::summarise(
      GR_AOC = caTools::trapz(x = cc, y = 1 - GR_mean),
      AUC = caTools::trapz(x = cc, y = rel_cell_mean) ) %>%
    dplyr::ungroup()
  
  data_grp_conc = suppressMessages(dplyr::full_join(data_grp_conc_all, data_grp_conc_max, 
                                                    by = c("experiment", groupingVariables) ) )
  
  data_grp_summ %<>% dplyr::group_by(experiment, !!!grp) %>%
    dplyr::mutate(
      ec50_low = log10(max(c(min(unlist(conc)) * 1e-4, 1e-7))) ,
      ec50_high = log10(min(c(max(unlist(conc)) * 1e2, 1e2))) 
    ) %>%
    dplyr::mutate(
      p_bi_GR = list ( tibble::tribble(
        ~parameter,                 ~lower,     ~prior,             ~upper,
        "GRinf_1",                     -.05,        0.1,                  1, 
        "log10_GEC50_1",          ec50_low, median(unlist(cc)),   log10(1),
        "h_GR_1",                       0.025,          2,                  5,
        "GRinf_2",                       -1,       -0.1,                0.5, 
        "log10_GEC50_2",        log10(0.3),   log10(1),          ec50_high,
        "h_GR_2",                       0.025,          2,                  5
      ) ),
      p_bi_rel = list ( tibble::tribble(
        ~parameter,                 ~lower,     ~prior,             ~upper,
        "Einf_1",                      0.5,       0.75,                  1, 
        "log10_EC50_1",          ec50_low, median(unlist(cc)),   log10(1),
        "h_1",                       0.025,          2,                  5,
        "Einf_2",                        0,       0.25,                0.5, 
        "log10_EC50_2",        log10(0.3),   log10(1),          ec50_high,
        "h_2",                       0.025,          2,                  5
      ) ),
      p_sig_GR = list( tibble::tribble(
        ~parameter,                 ~lower,     ~prior,    ~upper,
        "GRinf",                       -1,        0.1,         1, 
        "log10_GEC50", min(unlist(cc))-2, median(unlist(cc)), max(unlist(cc))+2,
        "h_GR",                         0.1,          2,         5
      ) ),
      p_sig_rel = list ( tibble::tribble(
        ~parameter,                 ~lower,     ~prior,    ~upper,
        "Einf",                        0,        0.1,         1, 
        "log10_EC50", min(unlist(cc))-2, median(unlist(cc)), max(unlist(cc))+2,
        "h",                         0.1,          2,         5
      ) )
    )
  #yexp_bi = data_exp$GRvalue
  ## Define the residual sum of squares functions
  sum_square_bi = function(x, y, p) {sum((y - opfct_bi(x, p))^2)}
  sum_square_sig = function(x, y, p) {sum((y - opfct_sig(x, p))^2)}
  #bi_controls = list(maxit = 500)
  ## Fit biphasic curve (GR)
  data_grp_summ$bi_fit_GR = lapply(1:dim(data_grp_summ)[1], function(i) {
    param_df = data_grp_summ$p_bi_GR[[i]]
    xx = data_grp_summ$concentration[[i]]
    yy = data_grp_summ$GRvalue[[i]]
    fit = try(optim(par = param_df$prior, 
                    function(p, x, y) sum_square_bi(x = xx, y = yy, p = p),
                    hessian = TRUE, method = "L-BFGS-B", 
                    lower = param_df$lower, upper = param_df$upper))
    fit$parameters = param_df$parameter
    fit$lower = param_df$lower
    fit$upper = param_df$upper
    fit$prior = param_df$prior
    return(fit)
  })
  ## Fit biphasic curve (relative cell count)
  data_grp_summ$bi_fit_rel = lapply(1:dim(data_grp_summ)[1], function(i) {
    param_df = data_grp_summ$p_bi_rel[[i]]
    xx = data_grp_summ$concentration[[i]]
    yy = data_grp_summ$rel_cell_count[[i]]
    fit = try(optim(par = param_df$prior, 
                    function(p, x, y) sum_square_bi(x = xx, y = yy, p = p),
                    hessian = TRUE, method = "L-BFGS-B", 
                    lower = param_df$lower, upper = param_df$upper))
    fit$parameters = param_df$parameter
    fit$lower = param_df$lower
    fit$upper = param_df$upper
    fit$prior = param_df$prior
    return(fit)
  })
  ## Fit low sigmoidal fit (GR)
  data_grp_summ$sig_low_fit_GR = lapply(1:dim(data_grp_summ)[1], function(i) {
    param_df = data_grp_summ$p_bi_GR[[i]][1:3,]
    ## remove "_1" and "_2" from parameter names
    param_df %<>% dplyr::mutate(parameter = gsub("_1$", "", parameter))
    xx = data_grp_summ$concentration[[i]]
    yy = data_grp_summ$GRvalue[[i]]
    fit = try(optim(par = param_df$prior, 
                    function(p, x, y) sum_square_sig(x = xx, y = yy, p = p),
                    hessian = TRUE, method = "L-BFGS-B", 
                    lower = param_df$lower, upper = param_df$upper))
    fit$parameters = param_df$parameter
    fit$lower = param_df$lower
    fit$upper = param_df$upper
    fit$prior = param_df$prior
    return(fit)
  })
  ## Fit high sigmoidal fit (GR)
  data_grp_summ$sig_high_fit_GR = lapply(1:dim(data_grp_summ)[1], function(i) {
    param_df = data_grp_summ$p_bi_GR[[i]][4:6,]
    ## remove "_1" and "_2" from parameter names
    param_df %<>% dplyr::mutate(parameter = gsub("_2$", "", parameter))
    xx = data_grp_summ$concentration[[i]]
    yy = data_grp_summ$GRvalue[[i]]
    fit = try(optim(par = param_df$prior, 
                    function(p, x, y) sum_square_sig(x = xx, y = yy, p = p),
                    hessian = TRUE, method = "L-BFGS-B", 
                    lower = param_df$lower, upper = param_df$upper))
    fit$parameters = param_df$parameter
    fit$lower = param_df$lower
    fit$upper = param_df$upper
    fit$prior = param_df$prior
    return(fit)
  })
  ## Fit low sigmoidal fit (relative cell count)
  data_grp_summ$sig_low_fit_rel = lapply(1:dim(data_grp_summ)[1], function(i) {
    param_df = data_grp_summ$p_bi_rel[[i]][1:3,]
    ## remove "_1" and "_2" from parameter names
    param_df %<>% dplyr::mutate(parameter = gsub("_1$", "", parameter))
    xx = data_grp_summ$concentration[[i]]
    yy = data_grp_summ$rel_cell_count[[i]]
    fit = try(optim(par = param_df$prior, 
                    function(p, x, y) sum_square_sig(x = xx, y = yy, p = p),
                    hessian = TRUE, method = "L-BFGS-B", 
                    lower = param_df$lower, upper = param_df$upper))
    fit$parameters = param_df$parameter
    fit$lower = param_df$lower
    fit$upper = param_df$upper
    fit$prior = param_df$prior
    return(fit)
  })
  ## Fit high sigmoidal fit (relative cell count)
  data_grp_summ$sig_high_fit_rel = lapply(1:dim(data_grp_summ)[1], function(i) {
    param_df = data_grp_summ$p_bi_rel[[i]][4:6,]
    ## remove "_1" and "_2" from parameter names
    param_df %<>% dplyr::mutate(parameter = gsub("_2$", "", parameter))
    xx = data_grp_summ$concentration[[i]]
    yy = data_grp_summ$rel_cell_count[[i]]
    fit = try(optim(par = param_df$prior, 
                    function(p, x, y) sum_square_sig(x = xx, y = yy, p = p),
                    hessian = TRUE, method = "L-BFGS-B",
                    lower = param_df$lower, upper = param_df$upper))
    fit$parameters = param_df$parameter
    fit$lower = param_df$lower
    fit$upper = param_df$upper
    fit$prior = param_df$prior
    return(fit)
  })
  
  ## Fit normal sigmoidal fit (GR)
  data_grp_summ$sig_fit_GR = lapply(1:dim(data_grp_summ)[1], function(i) {
    param_df = data_grp_summ$p_sig_GR[[i]]
    xx = data_grp_summ$concentration[[i]]
    yy = data_grp_summ$GRvalue[[i]]
    fit = try(optim(par = param_df$prior, 
                    function(p, x, y) sum_square_sig(x = xx, y = yy, p = p),
                    hessian = TRUE, method = "L-BFGS-B",
                    lower = param_df$lower, upper = param_df$upper))
    fit$parameters = param_df$parameter
    fit$lower = param_df$lower
    fit$upper = param_df$upper
    fit$prior = param_df$prior
    return(fit)
  })
  ## Fit normal sigmoidal fit (relative cell count)
  data_grp_summ$sig_fit_rel = lapply(1:dim(data_grp_summ)[1], function(i) {
    param_df = data_grp_summ$p_sig_rel[[i]]
    xx = data_grp_summ$concentration[[i]]
    yy = data_grp_summ$rel_cell_count[[i]]
    fit = try(optim(par = param_df$prior, 
                    function(p, x, y) sum_square_sig(x = xx, y = yy, p = p),
                    hessian = TRUE, method = "L-BFGS-B",
                    lower = param_df$lower, upper = param_df$upper))
    fit$parameters = param_df$parameter
    fit$lower = param_df$lower
    fit$upper = param_df$upper
    fit$prior = param_df$prior
    return(fit)
  })
  
  constraints_sig = list(normal = NULL, low = NULL, high = NULL)
  constraints_bi = list(normal = NULL)
  fit_types = list(sigmoid = constraints_sig, biphasic = constraints_bi)
  parameters = list(GR = fit_types, rel_cell = fit_types)
  
  for(x in c("sig_fit_rel", "sig_fit_GR", "sig_low_fit_rel", "sig_low_fit_GR",
             "sig_high_fit_rel", "sig_high_fit_GR", "bi_fit_rel", "bi_fit_GR")) {
    params = data_grp_summ[[x]][[1]]$parameters
    df = data_grp_summ %>% dplyr::select(experiment, !!!grp, ctrl_cell_doublings, treated_cell_doublings,
                                         concentration_points) %>% dplyr::ungroup()
    ## get fitted curve parameters
    vals = sapply(data_grp_summ[[x]], function(y) { 
      if(!class(y) == "try-error") { y$par } else { rep(NA, length(params)) }
    }) %>% t() %>% as.data.frame() %>%
      magrittr::set_colnames(params)
    df = cbind(df, vals)
    ## Get GRmax (Emax) and GR_AOC (AUC)
    if(grepl("_GR", x)) { 
      df =  suppressMessages(dplyr::full_join(df, data_grp_conc %>% dplyr::select(experiment, !!!grp, GRmax, GR_AOC),
                                              by = c("experiment", groupingVariables) ))
      if(grepl("bi_fit", x)) {
        df %<>%
          dplyr::mutate(GEC50_1 = 10^log10_GEC50_1,
                        GEC50_2 = 10^log10_GEC50_2)# %>%
        #dplyr::mutate(GR50 = GEC50*((1-GRinf)/(0.5-GRinf) - 1)^(1/h_GR) ) %>%
        #dplyr::mutate(log10_GR50 = log10(GR50) )
        # df %<>%
        #   ### notes: what about log_IC50 and log_EC50? ask Marc.
        #   dplyr::mutate(GEC50_1 = ifelse(fit == "flat", 0, GEC50_1),
        #                 GEC50_2 = ifelse(fit == "flat", 0, GEC50_2),
        #                 GRinf_1 = ifelse(fit == "flat", GRmax, GRinf_1),
        #                 GRinf_2 = ifelse(fit == "flat", GRmax, GRinf_2),
        #                 ### note: does the second ifelse statement work correctly here?
        #                 GR50 = ifelse(fit == "flat", ifelse(flat > .5, Inf, -Inf) , GR50),
        #                 h_GR_1 = ifelse(fit == "flat", 0.01, h_GR_1),
        #                 h_GR_2 = ifelse(fit == "flat", 0.01, h_GR_2)
        #   )
      } else {
        df %<>%
          dplyr::mutate(GEC50 = 10^log10_GEC50) %>%
          dplyr::mutate(GR50 = GEC50*((1-GRinf)/(0.5-GRinf) - 1)^(1/h_GR) ) %>%
          dplyr::mutate(log10_GR50 = log10(GR50) )
        ## adjust parameters for flat fits
        # df %<>%
        #   ### notes: what about log_GR50 and log_GEC50? ask Marc.
        #   dplyr::mutate(GEC50 = ifelse(fit == "flat", 0, GEC50),
        #                 GRinf = ifelse(fit == "flat", GRmax, GRinf),
        #                 ### note: does the second ifelse statement work correctly here?
        #                 GR50 = ifelse(fit == "flat", ifelse(flat > .5, Inf, -Inf) , GR50),
        #                 h_GR = ifelse(fit == "flat", 0.01, h_GR)
        #                 )
      }
    }
    if(grepl("_rel", x)) {
      df = suppressMessages(dplyr::full_join(df, data_grp_conc %>% dplyr::select(experiment, !!!grp, Emax, AUC),
                                             by = c("experiment", groupingVariables) ) )
      if(grepl("bi_fit", x)) {
        df %<>% 
          dplyr::mutate(EC50_1 = 10^log10_EC50_1,
                        EC50_2 = 10^log10_EC50_2) #%>%
        #dplyr::mutate(IC50 = EC50*((1-Einf)/(0.5-Einf) - 1)^(1/h) ) %>%
        #dplyr::mutate(log10_IC50 = log10(IC50) )
        ## adjust parameters for flat fits
        # df %<>%
        #   ### notes: what about log_IC50 and log_EC50? ask Marc.
        #   dplyr::mutate(EC50_1 = ifelse(fit == "flat", 0, EC50_1),
        #                 EC50_2 = ifelse(fit == "flat", 0, EC50_2),
        #                 Einf_1 = ifelse(fit == "flat", Emax, Einf_1),
        #                 Einf_2 = ifelse(fit == "flat", Emax, Einf_2),
        #                 ### note: does the second ifelse statement work correctly here?
        #                 IC50 = ifelse(fit == "flat", ifelse(flat > .5, Inf, -Inf) , IC50),
        #                 h_1 = ifelse(fit == "flat", 0.01, h_1),
        #                 h_2 = ifelse(fit == "flat", 0.01, h_2)
        #   )
      } else {
        df %<>% 
          dplyr::mutate(EC50 = 10^log10_EC50) %>%
          dplyr::mutate(IC50 = EC50*((1-Einf)/(0.5-Einf) - 1)^(1/h) ) %>%
          dplyr::mutate(log10_IC50 = log10(IC50) )
        ## adjust parameters for flat fits
        # df %<>%
        #   ### notes: what about log_IC50 and log_EC50? ask Marc.
        #   dplyr::mutate(EC50 = ifelse(fit == "flat", 0, EC50),
        #                 Einf = ifelse(fit == "flat", Emax, Einf),
        #                 ### note: does the second ifelse statement work correctly here?
        #                 IC50 = ifelse(fit == "flat", ifelse(flat > .5, Inf, -Inf) , IC50),
        #                 h = ifelse(fit == "flat", 0.01, h)
        #   )
      }
    }
    ## Get RSS1 
    if(grepl("_GR", x)) { df$RSS1 = data_grp_summ$RSS1_GR }
    if(grepl("_rel", x)) { df$RSS1 = data_grp_summ$RSS1_rel_cell }
    ## Get RSS2
    df$RSS2 = sapply(data_grp_summ[[x]], function(y) { 
      if(!class(y) == "try-error") { y$value } else { NA }
    })
    Npara_flat = 1
    Npara = length(params)
    df$df1 = Npara_flat
    df$df2 = Npara
    df$n =  sapply(data_grp_summ$GRvalue, function(y) return(length(na.omit(y))) )
    ## note: f_value same as before, just expressed differently
    df$f_value = with(df, ( (RSS1 - RSS2)/(df2 - df1) )/(RSS2/(n - df2) ) )
    df$f_pval = with(df, stats::pf(f_value, df1, df2, lower.tail = FALSE) )
    ## note: RSS1 = residual sum of squares of flat fit = total sum of squares
    df$R_square = with(df, 1 - RSS2/RSS1 )
    
    pcutoff = ifelse(force == FALSE, .05 , 1)
    # Flat or sigmoid fit for GR curve
    df %<>% dplyr::mutate(fit = ifelse(f_pval > pcutoff | is.na(f_pval), "flat", "curve" ))
    
    ### select columns to keep
    # df = df[,c('GR50', 'log10_GR50','GRmax','GR_AOC','GEC50', 'log10_GEC50',
    #            'GRinf','h_GR','r2_GR','pval_GR', 'fit_GR','flat_fit_GR', 
    #            'IC50', 'log10_IC50','Emax', 'AUC', 'EC50', 'log10_EC50',
    #            'Einf', 'h', 
    #            "Einf_1", "log10_GEC50_1", "h_1", "Einf_2", "log10_GEC50_2", "h_2",
    #            'r2_rel_cell', 'pval_rel_cell', 'fit_rel_cell',
    #            'flat_fit_rel_cell','experiment',
    #            'concentration_points', 'ctrl_cell_doublings')]
    
    if(grepl("_GR", x)) { df$flat = data_grp_summ$GR_mean }
    if(grepl("_rel", x)) { df$flat = data_grp_summ$rel_cell_mean }
    
    if(x == "sig_fit_GR") { parameters$GR$sigmoid$normal = df }
    if(x == "sig_low_fit_GR") { parameters$GR$sigmoid$low = df }
    if(x == "sig_high_fit_GR") { parameters$GR$sigmoid$high = df }
    if(x == "bi_fit_GR") { parameters$GR$biphasic$normal = df }
    
    if(x == "sig_fit_rel") { parameters$rel_cell$sigmoid$normal = df }
    if(x == "sig_low_fit_rel") { parameters$rel_cell$sigmoid$low = df }
    if(x == "sig_high_fit_rel") { parameters$rel_cell$sigmoid$high = df }
    if(x == "bi_fit_rel") { parameters$rel_cell$biphasic$normal = df }
  }
  
  return(parameters)
}