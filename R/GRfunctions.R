.GRcalculate = function(inputData, groupingVariables, cap = FALSE, case = "A",
                        initial_count){
  # declaring values NULL to avoid note on package check
  cell_count = NULL
  cell_count__time0 = NULL
  cell_count__ctrl = NULL
  duration = NULL
  division_time = NULL
  if(initial_count) {
    log2nn = with(inputData, log2(cell_count/cell_count__time0))
    log2nn_ctrl = with(inputData, log2(cell_count__ctrl/cell_count__time0))
    GR = 2^(log2nn/log2nn_ctrl) - 1
  } else {
    log2_rel = with(inputData, log2(cell_count/cell_count__ctrl))
    log2nn_ctrl = with(inputData, treatment_duration/division_time)
    GR = 2^(1 + log2_rel/log2nn_ctrl) - 1
  }
  rel_cell_count = with(inputData, cell_count/cell_count__ctrl)
  input_edited = inputData
  input_edited$log10_concentration = log10(input_edited$concentration)
  input_edited$GRvalue = GR
  input_edited$rel_cell_count = rel_cell_count
  input_edited$ctrl_cell_doublings = log2nn_ctrl
  input_edited$treated_cell_doublings = log2nn
  tmp<-input_edited[,groupingVariables, drop = FALSE]
  experimentNew = (apply(tmp,1, function(x) (paste(x,collapse=" "))))
  if(cap == TRUE) {
    input_edited$GRvalue[input_edited$GRvalue > 1] = 1
    input_edited$rel_cell_count[input_edited$rel_cell_count > 1] = 1
  }
  if(length(groupingVariables) > 0) {
    input_edited$experiment = as.factor(experimentNew)
  } else {
    input_edited$experiment = as.factor("All Data")
  }
  return(as.data.frame(input_edited))
}

.GRlogisticFit = function(inputData, groupingVariables, force = FALSE,
                          cap = FALSE) {
  # if(length(groupingVariables) > 0) {
  #   metadata = matrix(data = NA, ncol = length(groupingVariables),
  #                     nrow = length(experiments))
  #   metadata = as.data.frame(metadata)
  #   colnames(metadata) = groupingVariables
  # } else {
  #   metadata = NULL
  # }

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

.GRlogistic_3u = function(c, GRinf, GEC50, h_GR){
  GRinf + (1 - GRinf)/(1 + (c/GEC50)^h_GR)
}

.rel_cell_logistic_3u = function(c, Einf, EC50, h){
  Einf + (1 - Einf)/(1 + (c/EC50)^h)
}

.trim_mean = function(x, percent) {
  x = x[!is.na(x)]
  n = length(x)
  k = n*(percent/100)/2
  # round down if k is half an integer
  if(round(k) != k & round(k*2) == k*2) {
    lo = floor(k) + 1
    hi = n - lo + 1
  } else {
    lo = round(k) + 1
    hi = n - lo + 1
  }
  x = sort(x)[lo:hi]
  return(mean(x))
}

.check = function(inputData, case) {
  message = NULL # an error message, if applicable
  initial_count = TRUE # a logical for whether initial cell count is provided
  input_cols = colnames(inputData)
  caseA = c('concentration', 'cell_count', 'cell_count__ctrl',
            'cell_count__time0')
  caseA_div_time = c('concentration', 'cell_count','cell_count__ctrl',
                     'treatment_duration','division_time')
  if(case == "A") {
    col_check = caseA %in% input_cols
    col_check2 = caseA_div_time %in% input_cols
    # check for correct input columns
    if(sum(col_check) != 4 & sum(col_check2) != 5) {
        message = "There must be columns named 'concentration', 'cell_count',
          'cell_count__ctrl', and 'cell_count__time0' in inputData. If 
           initial cell count (cell_count__time0) is not available, the assay 
        duration and division time of cells can be used instead in columns 
        labeled 'treatment_duration' and 'division_time'"
        return(list(message, initial_count))
    }
    num_cols = intersect(input_cols, union(caseA, caseA_div_time))
    num_cols_data = inputData[,num_cols]
    num_cols_test = unlist(lapply(num_cols_data, is.numeric))
    # check that columns are numeric
    if(sum(!num_cols_test) > 0) {
      non_numeric_cols = toString(names(which(!num_cols_test)))
      message = paste("The following columns need to be numeric: ", non_numeric_cols)
      return(list(message, initial_count))
    }
    cond1 = 'cell_count__time0' %in% colnames(inputData)
    cond2 = length(intersect(colnames(inputData), 
                   c('treatment_duration','division_time'))) == 2
    if(cond1) {
      initial_count = TRUE
      if(cond2) {
        warning("Initial cell count given, ignoring columns 'treatment_duration' and 
                'division_time' for calculation of GR values.")
      }
    } else {
      if(!cond2) {
        message = "Need initial cell count or treatment_duration and division
        time for control cells."
      }
      initial_count = FALSE
    }
  }
  if(case == "C") {
    # check for correct input columns
    if(length(intersect(colnames(inputData), c('concentration', 'cell_count',
                                               'time'))) != 3) {
      message = "There must be columns named 'concentration', 'cell_count',
           and 'time' in inputData"
    } else {
      # check for time 0 cell counts
      if(sum(inputData$time == 0) == 0) {
        initial_count = FALSE
        if(length(intersect(colnames(inputData), c('treatment_duration',
                                                   'division_time'))) != 2) {
          message = "Need initial cell count or treatment_duration and division 
          time for control cells."
        }
      } else {
        if(length(intersect(colnames(inputData), c('treatment_duration',
                                                   'division_time'))) == 2) {
          message = "You have provided both time 0 cell counts and division 
          times. Please provide one or the other."
        }
      }
    }
  }
  return(list(message, initial_count))
}

.convert = function(inputData, case, initial_count) {
  if(case == "A") {
      return(inputData)
  } else if(case == "C") {
    delete_cols = which(colnames(inputData) %in% c('concentration',
                                                   'cell_count'))
    keys = colnames(inputData)[-delete_cols]
    time0 = inputData[inputData$time == 0, c(keys, 'cell_count')]
    ctrl = inputData[inputData$concentration == 0 & inputData$time > 0,
                     c(keys, 'cell_count')]
    data = inputData[inputData$concentration != 0 & inputData$time > 0, ]
    time0_keys = NULL
    ctrl_keys = NULL
    for(i in 1:length(keys)) {
      time0_keys[i] = length(intersect(time0[[ keys[i] ]],
                                         data[[ keys[i] ]])) > 0
      ctrl_keys[i] = length(intersect(ctrl[[ keys[i] ]],
                                      data[[ keys[i] ]])) > 0
    }
    ctrl_keys = keys[ctrl_keys]
    
    time0_keys = keys[time0_keys]

    temp = ctrl[, ctrl_keys]
    ctrl$key = apply(temp, 1, function(x) paste(x, collapse = ' '))

    temp = time0[, time0_keys]
    time0$key = apply(temp, 1, function(x) paste(x, collapse = ' '))

    temp = data[, ctrl_keys]
    data$key_ctrl = apply(temp, 1, function(x) paste(x, collapse = ' '))

    temp = data[, time0_keys]
    data$key_time0 = apply(temp, 1, function(x) paste(x, collapse = ' '))

    data$cell_count__ctrl = NA
    data$cell_count__time0 = NA
    for(key in unique(ctrl$key)) {
      trimmed_mean = .trim_mean(ctrl[ctrl$key == key,]$cell_count, 50)
      data[data$key_ctrl == key, 'cell_count__ctrl'] = trimmed_mean
    }

    for(key in unique(time0$key)) {
      trimmed_mean = .trim_mean(time0[time0$key == key,]$cell_count, 50)
      data[data$key_time0 == key, 'cell_count__time0'] = trimmed_mean
    }

    delete_cols = which(colnames(data) %in% c('key_ctrl', 'key_time0'))
    data = data[, -delete_cols]

    if(!initial_count) { data$cell_count__time0 = NULL }
    data = as.data.frame(data)
    row.names(data) = 1:dim(data)[1]
    inputData = data
    return(inputData)
  }
}

#' Extract GR parameters from a dataset
#'
#' This function takes in a dataset with information about concentration,
#' cell counts over time, and additional grouping variables for a dose-response
#' assay and calculates growth-rate inhibition (GR) metrics as well as 
#' traditional metrics (IC50, Emax, etc.) for each experiment
#' in the dataset. The data must be in a specific format: either that specified
#' by case "A" or case "C" described in the details below.
#'
#' @param inputData a data table in one of the specified formats (Case A or
#' Case C). See details below for description. See \code{data(inputCaseA)} or
#' \code{data(inputCaseC)} for example input data frames. See help files for
#' \code{\link{inputCaseA}} and \code{\link{inputCaseC}} for description of
#' these examples.
#' @param groupingVariables a vector of column names from inputData. All of the
#' columns in inputData except for those identified here will be averaged over.
#' @param case either "A" or "C", indicating the format of the input data. See
#' below for descriptions of these formats.
#' @param force a logical value indicating whether to attempt to "force" a
#' sigmoidal fit, i.e. whether to allow fits with F-test p-values greater
#' than .05
#' @param cap a logical value indicating whether to cap GR values (or 
#' relative cell counts) at 1. If true, all values greater than 1 will 
#' be set to 1.
#' @return A SummarizedExperiment object containing GR metrics 
#' (GR50, GRmax, etc.) and traditional metrics (IC50, Emax, etc.) 
#' as well as goodness of fit measures is returned. The
#' object also contains, in its metadata, a table of the original data
#' converted to the style of "Case A" (with calculated GR values and relative 
#' cell counts for each row) and a vector of the grouping variables used for 
#' the calculation.
#' @author Nicholas Clark
#' @details
#' Calculation of GR values is performed by the function \code{.GRcalculate}
#' according to the "Online Methods" section of Hafner and Niepel et al.
#' (2016, \url{http://dx.doi.org/10.1038/nmeth.3853}).
#'
#' The fitting of the logistic curve is performed by the \code{.GRlogisticFit}
#' function, which calls the \code{drm} function from the \code{drc} package
#' to solve for the curve parameters. The GR curve fit function is
#' given by f(c) = GRinf + (1 - GRinf)/(1 + (c/GEC50)^h_GR) where c is
#' concentration. The fit is performed under following constraints: h_GR 
#' in [.1, 5], GRinf in [-1, 1], and GEC50 in [min(c)*1e-2, max(c)*1e2] (c is
#' concentration). The initial conditions for the fitting algorithm are h_GR 
#' = 2, GRinf = 0.1 and GEC50 = median(c). The fitting of the 
#' traditional dose response curve is done using the same formula, 
#' replacing GRinf with Einf, GEC50 with EC50, and h_GR with h. The fit is 
#' performed on the relative cell counts instead of GR values. Also, since the 
#' traditional dose response curve is bounded between 0 and 1 whereas the 
#' GR dose response curve is bounded between -1 and 1, we restrict Einf to 
#' the range [0, 1].
#'
#' The parameters of the GR dose response curves (and traditional dose 
#' response curves) for each experiment are fitted separately. An
#' F-test is used to compare the sigmoidal fit to a flat line fit. If the
#' p-value of the F-test is less than .05, the sigmoidal fit is accepted. If
#' the p-value is greater than or equal to .05, a flat horizontal line fit is
#' given, with y equal to the mean of the GR values (or relative cell counts 
#' in the case of the traditional dose response curve). For each flat fit, 
#' GEC50 (or EC50) is set to 0, h_GR (or h) is set to 0.01, GRinf (or Einf) is
#' set to the y value of the flat fit, and GR50 (or IC50) is set to +/-Inf 
#' depending on whether GRinf (or Einf) is greater or less than .5. 
#'
#' The mandatory columns for inputData for Case "A" are the following as
#' well as other grouping columns.
#'
#' 1. concentration - column with concentration values (not log transformed)
#' of the perturbagen on which dose-response curves will be evaluated
#'
#' 2. cell_count - column with the measure of cell number (or a surrogate of
#' cell number) after treatment
#'
#' 3. cell_count__time0 - column with initial (Time 0) cell counts - the
#' measure of cell number in untreated wells grown in parallel until the
#' time of treatment
#'
#' 4. cell_count__ctrl - column with the Control cell count: the measure of
#' cell number in control (e.g. untreated or DMSO-treated) wells from the
#' same plate
#'
#' All other columns will be treated as additional keys on which the data
#' will be grouped (e.g. cell_line, drug, time, replicate)
#'
#' The mandatory columns for inputData for Case "C" are the following as
#' well as other grouping columns.
#'
#' 1. concentration - column with concentration values (not log transformed)
#' of the perturbagen on which dose-response curves will be evaluated
#'
#' 2. cell_count - column with the measure of cell number (or a surrogate of
#' cell number)
#'
#' 3. time - column with the time at which a cell count is observed
#'
#' All other columns will be treated as additional keys on which the data
#' will be grouped (e.g. cell_line, drug, replicate)
#' 
#' GR values and dose-response curves/metrics can also be computed using
#' division times for (untreated) cell lines in the place of time zero cell
#' counts, using the first formula in the Supplement of Hafner et al. (2017,
#' \url{http://dx.doi.org/10.1038/nbt.3882}).
#' 
#' To use division rate instead of initial cell count,
#' inputData should not have any initial cell counts (i.e. For Case "A", no 
#' "cell_count__time0" column. For Case "C", no values of 0 in the "time" 
#' column) and should instead have two columns "treatment_duration" and 
#' "division_time".
#' 
#' In the first column, "treatment duration", one should have the duration of 
#' the assay between time of treatment and the final cell counts (e.g. 72 for 
#' hours in a typical 3-day assay). In the second column, "division_time", one 
#' should have the time it takes for one cell doubling to occur in each
#' (untreated) cell line used under the conditions of the experiment. These 
#' two columns must contain numbers (no units), but need to refer to the same
#' units (e.g. hours). In most cases, all experiments of a particular cell 
#' line would have the same "division_time", however if the division rate of 
#' untreated cells varied on another parameter, for example seeding density,
#' it would be appropriate to measure and input division times based on 
#' cell line/seeding density pairs.
#' 
#' 
#' @note
#' To see the underlying code, use (\code{getAnywhere(.GRlogistic_3u)}), 
#' (\code{getAnywhere(.rel_cell_logistic_3u)}),
#' (\code{getAnywhere(.GRcalculate)}), and (\code{getAnywhere(.GRlogisticFit)})
#' @seealso See \code{\link{drm}} for the general logistic fit function that
#' solves for the parameters GRinf, GEC50, and h_GR. See
#' \code{\link{drmc}} for
#' options of this function. Use the functions \code{\link{GRdrawDRC}},
#' \code{\link{GRbox}}, and \code{\link{GRscatter}} to create visualizations
#' using the output from this function. For online GR calculator and browser,
#' see \url{http://www.grcalculator.org}.
#' @references Hafner, M., Niepel, M., Chung, M., and Sorger, P.K.,
#' "Growth Rate Inhibition Metrics Correct For Confounders In Measuring
#' Sensitivity To Cancer Drugs". \emph{Nature Methods} 13.6 (2016): 521-527.
#' \url{http://dx.doi.org/10.1038/nmeth.3853}
#' @references Hafner, M., Niepel, M., Sorger, P.K.,
#' "Alternative drug sensitivity metrics improve preclinical cancer
#' pharmacogenomics". \emph{Nature Biotechnology} 35.6 (2017): 500-502.
#' \url{http://dx.doi.org/10.1038/nbt.3882}
#'
#' @examples
#' # Load Case A (example 1) input
#' data("inputCaseA")
#' head(inputCaseA)
#' # Run GRfit function with case = "A"
#' output1 = GRfit(inputData = inputCaseA,
#' groupingVariables = c('cell_line','agent', 'perturbation','replicate',
#' 'time'))
#' # Overview of SummarizedExperiment output data
#' output1
#' \dontrun{
#' # View GR metrics table
#' View(GRgetMetrics(output1))
#' # View descriptions of each metric (or goodness of fit measure)
#' View(GRgetDefs(output1))
#' # View table of original data (converted to style of Case A) with GR values
#' # and relative cell counts
#' View(GRgetValues(output1))
#' # View vector of grouping variables used for calculation
#' GRgetGroupVars(output1)
#' }
#' # Load Case C (example 4) input
#' # Same data, different format
#' data("inputCaseC")
#' head(inputCaseC)
#' output4 = GRfit(inputData = inputCaseC,
#' groupingVariables = c('cell_line','agent', 'perturbation','replicate',
#' 'time'),
#' case = "C")
#' # Extract data tables and export to .tsv or .csv
#' \dontrun{
#' # Write GR metrics parameter table to tab-separated text file
#' write.table(GRgetMetrics(output1), file = "filename.tsv", quote = FALSE,
#' sep = "\t", row.names = FALSE)
#' # Write original data plus GR values to comma-separated file
#' write.table(GRgetValues(output1), file = "filename.csv", quote = FALSE,
#' sep = ",", row.names = FALSE)
#' }
#' @export

GRfit = function(inputData, groupingVariables, case = "A",
                 force = FALSE, cap = FALSE) {
  if('experiment' %in% colnames(inputData)) {
    stop("Change name of 'experiment' column.")
  }
  input_check = .check(inputData, case)
  message = input_check[[1]]
  initial_count = input_check[[2]]
  if(!is.null(message)) stop(message)
  inputData = .convert(inputData, case, initial_count)
  gr_table = .GRcalculate(inputData, groupingVariables, cap, case,
                          initial_count)
  GRlogfit = .GRlogisticFit(gr_table, groupingVariables, force, cap)
  #parameter_table = GRlogfit$parameters
  #GR_drc_list = GRlogfit$GR_drc_list
  #trad_drc_list = GRlogfit$trad_drc_list
  ##### temporary - return parameters ####
  output = list(
    assays = GRlogfit,
    #colData = colData,
    #rowData = rowData, 
    metadata = list(gr_table = gr_table, 
                    groupingVariables = groupingVariables))
  return(output)
  ########
  
  colData = parameter_table[ ,c(groupingVariables, 'fit_GR', 'fit_rel_cell',
                                'experiment', 'concentration_points')]
  rownames(colData) = colData$experiment
  colData = S4Vectors::DataFrame(colData)
  
  Metric = c('ctrl_cell_doublings','GR50','log10_GR50','GRmax','GR_AOC','GEC50','log10_GEC50',
             'GRinf','h_GR','r2_GR','pval_GR','flat_fit_GR', 
              'IC50', 'log10_IC50','Emax', 'AUC', 'EC50', 'log10_EC50','Einf', 'h', 
             "Einf_1", "log10_GEC50_1", "h_1", "Einf_2", "log10_GEC50_2", "h_2",
              'r2_rel_cell', 'pval_rel_cell', 'flat_fit_rel_cell')
  assays = parameter_table[ , Metric]
  rownames(assays) = parameter_table$experiment
  assays = t(assays)

  Description = c(
    "The number of cell doublings in the control population during the assay",
    "The concentration at which GR(c) = 0.5",
    "log10 value of GR50",
    "The maximal effect of the drug (minimal GR value)",
    "The 'Area Over the Curve' - The area between the line GR = 1 and the curve, similar to traditional AUC",
    "The concentration at half-maximal effect (growth rate normalized)",
    "The asymptotic effect of the drug (growth rate normalized)",
    "The Hill coefficient of the fitted (GR) curve, which reflects how steep the (GR) dose response curve is",
    "log10 value of GEC50",
    "The coefficient of determination - essentially how well the (GR) curve fits to the data points",
    "The p-value of the F-test comparing the fit of the (GR) curve to a horizontal line fit",
    "For data that doesn't significantly fit better to a curve than a horizontal line fit, the y value (GR) of the flat line", 
    "The concentration at which relative cell count = 0.5",
    "log10 value of IC50",
    "The maximal effect of the drug (minimal relative cell count value)",
    "The 'Area Under the Curve' - The area below the fitted (traditional) dose response curve",
    "The concentration at half-maximal effect (not growth rate normalized)",
    "The asymptotic effect of the drug (not growth rate normalized)",
    "The Hill coefficient of the fitted (traditional) dose response curve, which reflects how steep the (traditional) dose response curve is",
    rep("", 6),
    "log10 value of EC50",
    "The coefficient of determination - essentially how well the (traditional) curve fits to the data points",
    "The p-value of the F-test comparing the fit of the (traditional) curve to a horizontal line fit",
    "For data that doesn't significantly fit better to a curve than a horizontal line fit, the y value (relative cell count) of the flat line"
                  )
  rowData = cbind(Metric, Description)
  rownames(rowData) = Metric
  rowData = S4Vectors::DataFrame(rowData)
  rowData$Metric = as.character(rowData$Metric)
  rowData$Description = as.character(rowData$Description)

  output = SummarizedExperiment::SummarizedExperiment(assays = assays,
                                                      colData = colData,
            rowData = rowData, metadata = list(gr_table = gr_table, groupingVariables = groupingVariables,
                                               GR_drc_list = GR_drc_list, trad_drc_list = trad_drc_list))
  return(output)
}

