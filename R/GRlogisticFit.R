.GRlogisticFit = function(inputData, groupingVariables, force = FALSE,
                          cap = FALSE, case) {
  if(case == "static_vs_toxic") {
    fits = list(GR = list(sigmoid = list(normal = T, low = F, high = F, static = T, toxic = T), 
                          biphasic = list(normal = F)),
                rel_cell = list(sigmoid = list(normal = T, low = F, high = F), 
                                biphasic = list(normal = F))
    )
  } else if(case %in% c("A", "B")) {
    fits = list(GR = list(sigmoid = list(normal = T, low = T, high = T, static = F, toxic = F), 
                          biphasic = list(normal = T)),
                rel_cell = list(sigmoid = list(normal = T, low = T, high = T), 
                          biphasic = list(normal = T))
                )
  }
  parameters = list()
  ## Define the biphasic dose-response function
  ## x is concentration, p is the parameter vector
  opfct_bi = function(x, p) {
    term1 = 1 + (p[1] + (1 - p[1])/(1 + (10^x / (10^p[2])) ^ p[3]))
    term2 = 1 + (p[4] + (1 - p[4])/(1 + (10^x / (10^p[5])) ^ p[6]))
    2^( 0.5*( log2(term1) + log2(term2) ) ) - 1
  }
  ## Define the sigmoidal (or logistic) dose-response function
  opfct_sig = function(x, p) {
    p[1] + (1 - p[1])/(1 + (10^x / (10^p[2])) ^ p[3])
  }
  ## Define the sigmoidal (or logistic) dose-response function for GR_toxic
  opfct_sig_GR_toxic = function(x, p) {
    p[1] + (0 - p[1])/(1 + (10^x / (10^p[2])) ^ p[3])
  }
  ##### start static vs. toxic curve fitting #######
  if(case == "static_vs_toxic") {
    grp = dplyr::syms(groupingVariables)
    data_grp = inputData %>% dplyr::group_by(experiment, !!!grp)
    data_grp_summ = data_grp %>% 
      dplyr::filter(concentration > 0) %>%
      dplyr::summarise(
        GR_static_mean = mean(GR_static, na.rm = TRUE),
        GR_toxic_mean = mean(GR_toxic, na.rm = TRUE),
        #rel_cell_mean = mean(rel_cell_count, na.rm = TRUE),
        #ctrl_cell_doublings = mean(ctrl_cell_doublings, na.rm = TRUE),
        #treated_cell_doublings = mean(treated_cell_doublings, na.rm = TRUE),
        concentration_points = length(unique(concentration)),
        conc = list(unique(concentration)),
        cc = list( log10(unique(concentration)) ),
        GR_static = list ( GR_static ),
        GR_toxic = list ( GR_toxic ),
        #rel_cell_count = list ( rel_cell_count ),
        concentration = list ( concentration )
      )
    ## get RSS of flat fit for GR and relative cell count curve
    data_grp_summ %<>% dplyr::mutate(
      RSS1_GR_static = sum( (unlist(GR_static) - GR_static_mean)^2, na.rm = TRUE ),
      RSS1_GR_toxic = sum( (unlist(GR_toxic) - GR_toxic_mean)^2, na.rm = TRUE )
      )
    ## calculate GRmax (or Emax) and GR_AOC (or AUC) directly from unfitted data points
    data_grp_conc = inputData %>% 
      dplyr::filter(concentration > 0) %>%
      dplyr::group_by(experiment, !!!grp, log10_concentration) %>%
      dplyr::summarise(
        GR_static_mean = mean(GR_static, na.rm = TRUE),
        GR_toxic_mean = mean(GR_toxic, na.rm = TRUE),
        cc = unique(log10_concentration)
      ) %>%
      ## get maximum (and second largest) concentration for each group
      dplyr::mutate(max_conc = list(sort(cc, decreasing = T)[1:2]) )
    data_grp_conc_max = data_grp_conc %>% 
      ## look at only two highest concentrations
      dplyr::filter( cc %in%  unlist(max_conc) ) %>%
      ## GRmax = min of avg. GR value at the two highest concentrations
      dplyr::summarise(
        GR_static_max = min(GR_static_mean, na.rm = T),
        GR_toxic_max = min(GR_toxic_mean, na.rm = T)
      ) %>% dplyr::ungroup()
    data_grp_conc_all = data_grp_conc %>% dplyr::arrange(cc) %>%
      dplyr::summarise(
        GR_static_AOC = caTools::trapz(x = cc, y = 1 - GR_static_mean),
        GR_toxic_AOC = caTools::trapz(x = cc, y = - GR_toxic_mean) ) %>%
      dplyr::ungroup()
    
    data_grp_conc = suppressMessages(
      dplyr::full_join(data_grp_conc_all, data_grp_conc_max, 
                       by = c("experiment", groupingVariables) )
      )
    data_grp_summ %<>% dplyr::group_by(experiment, !!!grp) %>%
      dplyr::mutate(
        p_sig_GR_static = list( tibble::tribble(
          ~parameter,                 ~lower,     ~prior,    ~upper,
          "GRinf",                       0,        0.5,         1, 
          "log10_GEC50", min(unlist(cc))-2, median(unlist(cc)), max(unlist(cc))+2,
          "h_GR",                         0.1,          2,         5
        ) ),
        p_sig_GR_toxic = list( tibble::tribble(
          ~parameter,                 ~lower,     ~prior,    ~upper,
          "GRinf",                       -1,        -0.5,         0, 
          "log10_GEC50", min(unlist(cc))-2, median(unlist(cc)), max(unlist(cc))+2,
          "h_GR",                         0.1,          2,         5
        ) )
      )
    sum_square_sig_GR_static = function(x, y, p) {sum((y - opfct_sig(x, p))^2)}
    sum_square_sig_GR_toxic = function(x, y, p) {sum((y - opfct_sig_GR_toxic(x, p))^2)}
    if(fits$GR$sigmoid$static) {
      ## Fit GR_static sigmoidal curve
      data_grp_summ$sig_fit_GR_static = lapply(1:dim(data_grp_summ)[1], function(i) {
        param_df = data_grp_summ$p_sig_GR_static[[i]]
        xx = data_grp_summ$concentration[[i]]
        yy = data_grp_summ$GR_static[[i]]
        startVec = param_df$prior
        psVec <- abs(startVec)
        psVec[psVec < 1e-4] <- 1
        fit = suppressMessages(try(optim(par = param_df$prior, control = list(maxit = 500, parscale = psVec),
                                         function(p, x, y) sum_square_sig_GR_static(x = log10(xx), y = yy, p = p),
                                         hessian = TRUE, method = "L-BFGS-B",
                                         lower = param_df$lower, upper = param_df$upper)))
        fit$parameters = param_df$parameter
        fit$lower = param_df$lower
        fit$upper = param_df$upper
        fit$prior = param_df$prior
        return(fit)
      })
    }

    ## Fit GR_toxic sigmoidal curve
    if(fits$GR$sigmoid$toxic) {
      data_grp_summ$sig_fit_GR_toxic = lapply(1:dim(data_grp_summ)[1], function(i) {
        param_df = data_grp_summ$p_sig_GR_toxic[[i]]
        xx = data_grp_summ$concentration[[i]]
        yy = data_grp_summ$GR_toxic[[i]]
        startVec = param_df$prior
        psVec <- abs(startVec)
        psVec[psVec < 1e-4] <- 1
        fit = suppressMessages(try(optim(par = param_df$prior, control = list(maxit = 500, parscale = psVec),
                        function(p, x, y) sum_square_sig_GR_toxic(x = log10(xx), y = yy, p = p),
                        hessian = TRUE, method = "L-BFGS-B",
                        lower = param_df$lower, upper = param_df$upper)))
        fit$parameters = param_df$parameter
        fit$lower = param_df$lower
        fit$upper = param_df$upper
        fit$prior = param_df$prior
        return(fit)
      })
    }
    fit_types = list(static = NULL, toxic = NULL)
    #parameters = list(GR = fit_types)
    
    for(x in c("sig_fit_GR_static", "sig_fit_GR_toxic")) {
      params = data_grp_summ[[x]][[1]]$parameters
      df = data_grp_summ %>% 
        dplyr::select(experiment, !!!grp, #ctrl_cell_doublings, 
          #treated_cell_doublings, 
          concentration_points) %>% 
        dplyr::ungroup() %>%
        dplyr::mutate_if(is.factor, as.character)
      ## get fitted curve parameters
      vals = sapply(data_grp_summ[[x]], function(y) { 
        if(class(y$par) == "numeric") { y$par } else { rep(NA, length(params)) }
      }) %>% t() %>% as.data.frame() %>%
        magrittr::set_colnames(params) %>%
        dplyr::mutate_if(is.factor, function(z) as.numeric(as.character(z)))
      df = cbind(df, vals)
      ## Get GRmax (Emax) and GR_AOC (AUC)
      if(grepl("GR_static", x)) { 
        df =  suppressMessages(
          dplyr::left_join(df, 
            data_grp_conc %>% dplyr::select(experiment, !!!grp, GR_static_max, GR_static_AOC),
              by = c("experiment", groupingVariables) ))
        df %<>% dplyr::rename(GRmax = GR_static_max, GR_AOC = GR_static_AOC) %>%
          dplyr::mutate(GR_metric = "GR_static")
      } else if(grepl("GR_toxic", x)) {
        df =  suppressMessages(
          dplyr::left_join(df, 
            data_grp_conc %>% dplyr::select(experiment, !!!grp, GR_toxic_max, GR_toxic_AOC),
              by = c("experiment", groupingVariables) ))
        df %<>% dplyr::rename(GRmax = GR_toxic_max, GR_AOC = GR_toxic_AOC) %>%
          dplyr::mutate(GR_metric = "GR_toxic")
      }
      ## Get RSS1 
      if(grepl("GR_static", x)) { df$RSS1 = data_grp_summ$RSS1_GR_static }
      if(grepl("GR_toxic", x)) { df$RSS1 = data_grp_summ$RSS1_GR_toxic }
      
      ## Get RSS2
      df$RSS2 = sapply(data_grp_summ[[x]], function(y) {
        if("value" %in% names(y) && is.numeric(y$value)) { y$value } else { NA }
      })
      
      Npara_flat = 1
      Npara = length(params)
      df$df1 = Npara_flat
      df$df2 = Npara
      if(grepl("_GR_static", x)) {
        df$n =  sapply(data_grp_summ$GR_static, function(y) return(length(na.omit(y))) )
      }
      if(grepl("_GR_toxic", x)) {
        df$n =  sapply(data_grp_summ$GR_toxic, function(y) return(length(na.omit(y))) )
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
      
      if(grepl("GR_static", x)) {
        df$flat = data_grp_summ$GR_static_mean
        df %<>%
          dplyr::mutate(GEC50 = 10^log10_GEC50) %>%
          dplyr::mutate(GR50 = GEC50*((1-GRinf)/(0.5-GRinf) - 1)^(1/h_GR)) %>%
          dplyr::mutate(log10_GR50 = log10(GR50))
        df %<>%
          dplyr::mutate(
            GR50 = case_when(
              fit == "flat" & flat >= 0.5 ~ Inf,
              fit == "flat" & flat < 0.5 ~ 0,
              fit != "flat" ~ GR50
            ),
            log10_GR50 = case_when(
              fit == "flat" & flat >= 0.5 ~ Inf,
              fit == "flat" & flat < 0.5 ~ -Inf,
              fit != "flat" ~ log10_GR50
            )
          )
      }
      if(grepl("GR_toxic", x)) {
        df$flat = data_grp_summ$GR_toxic_mean
        toxic_GR50_val = -0.5
        df %<>%
          dplyr::mutate(GEC50 = 10^log10_GEC50) %>%
          ## GR50 is tentatively set to when GR = -0.5 for toxic curves
          dplyr::mutate(GR50 = GEC50*((-GRinf)/(toxic_GR50_val-GRinf) - 1)^(1/h_GR)) %>%
          dplyr::mutate(log10_GR50 = log10(GR50))
        df %<>%
          dplyr::mutate(
            GR50 = case_when(
              fit == "flat" & flat >= toxic_GR50_val ~ Inf,
              fit == "flat" & flat < toxic_GR50_val ~ 0,
              fit != "flat" ~ GR50
              ),
            log10_GR50 = case_when(
              fit == "flat" & flat >= toxic_GR50_val ~ Inf,
              fit == "flat" & flat < toxic_GR50_val ~ -Inf,
              fit != "flat" ~ log10_GR50
            )
          )
      }
      df %<>%
        dplyr::mutate(
          GEC50 = case_when(
            fit == "flat" ~ 0,
            fit != "flat" ~ GEC50),
          log10_GEC50 = case_when(
            fit == "flat" ~ -Inf,
            fit != "flat" ~ log10_GEC50
          ),
          GRinf = case_when(
            fit == "flat" ~ GRmax,
            fit != "flat" ~ GRinf
          ),
          h_GR = case_when(
            fit == "flat" ~ 0.01,
            fit != "flat" ~ h_GR
          )
        )
      df_cols = colnames(df)
      cols_order = c("experiment", groupingVariables, "fit",
                    "GR_metric","GR50", "log10_GR50",
                    "GEC50", "log10_GEC50",
                    "GRmax", "GRinf", "h_GR", "GR_AOC")
      cols_order = cols_order[cols_order %in% df_cols]
      cols_left = setdiff(df_cols, cols_order)
      df = df[, c(cols_order, cols_left)]
    # Add GR50 = +/-Inf for any curves that don't reach GR = 0.5
    df %<>% dplyr::mutate(
      GR50 = case_when(
        is.na(GR50) & flat >= .5 ~ Inf,
        is.na(GR50) & flat < .5 ~ 0,
        !is.na(GR50) ~ GR50
      ),
      log10_GR50 = case_when(
        is.na(log10_GR50) & flat >= .5 ~ Inf,
        is.na(log10_GR50) & flat < .5 ~ -Inf,
        !is.na(log10_GR50) ~ log10_GR50
      )
    )
      
      if(x == "sig_fit_GR_static") { parameters$GR$sigmoid$static = df }
      if(x == "sig_fit_GR_toxic") { parameters$GR$sigmoid$toxic = df }
    }
  }
  ##### end static vs. toxic curve fitting #######
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
      p_sig_GR_low = list( tibble::tribble(
        ~parameter,                 ~lower,     ~prior,    ~upper,
        "GRinf",                       -1,        0.1,         1, 
        "log10_GEC50", min(unlist(cc))-2, median(unlist(cc)), max(unlist(cc))+2,
        "h_GR",                         0.1,          2,         5
      ) ),
      p_sig_GR_high = list( tibble::tribble(
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
  if(fits$GR$biphasic$normal) {
    ## Fit biphasic curve (GR)
    data_grp_summ$bi_fit_GR = lapply(1:dim(data_grp_summ)[1], function(i) {
      param_df = data_grp_summ$p_bi_GR[[i]]
      xx = data_grp_summ$concentration[[i]]
      yy = data_grp_summ$GRvalue[[i]]
      startVec = param_df$prior
      psVec <- abs(startVec)
      psVec[psVec < 1e-4] <- 1
      fit = suppressMessages(try(optim(par = param_df$prior, control = list(maxit = 500, parscale = psVec),
                      function(p, x, y) sum_square_bi(x = log10(xx), y = yy, p = p),
                      hessian = TRUE, method = "L-BFGS-B", 
                      lower = param_df$lower, upper = param_df$upper)))
      fit$parameters = param_df$parameter
      fit$lower = param_df$lower
      fit$upper = param_df$upper
      fit$prior = param_df$prior
      return(fit)
    })
  }
  ## Fit biphasic curve (relative cell count)
  if(fits$rel_cell$biphasic$normal) {
    data_grp_summ$bi_fit_rel = lapply(1:dim(data_grp_summ)[1], function(i) {
      param_df = data_grp_summ$p_bi_rel[[i]]
      xx = data_grp_summ$concentration[[i]]
      yy = data_grp_summ$rel_cell_count[[i]]
      startVec = param_df$prior
      psVec <- abs(startVec)
      psVec[psVec < 1e-4] <- 1
      fit = suppressMessages(try(optim(par = param_df$prior, control = list(maxit = 500, parscale = psVec),
                      function(p, x, y) sum_square_bi(x = log10(xx), y = yy, p = p),
                      hessian = TRUE, method = "L-BFGS-B", 
                      lower = param_df$lower, upper = param_df$upper)))
      fit$parameters = param_df$parameter
      fit$lower = param_df$lower
      fit$upper = param_df$upper
      fit$prior = param_df$prior
      return(fit)
    })
  }
  if(fits$GR$sigmoid$low) {
    ## Fit low sigmoidal fit (GR)
    data_grp_summ$sig_low_fit_GR = lapply(1:dim(data_grp_summ)[1], function(i) {
      param_df = data_grp_summ$p_bi_GR[[i]][1:3,]
      ## remove "_1" and "_2" from parameter names
      param_df %<>% dplyr::mutate(parameter = gsub("_1$", "", parameter))
      xx = data_grp_summ$concentration[[i]]
      yy = data_grp_summ$GRvalue[[i]]
      startVec = param_df$prior
      psVec <- abs(startVec)
      psVec[psVec < 1e-4] <- 1
      fit = suppressMessages(try(optim(par = param_df$prior, control = list(maxit = 500, parscale = psVec),                    function(p, x, y) sum_square_sig(x = log10(xx), y = yy, p = p),
                      hessian = TRUE, method = "L-BFGS-B", 
                      lower = param_df$lower, upper = param_df$upper)))
      fit$parameters = param_df$parameter
      fit$lower = param_df$lower
      fit$upper = param_df$upper
      fit$prior = param_df$prior
      return(fit)
    })
  }
  
  if(fits$GR$sigmoid$high) {
    ## Fit high sigmoidal fit (GR)
    data_grp_summ$sig_high_fit_GR = lapply(1:dim(data_grp_summ)[1], function(i) {
      param_df = data_grp_summ$p_bi_GR[[i]][4:6,]
      ## remove "_1" and "_2" from parameter names
      param_df %<>% dplyr::mutate(parameter = gsub("_2$", "", parameter))
      xx = data_grp_summ$concentration[[i]]
      yy = data_grp_summ$GRvalue[[i]]
      startVec = param_df$prior
      psVec <- abs(startVec)
      psVec[psVec < 1e-4] <- 1
      fit = suppressMessages(try(optim(par = param_df$prior, control = list(maxit = 500, parscale = psVec),
                      function(p, x, y) sum_square_sig(x = log10(xx), y = yy, p = p),
                      hessian = TRUE, method = "L-BFGS-B", 
                      lower = param_df$lower, upper = param_df$upper)))
      fit$parameters = param_df$parameter
      fit$lower = param_df$lower
      fit$upper = param_df$upper
      fit$prior = param_df$prior
      return(fit)
    })
  }
  if(fits$rel_cell$sigmoid$low) {
    ## Fit low sigmoidal fit (relative cell count)
    data_grp_summ$sig_low_fit_rel = lapply(1:dim(data_grp_summ)[1], function(i) {
      param_df = data_grp_summ$p_bi_rel[[i]][1:3,]
      ## remove "_1" and "_2" from parameter names
      param_df %<>% dplyr::mutate(parameter = gsub("_1$", "", parameter))
      xx = data_grp_summ$concentration[[i]]
      yy = data_grp_summ$rel_cell_count[[i]]
      startVec = param_df$prior
      psVec <- abs(startVec)
      psVec[psVec < 1e-4] <- 1
      fit = suppressMessages(try(optim(par = param_df$prior, control = list(maxit = 500, parscale = psVec),                    function(p, x, y) sum_square_sig(x = log10(xx), y = yy, p = p),
                      hessian = TRUE, method = "L-BFGS-B", 
                      lower = param_df$lower, upper = param_df$upper)))
      fit$parameters = param_df$parameter
      fit$lower = param_df$lower
      fit$upper = param_df$upper
      fit$prior = param_df$prior
      return(fit)
    })
  }
  if(fits$rel_cell$sigmoid$high) {
    ## Fit high sigmoidal fit (relative cell count)
    data_grp_summ$sig_high_fit_rel = lapply(1:dim(data_grp_summ)[1], function(i) {
      param_df = data_grp_summ$p_bi_rel[[i]][4:6,]
      ## remove "_1" and "_2" from parameter names
      param_df %<>% dplyr::mutate(parameter = gsub("_2$", "", parameter))
      xx = data_grp_summ$concentration[[i]]
      yy = data_grp_summ$rel_cell_count[[i]]
      startVec = param_df$prior
      psVec <- abs(startVec)
      psVec[psVec < 1e-4] <- 1
      fit = suppressMessages(try(optim(par = param_df$prior, control = list(maxit = 500, parscale = psVec),                    function(p, x, y) sum_square_sig(x = log10(xx), y = yy, p = p),
                      hessian = TRUE, method = "L-BFGS-B",
                      lower = param_df$lower, upper = param_df$upper)))
      fit$parameters = param_df$parameter
      fit$lower = param_df$lower
      fit$upper = param_df$upper
      fit$prior = param_df$prior
      return(fit)
    })
  }
  
  if(fits$GR$sigmoid$normal) {
    ## Fit normal sigmoidal fit (GR)
    data_grp_summ$sig_fit_GR = lapply(1:dim(data_grp_summ)[1], function(i) {
      param_df = data_grp_summ$p_sig_GR[[i]]
      xx = data_grp_summ$concentration[[i]]
      yy = data_grp_summ$GRvalue[[i]]
      startVec = param_df$prior
      psVec <- abs(startVec)
      psVec[psVec < 1e-4] <- 1
      fit = suppressMessages(try(optim(par = param_df$prior, control = list(maxit = 500, parscale = psVec),                    function(p, x, y) sum_square_sig(x = log10(xx), y = yy, p = p),
                      hessian = TRUE, method = "L-BFGS-B",
                      lower = param_df$lower, upper = param_df$upper)))
      fit$parameters = param_df$parameter
      fit$lower = param_df$lower
      fit$upper = param_df$upper
      fit$prior = param_df$prior
      return(fit)
    })
  }
  if(fits$rel_cell$sigmoid$normal) {
    ## Fit normal sigmoidal fit (relative cell count)
    data_grp_summ$sig_fit_rel = lapply(1:dim(data_grp_summ)[1], function(i) {
      param_df = data_grp_summ$p_sig_rel[[i]]
      xx = data_grp_summ$concentration[[i]]
      yy = data_grp_summ$rel_cell_count[[i]]
      startVec = param_df$prior
      psVec <- abs(startVec)
      psVec[psVec < 1e-4] <- 1
      fit = suppressMessages(try(optim(par = param_df$prior, control = list(maxit = 500, parscale = psVec),                    function(p, x, y) sum_square_sig(x = log10(xx), y = yy, p = p),
                      hessian = TRUE, method = "L-BFGS-B",
                      lower = param_df$lower, upper = param_df$upper)))
      fit$parameters = param_df$parameter
      fit$lower = param_df$lower
      fit$upper = param_df$upper
      fit$prior = param_df$prior
      return(fit)
    })
  }
  
  constraints_sig = list(normal = NULL, low = NULL, high = NULL)
  constraints_bi = list(normal = NULL)
  fit_types = list(sigmoid = constraints_sig, biphasic = constraints_bi)
  #parameters = list(GR = fit_types, rel_cell = fit_types)
  
  for(x in c("sig_fit_rel", "sig_fit_GR", "sig_low_fit_rel", "sig_low_fit_GR",
             "sig_high_fit_rel", "sig_high_fit_GR", "bi_fit_rel", "bi_fit_GR")) {
    if(!x %in% colnames(data_grp_summ)) next
    params = data_grp_summ[[x]][[1]]$parameters
    df = data_grp_summ %>% dplyr::select(experiment, !!!grp, ctrl_cell_doublings, treated_cell_doublings,
                                         concentration_points) %>% dplyr::ungroup()
    ## get fitted curve parameters
    vals = sapply(data_grp_summ[[x]], function(y) { 
      if(class(y$par) == "numeric") { y$par } else { rep(NA, length(params)) }
    }) %>% t() %>% as.data.frame() %>%
      magrittr::set_colnames(params)
    df = cbind(df, vals)
    ## Get GRmax (Emax) and GR_AOC (AUC)
    if(grepl("_GR", x)) { 
      df =  suppressMessages(dplyr::full_join(df, data_grp_conc %>% dplyr::select(experiment, !!!grp, GRmax, GR_AOC),
                                              by = c("experiment", groupingVariables) ))
      if(grepl("bi_fit", x)) {
        # opfct_bi = function(x, p) {
        #   term1 = 1 + (p[1] + (1 - p[1])/(1 + (10^x / (10^p[2])) ^ p[3]))
        #   term2 = 1 + (p[4] + (1 - p[4])/(1 + (10^x / (10^p[5])) ^ p[6]))
        #   2^( 0.5*( log2(term1) + log2(term2) ) ) - 1
        # }
        #column with lists of parametersdata_grp_summ$bi_fit_GR
        bi_50 = function(xx,pp) { opfct_bi(x = xx, p = pp) - 0.5 }
        GR50_temp = sapply(data_grp_summ$bi_fit_GR, function(yy) {
          if(!is.numeric(yy$par)) return(NA)
          ret = try( uniroot(bi_50, pp = yy$par, interval = c(-8,6), tol = 10e-8)$root )
          if(class(ret) == "numeric") {
            return(ret)
          } else {
            return(NA)
          }
        })
        # check_GR50 = sapply(1:length(data_grp_summ$bi_fit_GR), function(ii) {
        #   ret = opfct_bi(x = GR50[ii], p = data_grp_summ$bi_fit_GR[[ii]]$par)
        #   if(class(ret) == "numeric") {
        #     return(ret)
        #   } else {
        #     return(NA)
        #   }
        # })
        df %<>%
          dplyr::mutate(GEC50_1 = 10^log10_GEC50_1,
                        GEC50_2 = 10^log10_GEC50_2,
                        GR50 = 10^GR50_temp,
                        log10_GR50 = GR50_temp
                        )
      } else {
        df %<>%
          dplyr::mutate(GEC50 = 10^log10_GEC50) %>%
          dplyr::mutate(GR50 = GEC50*((1-GRinf)/(0.5-GRinf) - 1)^(1/h_GR)) %>%
          dplyr::mutate(log10_GR50 = log10(GR50))
      }
    }
    if(grepl("_rel", x)) {
      df = suppressMessages(dplyr::full_join(df, data_grp_conc %>% dplyr::select(experiment, !!!grp, Emax, AUC),
                                             by = c("experiment", groupingVariables) ) )
      if(grepl("bi_fit", x)) {
        bi_50 = function(xx,pp) { opfct_bi(x = xx, p = pp) - 0.5 }
        IC50_temp = sapply(data_grp_summ$bi_fit_rel, function(yy) {
          ret = try( uniroot(bi_50, pp = yy$par, interval = c(-8,6), tol = 10e-8)$root )
          if(class(ret) == "numeric") {
            return(ret)
          } else {
            return(NA)
          }
        })
        
        df %<>% 
          dplyr::mutate(EC50_1 = 10^log10_EC50_1,
                        EC50_2 = 10^log10_EC50_2,
                        IC50 = 10^IC50_temp,
                        log10_IC50 = IC50_temp
                        )
      } else {
        df %<>% 
          dplyr::mutate(EC50 = 10^log10_EC50) %>%
          dplyr::mutate(IC50 = EC50*((1-Einf)/(0.5-Einf) - 1)^(1/h) ) %>%
          dplyr::mutate(log10_IC50 = log10(IC50) )
      }
    }
    ## Get RSS1 
    if(grepl("_GR", x)) { df$RSS1 = data_grp_summ$RSS1_GR }
    if(grepl("_rel", x)) { df$RSS1 = data_grp_summ$RSS1_rel_cell }
    ## Get RSS2
    df$RSS2 = sapply(data_grp_summ[[x]], function(y) { 
      if(class(y$value) == "numeric") { y$value } else { NA }
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
    
    if(grepl("_GR", x)) { df$flat = data_grp_summ$GR_mean }
    if(grepl("_rel", x)) { df$flat = data_grp_summ$rel_cell_mean }
    
    ## adjust parameters for flat fits
    if(grepl("_GR", x)) { ### GR fits
      if(grepl("bi_fit", x)) { ### biphasic fits
        df %<>%
          dplyr::mutate(
            GEC50_1 = case_when(
              fit == "flat" ~ 0,
              fit != "flat" ~ GEC50_1),
            GEC50_2 = case_when(
              fit == "flat" ~ 0,
              fit != "flat" ~ GEC50_2),
            log10_GEC50_1 = case_when(
              fit == "flat" ~ -Inf,
              fit != "flat" ~ log10_GEC50_1
            ),
            log10_GEC50_2 = case_when(
              fit == "flat" ~ -Inf,
              fit != "flat" ~ log10_GEC50_2
            ),
            GRinf_1 = case_when(
              fit == "flat" ~ GRmax,
              fit != "flat" ~ GRinf_1
            ),
            GRinf_2 = case_when(
              fit == "flat" ~ GRmax,
              fit != "flat" ~ GRinf_2
            ),
            GR50 = case_when(
              fit == "flat" & flat >= .5 ~ Inf,
              fit == "flat" & flat < .5 ~ 0,
              fit != "flat" ~ GR50
            ),
            log10_GR50 = case_when(
              fit == "flat" & flat >= .5 ~ Inf,
              fit == "flat" & flat < .5 ~ -Inf,
              fit != "flat" ~ log10_GR50
            ),
            h_GR_1 = case_when(
              fit == "flat" ~ 0.01,
              fit != "flat" ~ h_GR_1
            ),
            h_GR_2 = case_when(
              fit == "flat" ~ 0.01,
              fit != "flat" ~ h_GR_2
            )
          )
        df_cols = colnames(df)
        cols_order = c("experiment", groupingVariables, "fit", "GR50", "log10_GR50",
                       "GEC50_1", "log10_GEC50_1", "GEC50_2", "log10_GEC50_2",
                       "GRmax", "GRinf_1", "GRinf_2", "h_GR_1", "h_GR_2", "GR_AOC")
        cols_left = setdiff(df_cols, cols_order)
        df = df[, c(cols_order, cols_left)]
      } else { ### sigmoid fits
        df %<>%
          dplyr::mutate(
            GEC50 = case_when(
              fit == "flat" ~ 0,
              fit != "flat" ~ GEC50),
            log10_GEC50 = case_when(
              fit == "flat" ~ -Inf,
              fit != "flat" ~ log10_GEC50
            ),
            GRinf = case_when(
              fit == "flat" ~ GRmax,
              fit != "flat" ~ GRinf
            ),
            GR50 = case_when(
              fit == "flat" & flat >= .5 ~ Inf,
              fit == "flat" & flat < .5 ~ 0,
              fit != "flat" ~ GR50
            ),
            log10_GR50 = case_when(
              fit == "flat" & flat >= .5 ~ Inf,
              fit == "flat" & flat < .5 ~ -Inf,
              fit != "flat" ~ log10_GR50
            ),
            h_GR = case_when(
              fit == "flat" ~ 0.01,
              fit != "flat" ~ h_GR
            )
          )
        df_cols = colnames(df)
        cols_order = c("experiment", groupingVariables, "fit", "GR50", "log10_GR50",
                       "GEC50", "log10_GEC50",
                       "GRmax", "GRinf", "h_GR", "GR_AOC")
        cols_left = setdiff(df_cols, cols_order)
        df = df[, c(cols_order, cols_left)]
      }
      # Add GR50 = +/-Inf for any curves that don't reach GR = 0.5
      df %<>% dplyr::mutate(
        GR50 = case_when(
          is.na(GR50) & flat >= .5 ~ Inf,
          is.na(GR50) & flat < .5 ~ 0,
          !is.na(GR50) ~ GR50
          ),
        log10_GR50 = case_when(
          is.na(log10_GR50) & flat >= .5 ~ Inf,
          is.na(log10_GR50) & flat < .5 ~ -Inf,
          !is.na(log10_GR50) ~ log10_GR50
          )
        )
    } else { ### relative cell count fits
      if(grepl("bi_fit", x)) { ### biphasic fits
        df %<>%
          dplyr::mutate(
            EC50_1 = case_when(
              fit == "flat" ~ 0,
              fit != "flat" ~ EC50_1),
            EC50_2 = case_when(
              fit == "flat" ~ 0,
              fit != "flat" ~ EC50_2),
            log10_EC50_1 = case_when(
              fit == "flat" ~ -Inf,
              fit != "flat" ~ log10_EC50_1
            ),
            log10_EC50_2 = case_when(
              fit == "flat" ~ -Inf,
              fit != "flat" ~ log10_EC50_2
            ),
            Einf_1 = case_when(
              fit == "flat" ~ Emax,
              fit != "flat" ~ Einf_1
            ),
            Einf_2 = case_when(
              fit == "flat" ~ Emax,
              fit != "flat" ~ Einf_2
            ),
            IC50 = case_when(
              fit == "flat" & flat >= .5 ~ Inf,
              fit == "flat" & flat < .5 ~ 0,
              fit != "flat" ~ IC50
            ),
            log10_IC50 = case_when(
              fit == "flat" & flat >= .5 ~ Inf,
              fit == "flat" & flat < .5 ~ -Inf,
              fit != "flat" ~ log10_IC50
            ),
            h_1 = case_when(
              fit == "flat" ~ 0.01,
              fit != "flat" ~ h_1
            ),
            h_2 = case_when(
              fit == "flat" ~ 0.01,
              fit != "flat" ~ h_2
            )
          )
        df_cols = colnames(df)
        cols_order = c("experiment", groupingVariables, "fit", "IC50", "log10_IC50",
                       "EC50_1", "log10_EC50_1", "EC50_2", "log10_EC50_2",
                       "Emax", "Einf_1", "Einf_2", "h_1", "h_2", "AUC")
        cols_left = setdiff(df_cols, cols_order)
        df = df[, c(cols_order, cols_left)]
      } else { ### sigmoid fits
        df %<>%
          dplyr::mutate(
            EC50 = case_when(
              fit == "flat" ~ 0,
              fit != "flat" ~ EC50),
            log10_EC50 = case_when(
              fit == "flat" ~ -Inf,
              fit != "flat" ~ log10_EC50
            ),
            Einf = case_when(
              fit == "flat" ~ Emax,
              fit != "flat" ~ Einf
            ),
            IC50 = case_when(
              fit == "flat" & flat >= .5 ~ Inf,
              fit == "flat" & flat < .5 ~ 0,
              fit != "flat" ~ IC50
            ),
            log10_IC50 = case_when(
              fit == "flat" & flat >= .5 ~ Inf,
              fit == "flat" & flat < .5 ~ -Inf,
              fit != "flat" ~ log10_IC50
            ),
            h = case_when(
              fit == "flat" ~ 0.01,
              fit != "flat" ~ h
            )
          )
        df_cols = colnames(df)
        cols_order = c("experiment", groupingVariables, "fit", "IC50", "log10_IC50",
                       "EC50", "log10_EC50",
                       "Emax", "Einf", "h", "AUC")
        cols_left = setdiff(df_cols, cols_order)
        df = df[, c(cols_order, cols_left)]
      }
      # Add IC50 = +/-Inf for any curves that don't reach 0.5
      df %<>% dplyr::mutate(
        IC50 = case_when(
          is.na(IC50) & flat >= .5 ~ Inf,
          is.na(IC50) & flat < .5 ~ 0,
          !is.na(IC50) ~ IC50
        ),
        log10_IC50 = case_when(
          is.na(log10_IC50) & flat >= .5 ~ Inf,
          is.na(log10_IC50) & flat < .5 ~ -Inf,
          !is.na(log10_IC50) ~ log10_IC50
        )
      )
    }
    
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