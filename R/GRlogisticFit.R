.GRlogisticFit = function(inputData, groupingVariables, force = FALSE,
                          cap = FALSE) {
  # declaring values NULL to avoid note on package check
  #experiment = NULL
  # GR curve parameters
  GEC50 = NULL
  GRinf = NULL
  h_GR = NULL
  # Relative cell count curve parameters
  EC50 = NULL
  Einf = NULL
  h = NULL
  experiments = levels(inputData$experiment)
  
  parameters = matrix(data = NA, ncol = 3, nrow = length(experiments))
  parameters = as.data.frame(parameters)
  
  parameters2 = matrix(data = NA, ncol = 3, nrow = length(experiments))
  parameters2 = as.data.frame(parameters)
  
  colnames(parameters) = c('h_GR','GRinf','GEC50')
  colnames(parameters2) = c('h', 'Einf', 'EC50')
  if(length(groupingVariables) > 0) {
    metadata = matrix(data = NA, ncol = length(groupingVariables),
                      nrow = length(experiments))
    metadata = as.data.frame(metadata)
    colnames(metadata) = groupingVariables
  } else {
    metadata = NULL
  }
  # More GR curve parameters
  pval_GR = NULL
  GRmax = NULL
  GR_mean = NULL
  AOC = NULL
  R_square_GR = NULL
  
  # More Relative cell count curve parameters
  pval_rel_cell = NULL
  Emax = NULL
  rel_cell_mean = NULL
  AUC = NULL
  R_square_rel_cell = NULL
  
  # other curve parameters
  concentration_points = NULL
  ctrl_cell_doublings = NULL
  for(i in 1:length(experiments)) {
    # print(i)
    data_exp = inputData[inputData$experiment == experiments[i], ]
    
    if(!is.null(metadata)) {
      metadata[i,] = data_exp[1,groupingVariables, drop = FALSE]
    }
    GR_mean[i] = mean(data_exp$GRvalue, na.rm = TRUE)
    rel_cell_mean[i] = mean(data_exp$rel_cell_count, na.rm = TRUE)
    # calculate avg number of cell doublings in control and treated experiments
    ctrl_cell_doublings[i] = mean(data_exp$ctrl_cell_doublings, na.rm = TRUE)
    #===== constrained fit GR curve ============
    c = unique(data_exp$concentration)
    priors = c(2, 0.1, stats::median(c))
    lower = c(.1, -1, min(c)*1e-2)
    upper = c(5, 1, max(c)*1e2)
    #===== constrained fit Relative cell count curve ============
    priors_rel_cell = c(2, 0.1, stats::median(c))
    lower_rel_cell = c(.1, 0, min(c)*1e-2)
    upper_rel_cell = c(5, 1, max(c)*1e2)
    if(dim(data_exp)[1] > 1) {
      controls = drc::drmc()
      controls$relTol = 1e-06
      controls$errorm = FALSE
      controls$noMessage = TRUE
      controls$rmNA = TRUE
      # GR curve fitting
      output_model_new = try(drc::drm(
        GRvalue~log10_concentration, experiment, data=data_exp, logDose = 10,
        fct=drc::LL.3u(names=c('h_GR','GRinf','GEC50')), start = priors,
        lowerl = lower, upperl = upper, control = controls,
        na.action = na.omit))
      if(class(output_model_new) == "drc" && 
         !is.null(stats::coef(output_model_new)) && 
         !is.null(stats::residuals(output_model_new))
         ) {
        parameters[i,] = c(as.numeric(stats::coef(output_model_new)))
        # F-test for the significance of the sigmoidal fit
        Npara = 3 # N of parameters in the growth curve
        Npara_flat = 1 # F-test for the models
        RSS2 = sum(stats::residuals(output_model_new)^2, na.rm = TRUE)
        RSS1 = sum((data_exp$GRvalue - mean(data_exp$GRvalue, na.rm = TRUE))^2,
                   na.rm = TRUE)
        df1 = (Npara - Npara_flat)
        df2 = (length(na.omit(data_exp$GRvalue)) - Npara + 1)
        f_value = ((RSS1-RSS2)/df1)/(RSS2/df2)
        f_pval = stats::pf(f_value, df1, df2, lower.tail = FALSE)
        pval_GR[i] = f_pval
        R_square_GR[i] = 1 - RSS2/RSS1
      } else {
        parameters[i,] = NA
        pval_GR[i] = NA
        R_square_GR[i] = NA
      }
      # Relative cell count curve fitting
      output_model_new_rel_cell = try(drc::drm(
        rel_cell_count~log10_concentration, experiment, data=data_exp, logDose = 10,
        fct=drc::LL.3u(names=c('h','Einf','EC50')), start = priors_rel_cell,
        lowerl = lower_rel_cell, upperl = upper_rel_cell, control = controls,
        na.action = na.omit))
      if(class(output_model_new_rel_cell) == "drc" && 
         !is.null(stats::coef(output_model_new_rel_cell)) && 
         !is.null(stats::residuals(output_model_new_rel_cell))
         ) {
        parameters2[i,] = c(as.numeric(stats::coef(output_model_new_rel_cell)))
        # F-test for the significance of the sigmoidal fit
        Npara = 3 # N of parameters in the growth curve
        Npara_flat = 1 # F-test for the models
        RSS2 = sum(stats::residuals(output_model_new_rel_cell)^2, na.rm = TRUE)
        RSS1 = sum((data_exp$rel_cell_count - mean(data_exp$rel_cell_count,
                   na.rm = TRUE))^2, na.rm = TRUE)
        df1 = (Npara - Npara_flat)
        df2 = (length(na.omit(data_exp$rel_cell_count)) - Npara + 1)
        f_value = ((RSS1-RSS2)/df1)/(RSS2/df2)
        f_pval = stats::pf(f_value, df1, df2, lower.tail = FALSE)
        pval_rel_cell[i] = f_pval
        R_square_rel_cell[i] = 1 - RSS2/RSS1
      } else {
        parameters2[i,] = NA
        pval_rel_cell[i] = NA
        R_square_rel_cell[i] = NA
      }
    }
    #==================================

    # Trapezoid rule for integration of GR_AOC
    GRavg = NULL
    rel_cell_avg = NULL
    
    concs = sort(unique(data_exp$concentration))
    l = length(concs)
    concentration_points[i] = l
    
    for(j in 1:l) {
      data_trapz = data_exp[data_exp$concentration == concs[j],]
      GRavg[j] = mean(data_trapz$GRvalue, na.rm = TRUE)
      rel_cell_avg[j] = mean(data_trapz$rel_cell_count, na.rm = TRUE)
    }
    
    GRmax[i] = min(GRavg[c(l,l-1)], na.rm = TRUE)
    Emax[i] = min(rel_cell_avg[c(l,l-1)], na.rm = TRUE)
    
    diff_vector = diff(log10(concs), lag = 1)
    conc_range = log10(concs[length(concs)]) - log10(concs[1])
    AOC[i] = sum((1 - (GRavg[1:(length(GRavg)-1)]+GRavg[2:length(GRavg)])/2)*
                   diff_vector, na.rm = TRUE)/conc_range
    AUC[i] = sum(((rel_cell_avg[1:(length(rel_cell_avg)-1)]+rel_cell_avg[2:length(rel_cell_avg)])/2)*
                   diff_vector, na.rm = TRUE)/conc_range
  }

  parameters = cbind(parameters, parameters2)
  parameters$ctrl_cell_doublings = ctrl_cell_doublings
  parameters$concentration_points = concentration_points
  # Calculate GR50 from parameters
  parameters$GR50 = with(parameters,GEC50*((1-GRinf)/(0.5-GRinf) - 1)^(1/h_GR))
  parameters$GRmax = GRmax
  parameters$GR_AOC = AOC
  parameters$r2_GR = R_square_GR
  # Calculate IC50 from parameters
  parameters$IC50 = with(parameters,EC50*((1-Einf)/(0.5-Einf) - 1)^(1/h))
  parameters$Emax = Emax
  parameters$AUC = AUC
  parameters$r2_rel_cell = R_square_rel_cell
  
  if(is.null(pval_GR)) {pval_GR = NA}
  if(is.null(pval_rel_cell)) {pval_rel_cell = NA}
  parameters$pval_GR = pval_GR
  parameters$pval_rel_cell = pval_rel_cell
  # Re-order rows to match reference_output
  parameters$experiment = experiments
  # Threshold for F-test pval_GR
  pcutoff = ifelse(force == FALSE, .05 , 1)
  for(i in 1:dim(parameters)[1]) {
    # Flat or sigmoid fit for GR curve
    if(!is.na(parameters$pval_GR[i])) {
      parameters$fit_GR[i] = ifelse(parameters$pval_GR[i] >= pcutoff |
                                 is.na(parameters$GEC50[i]), "flat","sigmoid")
    } else {
      parameters$fit_GR[i] = ifelse(is.na(parameters$GEC50[i]), "flat",
                                    "sigmoid")
    }
    # Flat or sigmoid fit for relative cell count curve
    if(!is.na(parameters$pval_rel_cell[i])) {
      parameters$fit_rel_cell[i] = ifelse(parameters$pval_rel_cell[i] >= pcutoff |
                                      is.na(parameters$EC50[i]), "flat",
                                    "sigmoid")
    } else {
      parameters$fit_rel_cell[i] = ifelse(is.na(parameters$EC50[i]), "flat",
                                    "sigmoid")
    }
  }
  # changed to above code to deal with NAs
  # Add values for flat fits: GEC50 = 0, h_GR = 0.01 and GR50 = +/- Inf
  
  parameters$flat_fit_GR = GR_mean
  parameters$flat_fit_rel_cell = rel_cell_mean
  for(i in 1:dim(parameters)[1]) {
    if(parameters$fit_GR[i] == "flat") {
      parameters$GEC50[i] = 0
      parameters$h_GR[i] = 0.01
      parameters$GR50[i] = ifelse(parameters$flat_fit_GR[i] > .5, Inf, -Inf)
      parameters$GRinf[i] = parameters$GRmax[i]
    }
  }
  for(i in 1:dim(parameters)[1]) {
    if(parameters$fit_rel_cell[i] == "flat") {
      parameters$EC50[i] = 0
      parameters$h[i] = 0.01
      parameters$IC50[i] = ifelse(parameters$flat_fit_rel_cell[i] > .5, Inf, -Inf)
      parameters$Einf[i] = parameters$Emax[i]
    }
  }
  # Add GR50 = +/-Inf for any curves that don't reach GR = 0.5
  for(i in 1:dim(parameters)[1]) {
    if(is.na(parameters$GR50[i])) {
      parameters$GR50[i] = ifelse(parameters$flat_fit_GR[i] > .5, Inf, -Inf)
    }
    if(is.na(parameters$IC50[i])) {
      parameters$IC50[i] = ifelse(parameters$flat_fit_rel_cell[i] > .5, Inf, -Inf)
    }
  }
  for(i in 1:dim(parameters)[1]) {
    if(parameters$fit_GR[i] == "sigmoid") {
      parameters$flat_fit_GR[i] = NA
    }
    if(parameters$fit_rel_cell[i] == "sigmoid") {
      parameters$flat_fit_rel_cell[i] = NA
    }
  }
  parameters = parameters[,c('GR50','GRmax','GR_AOC','GEC50','GRinf','h_GR',
                             'r2_GR','pval_GR', 'fit_GR',
                             'flat_fit_GR', 'IC50', 'Emax', 'AUC', 'EC50',
                             'Einf', 'h', 'r2_rel_cell', 'pval_rel_cell', 'fit_rel_cell',
                             'flat_fit_rel_cell','experiment',
                             'concentration_points', 'ctrl_cell_doublings')]
  if(!is.null(metadata)) {
    parameters = cbind(metadata, parameters)
  }
  return(parameters)
}