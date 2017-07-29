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
  input_edited$GR = GR
  input_edited$rel_cell_count = rel_cell_count
  input_edited$ctrl_cell_doublings = log2nn_ctrl
  tmp<-input_edited[,groupingVariables, drop = FALSE]
  experimentNew = (apply(tmp,1, function(x) (paste(x,collapse=" "))))
  if(cap == TRUE) {
    input_edited$GR[input_edited$GR > 1] = 1
    input_edited$rel_cell_count[input_edited$rel_cell_count > 1] = 1
  }
  if(length(groupingVariables) > 0) {
    input_edited$experiment = as.factor(experimentNew)
  } else {
    input_edited$experiment = as.factor("All Data")
  }
  return(input_edited)
}

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
    GR_mean[i] = mean(data_exp$GR, na.rm = TRUE)
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
        GR~log10_concentration, experiment, data=data_exp, logDose = 10,
        fct=drc::LL.3u(names=c('h_GR','GRinf','GEC50')), start = priors,
        lowerl = lower, upperl = upper, control = controls,
        na.action = na.omit))
      if(class(output_model_new)!="try-error") {
        parameters[i,] = c(as.numeric(stats::coef(output_model_new)))
        # F-test for the significance of the sigmoidal fit
        Npara = 3 # N of parameters in the growth curve
        Npara_flat = 1 # F-test for the models
        RSS2 = sum(stats::residuals(output_model_new)^2, na.rm = TRUE)
        RSS1 = sum((data_exp$GR - mean(data_exp$GR, na.rm = TRUE))^2,
                   na.rm = TRUE)
        df1 = (Npara - Npara_flat)
        df2 = (length(na.omit(data_exp$GR)) - Npara + 1)
        f_value = ((RSS1-RSS2)/df1)/(RSS2/df2)
        f_pval = stats::pf(f_value, df1, df2, lower.tail = FALSE)
        pval_GR[i] = f_pval
        R_square_GR[i] = 1 - RSS2/RSS1
      }
      # Relative cell count curve fitting
      output_model_new_rel_cell = try(drc::drm(
        rel_cell_count~log10_concentration, experiment, data=data_exp, logDose = 10,
        fct=drc::LL.3u(names=c('h','Einf','EC50')), start = priors_rel_cell,
        lowerl = lower_rel_cell, upperl = upper_rel_cell, control = controls,
        na.action = na.omit))
      if(class(output_model_new_rel_cell)!="try-error") {
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
      GRavg[j] = mean(data_trapz$GR, na.rm = TRUE)
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
    if(initial_count) {
      time0 = inputData[inputData$time == 0, c(keys, 'cell_count')]
    }
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
  parameter_table = .GRlogisticFit(gr_table, groupingVariables, force, cap)

  colData = parameter_table[ ,c(groupingVariables, 'fit_GR', 'fit_rel_cell',
                                'experiment', 'concentration_points')]
  rownames(colData) = colData$experiment
  colData = S4Vectors::DataFrame(colData)
  
  Metric = c('ctrl_cell_doublings','GR50','GRmax','GR_AOC','GEC50','GRinf',
             'h_GR','r2_GR','pval_GR','flat_fit_GR', 
              'IC50', 'Emax', 'AUC', 'EC50','Einf', 'h', 
              'r2_rel_cell', 'pval_rel_cell', 'flat_fit_rel_cell')
  assays = parameter_table[ , Metric]
  rownames(assays) = parameter_table$experiment
  assays = t(assays)

  Description = c(
    "The number of cell doublings in the control population during the assay",
    "The concentration at which GR(c) = 0.5",
    "The maximal effect of the drug (minimal GR value)",
    "The 'Area Over the Curve' - The area between the line GR = 1 and the curve, similar to traditional AUC",
    "The concentration at half-maximal effect (growth rate normalized)",
    "The asymptotic effect of the drug (growth rate normalized)",
    "The Hill coefficient of the fitted (GR) curve, which reflects how steep the (GR) dose response curve is",
    "The coefficient of determination - essentially how well the (GR) curve fits to the data points",
    "The p-value of the F-test comparing the fit of the (GR) curve to a horizontal line fit",
    "For data that doesn't significantly fit better to a curve than a horizontal line fit, the y value (GR) of the flat line", 
    "The concentration at which relative cell count = 0.5",
    "The maximal effect of the drug (minimal relative cell count value)",
    "The 'Area Under the Curve' - The area below the fitted (traditional) dose response curve",
    "The concentration at half-maximal effect (not growth rate normalized)",
    "The asymptotic effect of the drug (not growth rate normalized)",
    "The Hill coefficient of the fitted (traditional) dose response curve, which reflects how steep the (traditional) dose response curve is",
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
            rowData = rowData, metadata = list(gr_table, groupingVariables))
  return(output)
}

