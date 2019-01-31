.check = function(inputData, case) {
  message = NULL # an error message, if applicable
  initial_count = TRUE # a logical for whether initial cell count is provided
  input_cols = colnames(inputData)
  caseA = c('concentration', 'cell_count', 'cell_count__ctrl',
            'cell_count__time0')
  caseA_div_time = c('concentration', 'cell_count','cell_count__ctrl',
                     'treatment_duration','division_time')
  if(case == "static_vs_toxic") {
    ### check input
    counts = c("dead_count__time0", "cell_count__time0", "cell_count", 
               "dead_count", "dead_count__ctrl", "cell_count__ctrl")
    counts = c(counts, "concentration", "time")
    if(sum(!counts %in% colnames(inputData)) != 0) {
      missing = counts[!counts %in% colnames(inputData)]
      missing = paste0(missing, collapse = ", ")
      stop(paste0("Missing columns in inputData: ", missing))
    }
  }
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
  if(case == "B") {
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