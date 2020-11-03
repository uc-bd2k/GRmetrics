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
  return(input_edited)
}