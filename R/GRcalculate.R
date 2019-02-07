.GRcalculate = function(inputData, groupingVariables, cap = FALSE, case = "A",
                        initial_count){
  # declaring values NULL to avoid note on package check
  if(case %in% c("A", "B")) {
    cell_count = NULL
    cell_count__time0 = NULL
    cell_count__ctrl = NULL
    duration = NULL
    division_time = NULL
    if(initial_count) {
      log2nn = with(inputData, log2(cell_count/cell_count__time0))
      log2nn_ctrl = with(inputData, log2(cell_count__ctrl/cell_count__time0))
      GR = 2^(log2nn/log2nn_ctrl) - 1
    } 
    if(!initial_count) {
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
    if(initial_count) {
      input_edited$treated_cell_doublings = log2nn
    }
    if(!initial_count) {
      #input_edited$treated_cell_doublings = log2_rel + log2nn_ctrl
      input_edited$treated_cell_doublings = NA
    }
    tmp<-input_edited[,groupingVariables, drop = FALSE]
    experimentNew = (apply(tmp,1, function(x) (paste(x,collapse=" "))))
    if(cap == TRUE) {
      input_edited$GRvalue[input_edited$GRvalue > 1] = 1
      input_edited$rel_cell_count[input_edited$rel_cell_count > 1] = 1
    }
    if(length(groupingVariables) > 0) {
      input_edited$experiment = experimentNew
    } else {
      input_edited$experiment = "All Data"
    }
    return(as.data.frame(input_edited))
  }
  if(case == "static_vs_toxic") {
    ### see lines 152-213 in https://github.com/datarail/DrugResponse/blob/1d0fd1898e797b74f611eb61dbd6926ac715bb83/MATLAB/import_columbus/Process_CellCountData2.m
    ### Note: what about NA values in dead_count and cell_count?? these will propagate currently, but we may want to set them to 0 instead?
    inputData %<>% 
      dplyr::mutate(too_few = dead_count + cell_count < .95*(cell_count__time0 + dead_count__time0)) %>%
      dplyr::mutate(too_many = concentration > 0 &
        dead_count + cell_count > 1.15*(cell_count__ctrl + dead_count__ctrl))
    ## if missing cells, assign to dead cells
    inputData %<>% 
      dplyr::mutate( dead_count = ifelse(!too_few, dead_count, 
        pmax(dead_count, cell_count__time0 + dead_count__time0 - floor(.95*cell_count) + 1, na.rm = T)))
    # if too many dead cells -> subtract dead cells
    inputData %<>% 
      dplyr::mutate( dead_count = ifelse(!too_many, dead_count, 
        pmin(dead_count, cell_count__ctrl + dead_count__ctrl - ceiling(1.15*cell_count) - 1, na.rm = T)))
    
    num_cond_true = sum(inputData$too_few == T, na.rm = T)
    num_cond_missing = sum(is.na(inputData$too_few))
    if(num_cond_true != 0) {
      msg = paste0(num_cond_true, " conditions with too few cells -> assigned as dead cells")
      warning(msg)
    }
    if(num_cond_missing != 0) {
      msg = paste0(num_cond_missing, " conditions with missing cell counts")
      warning(msg)
    }
    
    ### how should we treat NAs here?? na.rm or not?
    inputData %<>% dplyr::mutate(
      # Dratio = (increase in dead cells)/(increase in alive cells)
      # Dratio_ctrl = same for ctrl cells
      # Note: minimum numerator allowed is 1, meaning the case of more dead cells at time0
      # is not allowed. Therefore numerator is always positive, so the denominator, that is
      # the increase or decrease in alive cells, controls the sign.
      Dratio = pmax(dead_count - dead_count__time0, 1, na.rm = T)/(cell_count - cell_count__time0),
      Dratio_ctrl = pmax(dead_count__ctrl - dead_count__time0, 1, na.rm = T)/(cell_count__ctrl - cell_count__time0),
      gr = log2(cell_count/cell_count__time0),
      gr_ctrl = log2(cell_count__ctrl/cell_count__time0)
    )
    ## These are slightly different from the formulas in the slide-deck, which are correct?
    inputData %<>% dplyr::mutate(
      GR_s = 2^( (1 + Dratio)*gr/( (1 + Dratio_ctrl)*gr_ctrl) ) - 1,
      GR_d = 2^( ( (Dratio_ctrl)*gr_ctrl - (Dratio)*gr )/time ) - 1
    )
    
    inputData %<>% dplyr::mutate(
      #diff = abs(cell_count - cell_count__time0)/cell_count,
      #Dratio_times_gr = gr*Dratio,
      Eqidx = abs(cell_count - cell_count__time0)/cell_count < 1e-10
    )
    if(sum(inputData$Eqidx, na.rm = TRUE) > 0) {
      warning("issue when cell_count == cell_count__time0 (approx.), denominator of Dratio is approx. zero -> use L'hospital's rule on Dratio*gr")
      # a = cell_count
      # b = cell_count__time0
      # d = pmax(dead_count - dead_count__time0, 1, na.rm = T)
      # Dratio*gr = d*log2(a/b)/(a-b)  ===> consider f(x) = d*log2(x/b)/(x-b) as x->b
      # lim x->b f(x) = lim x->b { d*log2(x/b)/(x-b) } = d/(b*log(2)) by L'hospital's rule
      inputData %<>% dplyr::mutate(
        Dratio_gr = ifelse(!Eqidx, NA,
                           pmax(dead_count - dead_count__time0, 1, na.rm = T) /(cell_count__time0*log(2)) )
      )
      inputData %<>% dplyr::mutate(
        GR_s = ifelse(!Eqidx, GR_s,
                      2^((gr + Dratio_gr)/((1 + Dratio_ctrl)*gr_ctrl))-1),
        GR_d = ifelse(!Eqidx, GR_d,
                      2^(((Dratio_ctrl)*gr_ctrl - Dratio_gr)/time)-1)
      )
    }
    
    inputData %<>% dplyr::mutate(
      GR_combined = GR_s + GR_d,
      GR_naive = 2^(gr/gr_ctrl) - 1
    )
    inputData$log10_concentration = log10(inputData$concentration)
    tmp<-inputData[,groupingVariables, drop = FALSE]
    experimentNew = (apply(tmp,1, function(x) (paste(x,collapse=" "))))
    # if(cap == TRUE) {
    #   inputData$GRvalue[inputData$GRvalue > 1] = 1
    #   inputData$rel_cell_count[inputData$rel_cell_count > 1] = 1
    # }
    if(length(groupingVariables) > 0) {
      inputData$experiment = experimentNew
    } else {
      inputData$experiment = "All Data"
    }
    return(inputData)
  }
}