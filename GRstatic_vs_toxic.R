.GRcalculateV2 = function(inputData, groupingVariables, cap = FALSE, case = "A",
                        initial_count) {
  ### see lines 152-213 in https://github.com/datarail/DrugResponse/blob/1d0fd1898e797b74f611eb61dbd6926ac715bb83/MATLAB/import_columbus/Process_CellCountData2.m
  
  #inputData = readr::read_csv("/Users/nicholasclark/Desktop/Git/grmetrics_resources/GR_v2_all.csv")
  ### Note: what about NA values in Deadcount and Cellcount?? these will propagate currently, but we may want to set them to 0 instead?
  inputData %<>% 
    dplyr::mutate(too_few = Deadcount + Cellcount < .95*(Day0Cnt + Day0DeadCnt)) %>%
    dplyr::mutate(too_many = pert_type != 'ctl_vehicle' & 
        Deadcount + Cellcount > 1.15*(Ctrlcount + Ctrl_Deadcount))
  ## if missing cells, assign to dead cells
  inputData %<>% dplyr::mutate( Deadcount = ifelse(!too_few, Deadcount, 
    pmax(Deadcount, Day0Cnt + Day0DeadCnt - floor(.95*Cellcount) + 1, na.rm = T)))
  # if too many dead cells -> subtract dead cells
  inputData %<>% dplyr::mutate( Deadcount = ifelse(!too_many, Deadcount, 
    pmin(Deadcount, Ctrlcount + Ctrl_Deadcount - ceiling(1.15*Cellcount) - 1, na.rm = T)))
  
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
  
  # Transform time to days for normalization of GR_d #### change this to time = "days" or "hours" # default = not changing
  # if(max(inputData$Time, na.rm = TRUE) > 14) {
  #   inputData %<>% dplyr::mutate( Time = Time/24 )
  #   warning("Time transformed to days for normalization of GR_d")
  # }
  ### how should we treat NAs here?? na.rm or not?
  inputData %<>% dplyr::mutate(
    # Dratio = (increase in dead cells, maxed out at 100%)/(increase in alive cells)
    Dratio = pmax(Deadcount - Day0DeadCnt, 1, na.rm = T)/(Cellcount - Day0Cnt),
    Dratio_ctrl = pmax(Ctrl_Deadcount - Day0DeadCnt, 1, na.rm = T)/(Ctrlcount - Day0Cnt),
    gr = log2(Cellcount/Day0Cnt),
    gr_ctrl = log2(Ctrlcount/Day0Cnt)
    )
  ## These are slightly different from the formulas in the slide-deck, which are correct?
  inputData %<>% dplyr::mutate(
    GR_s = 2^( (1 + Dratio)*gr/( (1 + Dratio_ctrl)*gr_ctrl) ) - 1,
    GR_d = 2^( ( (Dratio_ctrl)*gr_ctrl - (Dratio)*gr )/Time ) - 1
  )

  inputData %<>% dplyr::mutate(
    #diff = abs(Cellcount - Day0Cnt)/Cellcount,
    #Dratio_times_gr = gr*Dratio,
    Eqidx = abs(Cellcount - Day0Cnt)/Cellcount < 1e-10
  )
  if(sum(inputData$Eqidx, na.rm = TRUE) > 0) {
    warning("issue when Cellcount == Day0Cnt (approx.), denominator of Dratio is approx. zero -> use L'hospital's rule on Dratio*gr")
    # a = Cellcount
    # b = Day0Cnt
    # d = pmax(Deadcount - Day0DeadCnt, 1, na.rm = T)
    # Dratio*gr = d*log2(a/b)/(a-b)  ===> consider f(x) = d*log2(x/b)/(x-b) as x->b
    # lim x->b f(x) = lim x->b { d*log2(x/b)/(x-b) } = d/(b*log(2)) by L'hospital's rule
    inputData %<>% dplyr::mutate(
      Dratio_gr = ifelse(!Eqidx, NA,
        pmax(Deadcount - Day0DeadCnt, 1, na.rm = T) /(Day0Cnt*log(2)) )
    )
    inputData %<>% dplyr::mutate(
      GR_s = ifelse(!Eqidx, GR_s,
                    2^((gr + Dratio_gr)/((1 + Dratio_ctrl)*gr_ctrl))-1),
      GR_d = ifelse(!Eqidx, GR_d,
                    2^(((Dratio_ctrl)*gr_ctrl - Dratio_gr)/Time)-1)
    )
  }
  
  inputData %<>% dplyr::mutate(
    GR_combined = GR_s + GR_d,
    GR_naive = 2^(gr/gr_ctrl) - 1
  )
  return(inputData)

  ##### taylor series approximation, translated from matlab code. 
  # if(sum(inputData$Eqidx, na.rm = TRUE) > 0) {
  #   warning("issue when Cellcount == Day0Cnt (approx.), denominator of Dratio is approx. zero -> Taylor series")
  #   ##### not sure about the correctness of this part
  #   # a = Cellcount
  #   # b = Day0Cnt
  #   inputData %<>% dplyr::mutate(
  #     Dratio_gr = #ifelse(!Eqidx, NA,
  #       pmax(Deadcount - Day0DeadCnt, 1, na.rm = T) *
  #         (1/Day0Cnt + rowSums( sapply(1:35, function(x) {
  #           Day0Cnt*x*((Day0Cnt-Cellcount)^x)
  #           # b*n*(b-a)^n
  #         }) ) ) #)
  #   )
  #   inputData %<>% dplyr::mutate(
  #     GR_s = ifelse(!Eqidx, GR_s,
  #                   2^((gr + Dratio_gr)/((1 + Dratio_ctrl)*gr_ctrl))-1),
  #     GR_d = ifelse(!Eqidx, GR_d,
  #                   2^(((Dratio_ctrl)*gr_ctrl - Dratio_gr)/Time)-1)
  #   )
  # }

}  