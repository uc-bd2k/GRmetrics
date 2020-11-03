.convert_new_grp = function(inputData, case, initial_count, groupingVariables) {
  ## if case == "A", don't do anything
  if (case == "A") { return(inputData) }
  ## accept "B" or "C" for long-format case
  else if (case %in% c("B","C")) {
    ### for each combination of grouping variables
    all_keys = groupingVariables %>% union(c("time", "concentration"))
    ### warn about ignored columns
    delete_cols = colnames(inputData) %>% setdiff(c(all_keys, "cell_count"))
    if(length(delete_cols) > 0) {message(paste0("The following columns are not specified as grouping variables and will be ignored: ", paste0(delete_cols, collapse = ", ")))}
    ### Note: concentration should be zero and treatment = "-" for all control cell counts
    ctrl_keys = groupingVariables %>% union(c("time")) %>% setdiff("treatment")
    ### Note: concentration and time should be zero and treatment = "-" for all time zero cell counts
    time0_keys = groupingVariables %>% setdiff("treatment")
    ### get control data and average over replicates
    ctrl_df = inputData %>% dplyr::filter(concentration == 0, time > 0) %>%
      dplyr::group_by(across(all_of(ctrl_keys)))
    ctrl_df_avg = ctrl_df %>% dplyr::select(-treatment) %>%
      dplyr::summarise(cell_count__ctrl = GRmetrics:::.trim_mean(cell_count, percent = 50),
                       .groups = "keep")
    ### get time zero data and average over replicates
    time0_df = inputData %>% dplyr::filter(time == 0) %>%
      dplyr::group_by(across(all_of(time0_keys)))
    time0_df_avg = time0_df %>% 
      dplyr::summarise(cell_count__time0 = GRmetrics:::.trim_mean(cell_count, percent = 50),
                       .groups = "keep")
    ### get treated (non-control, non-time-zero) data
    trt_df = inputData %>% dplyr::filter(time > 0, concentration > 0)
    ### un-comment next line to delete extra columns not specified as grouping variables
    #trt_df = trt_df %>% dplyr::select(-all_of(delete_cols))
    ### join data frames
    caseA_output = trt_df %>% dplyr::full_join(ctrl_df_avg, by = ctrl_keys) %>%
      dplyr::full_join(time0_df_avg, by = time0_keys) %>%
      dplyr::arrange(across(all_of(ctrl_keys)))
    ### delete cell_count__time0 column if case initial_count == FALSE
    if (!initial_count) { caseA_output$cell_count__time0 = NULL }
    ### convert from tibble to data frame and return the output
    caseA_output = as.data.frame(caseA_output)
    row.names(caseA_output) = 1:dim(caseA_output)[1]
    return(caseA_output)
  }
}
