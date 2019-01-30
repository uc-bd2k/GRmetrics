GRdrawDRCV2 = function(fitData, points = c("average", "all", "none"),
                       experiments = "all",
                       plot_type = c("static", "interactive"),
                       output_type = c("together", "separate"),
                       ) {
  #### check input #####
  # make all inputs length 1
  #metric = metric[1]
  points = points[1]
  # curves = curves[1]
  # bars = bars[1]
  # xrug = xrug[1]
  # yrug = yrug[1]
  # theme = theme[1]
  # palette = palette[1]
  # plot_type = plot_type[1]
  # output_type = output_type[1]
  #########
  
  fit_df = fitData$metadata$gr_table
  # get grouping variables
  fit_groups = GRmetrics::GRgetGroupVars(fitData)
  #fit_groups = fitData$metadata$groupingVariables
  id = c(fit_groups, "Time", "Conc")
  fit_df_melt = fit_df %>%
    reshape2::melt(id.vars = id,
                   measure.vars = c("GR_s", "GR_d", "GR_naive", "GR_combined"),
                   value.name = "GRvalue", variable.name = "GR_metric"
                    )
  ### what should I do with GR values with concentration zero?
  #ctrl_df = fit_df_melt %>% dplyr::filter(Conc == 0)
  grp = syms(c(id, "GR_metric"))
  exp_sym = syms(c(fit_groups, "Time"))
  if(points == "average") {
    df = fit_df_melt %>%
      dplyr::filter(CellLine == "HCC1806") %>% 
      #dplyr::select(CellLine, pert_type, DrugName, Conc, value, variable) %>%
      dplyr::group_by(!!!grp) %>%
      dplyr::summarize(GRvalue = mean(GRvalue, na.rm = T)) %>%
      dplyr::ungroup() 
  } else if(points == "all") {
    df = fit_df_melt
  } else if(points == "none") {
    ### fill in later
  }
  df %<>%
    #dplyr::filter(CellLine == "HCC1806") %>% 
    #dplyr::select(CellLine, pert_type, DrugName, Conc, value, variable) %>%
    dplyr::mutate(exp = paste(!!!exp_sym)) %>%
    dplyr::filter(Conc > 0) #%>%
    #dplyr::filter(variable %in% c("GR_d", "GR_s"))
  n_plots = length(unique(df$exp))
  if(n_plots > 25) warning("Too many plots [change warning message later]")
  g = ggplot2::ggplot(df) + 
    ggplot2::geom_point(aes(x = log10(Conc), y = GRvalue, colour = GR_metric, group = exp)) + 
    #ggplot2::theme(legend.position = "none") + 
    ggplot2::facet_wrap(~exp)
  return(g)
}