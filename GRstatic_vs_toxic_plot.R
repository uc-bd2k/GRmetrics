.GRdrawSvT = function(inputData, groupingVariables) {
  #df = inputData %>% dplyr::group_by(CellLine, pert_type, DrugName)
  df = inputData %>%
    reshape2::melt(id.vars = c("CellLine", "pert_type", "DrugName", "Conc"),
                   measure.vars = c("GR_s", "GR_d", "GR_naive", "GR_combined")) %>%
    dplyr::filter(CellLine == "HCC1806") %>% 
    dplyr::select(CellLine, pert_type, DrugName, Conc, value, variable) %>%
    dplyr::group_by(CellLine, pert_type, DrugName, Conc, variable) %>%
    dplyr::summarize(value = mean(value, na.rm = T)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(exp = paste(CellLine, pert_type, DrugName)) %>%
    dplyr::filter(Conc > 0) %>%
    dplyr::filter(variable %in% c(GR_d, GR_s))
  g = ggplot2::ggplot(df) + 
    ggplot2::geom_point(aes(x = log10(Conc), y = value, colour = variable, group = exp)) + 
    #ggplot2::theme(legend.position = "none") + 
    ggplot2::facet_wrap(~exp)
  return(g)
}