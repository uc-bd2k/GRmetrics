GRfitV2 = function(inputData, groupingVariables, cap = FALSE, case = "A",
                   initial_count) {
  #### some sort of input checking if necessary #######
  
  ##########
  
  #### send to .GRcalculateV2 function to calculate GR values
  gr_table = .GRcalculateV2(inputData, groupingVariables, cap, case,
                          initial_count)
  #### send results of this function to .GRlogisticFitV2 function to 
  #### fit logistic curves to GR values
  #GRlogfit = .GRlogisticFit(gr_table, groupingVariables, force, cap)
  
  #### package both tables and other necessary data into a list object
  #### decide on output format later... ditch SummarizedExperiment?
  output = list(
    #assays = GRlogfit,
    #colData = colData,
    #rowData = rowData, 
    metadata = list(gr_table = gr_table, 
                    groupingVariables = groupingVariables))
  return(output)
}
