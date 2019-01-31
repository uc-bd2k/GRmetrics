test = read.csv("resources/GR_v2_all_from_R.csv")
colnames(test)[colnames(test) == "Day0DeadCnt"] = "dead_count__time0"
colnames(test)[colnames(test) == "Day0Cnt"] = "cell_count__time0"
colnames(test)[colnames(test) == "Cellcount"] = "cell_count"
colnames(test)[colnames(test) == "Deadcount"] = "dead_count"
colnames(test)[colnames(test) == "Ctrl_Deadcount"] = "dead_count__ctrl"
colnames(test)[colnames(test) == "Ctrlcount"] = "cell_count__ctrl"
colnames(test)[colnames(test) == "Conc"] = "concentration"
colnames(test)[colnames(test) == "Time"] = "time"
colnames(test)[colnames(test) == "CellLine"] = "cell_line"

counts = c("dead_count__time0", "cell_count__time0", "cell_count", "dead_count", "dead_count__ctrl", "cell_count__ctrl")
groups = c("cell_line", "time", "DesignNumber", "DrugName", "pert_type")
test_input = test[, colnames(test) %in% c(counts, groups, "concentration")]

fit = GRfitV2(test_input, groups[groups != "time"], cap = FALSE, case = "static_vs_toxic",
  initial_count = T)
GRdrawDRCV2(fit, points = "average",
           experiments = list(cell_line = c("HCC1806"), DrugName = "Paclitaxel"),
           plot_type = c("static", "interactive"),
           output_type = c("together", "separate") )
fit = .GRcalculate(test_input, groups[groups != "time"], cap = FALSE, case = "static_vs_toxic",
                   initial_count = T)
fitData$metadata$gr_table = fit
inputData = fit
fit = GRfit(test_input, groups[groups != "time"], cap = FALSE, case = "static_vs_toxic")


inputData = fit$metadata$gr_table
groupingVariables = fit$metadata$groupingVariables
