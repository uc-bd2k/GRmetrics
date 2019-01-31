
devtools::load_all()
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

fit = GRfit(test_input, groups[groups != "time"], cap = FALSE, case = "static_vs_toxic")
points = "average"
curves = "sigmoid"
plot_type = "static"
output_type = "separate"
fitData = fit
gg = GRdrawDRCV2(fit,
               points = "all",
               curves = c("sigmoid", "line", "none"),
               experiments = list(cell_line = c("HCC1806"), DrugName = c("Paclitaxel", "Abemaciclib/LY2835219") ),
               #experiments = list(cell_line = "HCC1806"),
               plot_type = c("static", "interactive"),
               #output_type = c("together", "separate")
               output_type = "separate"
)
