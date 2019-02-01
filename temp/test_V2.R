
# devtools::load_all()
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
colnames(test)[colnames(test) == "DrugName"] = "treatment"
colnames(test)[colnames(test) == "DesignNumber"] = "design"



counts = c("dead_count__time0", "cell_count__time0", "cell_count", "dead_count", "dead_count__ctrl", "cell_count__ctrl")
groups = c("cell_line", "time", "design", "treatment", "pert_type")
test_input = test[, colnames(test) %in% c(counts, groups, "concentration")]

test_input = test_input[,c("cell_line", "treatment", "pert_type", "design", "time",
                           "concentration", "cell_count__time0", "dead_count__time0",
                           "cell_count", "dead_count", "cell_count__ctrl", "dead_count__ctrl")]
test_small = test_input %>% dplyr::filter(cell_line == "BT20",
                                          concentration > 0)

test_med = test_input[1:50000,]
write.csv(test_med, file = "temp/gr_static_vs_toxic_input_med.csv", quote = T, row.names = F )
write.csv(test_small, file = "temp/gr_static_vs_toxic_input_small.csv", quote = T, row.names = F)

write.csv(test_input, file = "temp/gr_static_vs_toxic_input_edited.csv", quote = T, row.names = F)
test_small2 = read_csv("temp/gr_static_vs_toxic_input_small.csv")

test_input = read.csv("temp/gr_static_vs_toxic_input.csv")

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
