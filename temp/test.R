library(magrittr)

devtools::load_all()
inputData = inputCaseA
groupingVariables = c("cell_line", "agent")
case = "A"
force = FALSE
cap = FALSE
initial_count = T
inputData = .GRcalculate(inputData, groupingVariables, cap, case,
                        initial_count)

drc_output = GRfit(inputData, groupingVariables = c('cell_line','agent'))
fitData = drc_output
metric = c("GR", "rel_cell")
experiments = "all"
color = "experiment"
points = c("average", "all", "none")
curves = c("sigmoid", "line", "biphasic", "none")
bars = c("none", "sd", "se")
xrug = c("none", "GR50", "GEC50", "IC50", "EC50")
yrug = c("none", "GRinf", "GRmax", "Einf", "Emax")
theme = c("classic", "minimal", "bw")
palette = c("default","npg", "aaas")
facet = "none"
plot_type = c("static", "interactive")
min = "auto"; max = "auto"
