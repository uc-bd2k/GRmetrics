.rel_cell_logistic_3u = function(c, Einf, EC50, h){
  Einf + (1 - Einf)/(1 + (c/EC50)^h)
}