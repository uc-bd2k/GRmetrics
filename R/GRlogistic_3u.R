.GRlogistic_3u = function(c, GRinf, GEC50, h_GR){
  GRinf + (1 - GRinf)/(1 + (c/GEC50)^h_GR)
}