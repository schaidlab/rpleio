## internal function for pleiotropy functions, center a quantitative variable
## Author: DJ Schaid
## Date: 7/27/2016

center <- function(y){
  mn <- apply(y, 2, mean)
  y.centered <- t(t(y) - mn)
}
