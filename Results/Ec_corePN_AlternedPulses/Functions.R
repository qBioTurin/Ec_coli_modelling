
init.gen <- function(iG, iL) {
  
  yini.names <- c("glc_e", "lcts_e")
  y_ini <- c(iG, iL)
  
  names(y_ini) = yini.names
  y_ini = y_ini[yini.names]
  
  return(y_ini)
  
}

MatrixGeneration <- function(frame){
  data = read.table(file = frame, sep = ',', header = F)
  return(data)
}