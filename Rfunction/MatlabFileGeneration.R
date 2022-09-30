#install.packages("R.matlab")
library(R.matlab)
Matlab.file <- readMat(model)

save(Matlab.file,file = "Data/FBAmodel.RData")
