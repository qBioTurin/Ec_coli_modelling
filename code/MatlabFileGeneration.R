
library(sybilSBML)

# source("./Rfunction/readMATmod.R")
# Matlab.file <- readMATmod(model.path)

Matlab.file <- readSBMLmod(model.path)

save(Matlab.file, file = "Data/FBAmodel.RData")
