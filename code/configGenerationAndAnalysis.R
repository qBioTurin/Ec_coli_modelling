
library(readr)

system("cd /home/riccardo/Documents/AntibioticResistance/")
setwd("/home/riccardo/Documents/AntibioticResistance/")

file.remove("Net/fluxes.trace")
file.remove("trace.log")
system("find -type f -name '*error.log*' -delete")

# model.path = "./Data/Escherichia_coli_str_K_12_substr_MG1655.mat"
# model.path = "./Data/Escherichia_coli_SE11.mat"
# model.path = "./Data/Escherichia_coli_O157_H7_str_Sakai.mat"
# model.path = "./Data/Escherichia_coli_UTI89_UPEC.mat"
model.path = "./Data/coli.xml"

# source("Rfunction/MatlabFileGeneration.R")
# source("Rfunction/sybilRoutine.R")
# source("Rfunction/FinalMatrixGeneration.R")

# source("Rfunction/FinalMatrixGenerationChangeObj.R")

system(paste0("unfolding2 ./Net/PNFBA -long-names "))
system(paste0("PN2ODE.sh ./Net/PNFBA -C ./Cpp/transitions.cpp -H "))

model_analysis(solver_fname = "Net/PNFBA.solver",
               functions_fname = "Rfunction/Functions.R",
               parameters_fname = "Input/ParametersList4.csv",
               f_time = 40,
               s_time = 0.1,
               solver_type = "LSODA",
               n_run = 1)

config <- readRDS("results_model_analysis/params_PNFBA-analysis.RDS")

setwd("./Net")

for(i in 1:length(config[["config"]]))
  write.table(x = config[["config"]][[i]][[1]][[3]], 
              file = config[["config"]][[i]][[1]][[1]],
              col.names = FALSE, row.names = FALSE, sep = ",")

cmd<-paste0(" .", .Platform$file.sep, basename(config$files$solver_fname), " ",
            config$out_fname,
            " -stime 0.1 ",
            " -itime 0 ",
            " -ftime 40 ",
            " -type LSODA  > ../trace.log")

system(cmd)

setwd("/home/riccardo/Documents/ESPNFBA/")

file.remove("dockerID")
file.remove("ExitStatusFile")

