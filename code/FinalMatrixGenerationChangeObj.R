
library(foreach)
library(doParallel)

# write(paste0(nrow," ; ",ncol," ; ","GLP_MAX" ),file = fileName )
# write(paste(obj,collapse = " "),file = fileName,append = T )

FinalMatrix <- paste0(nrow," ; ",ncol," ; ","GLP_MAX" )

# The objective  flux  needs  to  be  defined,  which  is  optimized, hence
# minimized or maximized. Intuitively, this objective of optimization must meet 
# cellular  purposes. Most  often,  biomass  production  is chosen  as  objective, 
# since  they  strive  for  optimal  growth. However,  also  other  objectives  
# might  be  chosen. For  example  for  finding  optimal  sets  of  reactions for
# energy efficiency such as minimization of ATP consumption [doi: 10.1038/msb4100162]

# sets the objective function to the ATP maintenance reaction:

Matlab.file <- changeObjFunc(Matlab.file, c("ATPM", "BIOMASS"), c(1, 1))
printObjFunc(Matlab.file)
obj_new = c(Matlab.file@obj_coef)

FinalMatrix <- paste0(nrow," ; ",ncol," ; ","GLP_MAX" )
FinalMatrix <- c(FinalMatrix, paste(obj_new, collapse = " "))

for(i in 1:nrow){
  #write(paste0("GLP_F"," ; ",paste0(rb[i,],collapse = " ; ") ),file = fileName ,append = T)
  FinalMatrix <- c(FinalMatrix, paste0("GLP_FX"," ; ", paste0(rb[i, ], collapse = " ; ") ))
}

for(j in 1:ncol){
  #write(paste0("GLP_DB"," ; ",paste0(cb[j,],collapse = " ; ")),file = fileName ,append = T)
  if( diff(cb[j, ]) == 0 )
    FinalMatrix <- c(FinalMatrix, paste0("GLP_FX"," ; ", paste0(cb[j, ], collapse = " ; ") ))
  else
    FinalMatrix <- c(FinalMatrix,paste0("GLP_DB"," ; ",paste0(cb[j,],collapse = " ; ")) )
}

cores=detectCores()
# not to overload the machine
# cl <- makeCluster(cores[1]-1)
cl <- makeCluster(cores)
registerDoParallel(cl)

for(i in 1:nrow){
  print(i)
  Matrix <- foreach(j=1:ncol, .combine='c') %dopar% {
    paste0(i," ; ", j, " ; ", S[i,j])
  }
  FinalMatrix <- c(FinalMatrix, Matrix)
}

#  for(i in 1:nrow){
#    for(j in 1:ncol){
# write(paste0(i," ; ", j, " ; ", S[i,j]),file = fileName ,append = T)
#      FinalMatrix <- c(FinalMatrix, paste0(i," ; ", j, " ; ", S[i,j]))
#    }
#  }

stopCluster(cl)

save(FinalMatrix, file = "Data/FinalMatrix_ATPM.RData")
