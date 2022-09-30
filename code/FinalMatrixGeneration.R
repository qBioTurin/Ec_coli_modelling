
library(foreach)
library(doParallel)

# write(paste0(nrow," ; ",ncol," ; ","GLP_MAX" ),file = fileName )
# write(paste(obj,collapse = " "),file = fileName,append = T )

FinalMatrix <- paste0(nrow," ; ",ncol," ; ","GLP_MAX" )
FinalMatrix <- c(FinalMatrix, paste(obj, collapse = " "))

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

save(FinalMatrix, file = "Data/FinalMatrix.RData")
