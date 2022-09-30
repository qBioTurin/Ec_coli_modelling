
ExchangesFromReactions = function(model, target.react, ex) {
  
  #' It saves the output of the RStudio console to a log file in the R programming language
  #' rep
  #' 
  #' USAGE:
  #' 
  #' ExchangesFromReactions(model, target.react, ex)
  #' 
  #' INPUTS:
  #' 
  #' model:             A single character string indicating model name
  #' target.react:      A single character string indicating the target reactions
  #' ex:                An object of class reactId_Exch
  #
  # OUTPUTS:
  #     List of Ex reactions stored as console output "txt"
  
  # loading C. difficile metabolic model
  load(paste("./data/", model, "/", model, ".RData", sep = ""))
  
  ex = sybil::findExchReact(model.mat)
  
  BIGGdata__react =
    readr::read_delim(paste("data/", 
                            model, "/BIGGdata_", model, "_react.tsv", sep = ""), 
                      delim = "\t", escape_double = FALSE, 
                      trim_ws = TRUE)
  
  # open Matlab server
  R.matlab::Matlab$startServer()
  matlab <- R.matlab::Matlab()
  isOpen <- open(matlab)
  
  # Importing model into the MATLAB environment
  R.matlab::evaluate(matlab, "cd /home/riccardo/cobratoolbox")
  R.matlab::evaluate(matlab, "initCobraToolbox")
  R.matlab::evaluate(matlab, paste("cd /home/riccardo/Documents/FBAandPN/data/", model, sep = ""))
  R.matlab::evaluate(matlab, paste("model = importdata('", model, ".mat');", sep = ""))
  
  # File name of output log
  ExchangesFromReactions <- file(paste("./data/", model, 
                                       "/ExchangesFrom", 
                                       target.react, ".txt", sep = ""))
  
  
  R.matlab::evaluate(matlab, 
                     (paste0("[metList, stoichiometries] = findMetsFromRxns(model, ",
                             paste("'"), target.react, paste("'"), paste(");"), 
                             sep = "")))
  
  R.matlab::evaluate(matlab, ("list = metList{1, 1};"))
  
  list = R.matlab::getVariable(matlab, "list")
  list = gsub("\\c]", "\\e]", as.character(unlist(list)))
  
  # Writing console output to log file
  sink(ExchangesFromReactions, append = TRUE, type = "output")
  sink(ExchangesFromReactions, append = TRUE, type = "message")
  
  cat(" Metabolites involved in", paste(target.react), ":")
  cat("\n")
  print(list)

  
  for(i in 1:length(list)) {
    cat("Exchanges for", paste(list[i], "= "))
    (print(ex@react_id[which(ex@met_id == list[i])]))
  }
  
  # close the MATLAB server
  close(matlab)
  
  # Close connection to log file
  closeAllConnections() 
  
  file.remove("./InputStreamByteWrapper.class")
  file.remove("./MatlabServer.m")
  
}