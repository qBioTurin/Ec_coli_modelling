#' Find biomass reaction in model
#' @description Helper function to search for biomass reaction in available reactions of a model
#' 
#' @param model Object of class sybil::modelorg containging the genome sclae metabolic model
#' @param keys Vector with strings which are used to find biomass reaction in model
#' @return Vector with reaction ids for biomass reaction(s)

findrBiomass <- function(model, keys=c("biom", "cpd11416")){
  ex     <- sybil::findExchReact(model)
  ex_pos <- ex@react_pos
  ex_biom<- c(grep(paste0(keys, collapse = "|"), ex@met_id, ignore.case = TRUE),
              grep(paste0(keys, collapse = "|"), sybil::met_name(model)[ex@met_pos], ignore.case = TRUE))
  rbio <- vector()
  for(k in keys){
    idx <- grep(k, sybil::react_id(model), ignore.case = TRUE)
    if(length(idx)==0) idx <- grep(k, sybil::react_name(model), ignore.case = TRUE)
    if(length(idx)>0) rbio <- c(rbio, idx)
  }
  if( length(ex_biom) > 0 ){ # take care of biomass metabolite
    rbio   <- c(rbio, ex@react_pos[ex_biom])
    ex_pos <- setdiff(ex_pos, ex@react_pos[ex_biom]) 
  } 
  if(length(rbio)==0) return(NULL)
  
  rbio <- setdiff(rbio, ex_pos) # exclude exchange reactions
  return(sybil::react_id(model)[rbio])
}