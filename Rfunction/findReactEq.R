findReactEq <- function(model, reagent) {
  
  # Returns a list of equations (formulas) of reactions in which the 
  # metabolite (reagent)is involved
  #
  # USAGE:
  #
  #   findReactEq(model, reagent)
  #
  # INPUTS:
  #    model:             Model structure
  #    metList:           Metabolite
  #
  # OUTPUTS:
  #     List of reactions + Reaction formulas corresponding

  ReactionsNames = unlist(Matlab.file@react_id)
  ReagentsNames = unlist(Matlab.file@met_id)
  S = as.matrix(Matlab.file@S)
  
  ReagentsIndex = which(ReagentsNames == reagent)
  
  entries = S[ReagentsIndex, ]
  
  ReactionIndex = which(entries!=0)
  reactList = ReactionsNames[ReactionIndex]
  
  sybil::printReaction(model, react = reactList, printOut = TRUE)
}
