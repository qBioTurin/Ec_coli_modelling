
PRIME <- function(model, Regulator, Target, Magnitude, GeneExpressionIDs, 
                  GeneExpression)
  {
  
  # PRIME is an algorithm used to integrate regulatory network onto GSMMs
  #
  # It predicts the phenotype of TF and gene knockouts. 
  # It also gives the perturbed models as output.
  #
  # PRIME calculates the magnitude of TF influences on the genes and their 
  # associated reactions
  #
  #
  # USAGE:
      # c(RegulatedModel, F_KO, TFs) = PRIME (model, Regulator, 
      #                                       Target, Magnitude, 
      #                                       GeneExpressionIDs, 
      #                                       GeneExpression)
      # 
  # INPUT:
       #
       # model:                 Metabolic network as COBRA structure (with expression data 
       #                        applied using iMAT/GIMME) e.g., iEK1011
       # Regulator and Target:  Cell arrays listing the TF-target gene interactions. 
       #                        Regulators = {list of TFs}; Targets = {list of gene targets};
       #                        Regulators and Targets have the same number of rows; 
       #                        each row describes a separate interaction.
       # Magnitude:             Weights (Beta) from INFERELATOR run of EGRIN. 
       #                        Magnitude has the same number of rows as Regulators and Targets
       # GeneExpressionIDs:     Gene ID of GE data (present study) e.g., isoniazid treatment
       # Gene Expression:       GE data corresponding to the GeneExpressionIDs (present study)
       #
  # OUTPUT:
       # 
       # RegulatedModel:        TF pertubed model (Models of TF KOs)
       # TFs:                   List of TFs being knocked out or over expressed
  
  threshold = 10e-06
  uniqueTFs = unique(Regulator)
  
  Matlab$startServer()
  matlab <- Matlab()
  isOpen <- open(matlab)
  
  if (!isOpen) throw("MATLAB server is not running: waited 30 seconds.")
  
  evaluate(matlab, "cd /home/riccardo/cobratoolbox")
  evaluate(matlab, "initCobraToolbox")
  evaluate(matlab, "cd /home/riccardo/Documents/AntibioticResistance/models")
  evaluate(matlab, "load('icdf838.mat');")
  
  evaluate(matlab,"icdf838 = buildRxnGeneMat(icdf838)")
  evaluate(matlab, ("[minFlux, maxFlux] = fluxVariability(icdf838)"))
  
  minFlux = getVariable(matlab, "minFlux")
  maxFlux = getVariable(matlab, "maxFlux")
  
  # default minFlux
  # minFlux = model.lb
  # default maxFlux
  # maxFlux = model.ub
  
  minFlux = minFlux[[1]]
  maxFlux = maxFlux[[1]]
  
  for (i in 1:length(minFlux))  {
    if(abs(i) < threshold) {
      minFlux[i] = 0
    }
  }
  
  for (i in 1:length(maxFlux))  {
    if(abs(i) < threshold) {
      maxFlux[i] = 0
    }
  }
  
  # rank the magnitude value to compare the models across conditions
  Magnitude = (Magnitude-min(Magnitude))/(max(Magnitude)-min(Magnitude))
  
  RxnNumber = dim(icdf838@S)[2]
  TFs = uniqueTFs
  genelist = Target %in% icdf838@allGenes
  nGenes = length(icdf838@allGenes)
  nRxns = icdf838@react_num
  print("calculating TF influences...")
  
  fba = sybil::optimizeProb(icdf838, algorithm = "fba",solver = "glpkAPI")
  WTGr = fba@lp_obj
  fbaFlux = fba@fluxdist@fluxes@x
  Refine = replicate(RxnNumber, 1)
  
  evaluate(matlab, "[ReactionPosition, genePosition] = find(icdf838.rxnGeneMat);")
  ReactionPosition = getVariable(matlab, "ReactionPosition")
  ReactionPosition = ReactionPosition[[1]]
  genePosition = getVariable(matlab, "genePosition")
  genePosition = genePosition[[1]]
  
  sumgeneBetas = replicate(nGenes, 1)
  
  for (n in 1:length(icdf838@allGenes)){
    genePos = Target %in% icdf838@allGenes[n]
    geneMagnitude = Magnitude[genePos]
    sumgeneBetas[n] = sum(abs(geneMagnitude))
  }
  
  print('mapping ReFInE to maxFlux of Reactions...')
  print('running Flux Balance Ananlysis on the perturbed model')
  
  for (m in 1:length(TFs)){
    
    TFposition = Regulator %in% TFs[m]
    Targets = Target[TFposition]
    TargetMagnitude = Magnitude[TFposition]
    nReactions = length(ReactionPosition)
    calculatedBetas = replicate(nReactions, 1)
    Factor = replicate(nReactions, 1)
    FactorMultiplied = replicate(nReactions, 1)
    GEdata = GeneExpression
    TFpositionGE = GeneExpressionIDs %in% TFs[m]
    TFinGE = which(TFpositionGE)
    
    TFactivity = c()
    
    if (is.null(TFinGE) == F) {
      TFactivity[m] = GEdata[TFpositionGE]
    } else {
      TFactivity[m] = min(!is.na(GEdata))
    }
    
    for (i in 1:length(ReactionPosition)){
      gene = icdf838@allGenes[genePosition[i]]
      geneID = which(Targets %in% gene)
      geneIDmodel = which(icdf838@allGenes %in% gene)
      if (is.null(geneID) == T) {
        Beta = TargetMagnitude[geneID]
        calculatedBetas[i] = 1 - ((Beta/sumgeneBetas[geneIDmodel]) * TFactivity[m])
      }
      if(calculatedBetas[i] == 0) {
        calculatedBetas[i] = 1-(Beta*TFactivity[m])
      }
    }
    
    Factor = calculatedBetas
    # all Beta associated to the reaction
    FactorMultiplied = FactorMultiplied*Factor
    GPRreactions = unique(ReactionPosition)
    Reactions = icdf838@react_id[GPRreactions]
    
    for  (k in 1:length(GPRreactions)) {
      GPRreactionPosition = ReactionPosition %in% GPRreactions[k]
      ReactionIndex = which(GPRreactionPosition)
      # consider the min of all genes involved in the reaction
      FinalFactor = min(FactorMultiplied[ReactionIndex])
      ModelReactionIndex = icdf838@react_id %in% Reactions[k]
      ModelReactionPosition = which(ModelReactionIndex)
      # vector of Refine; assigned to the reaction position
      Refine[ModelReactionPosition] = FinalFactor
    }
    
    inhibit = numeric(length(icdf838@react_id))
    # new weight factor
    c = inhibit*min(abs(Refine))
    
    sybil::changeBounds(icdf838, 
                 checkReactId(icdf838, icdf838@react_id)[icdf838@react_id], 
                 lb = as.vector(minFlux * Refine), 
                 ub = as.vector(maxFlux * Refine))
    
    ubound = icdf838@uppbnd
    lbound = icdf838@lowbnd
    
    for (x in 1:RxnNumber) {
      # maintain stoichiometry consistency to be lb<ub
      ubound[x] = max(ubound[x],fbaFlux[x])
      # maintain stoichiometry consistency to be lb<ub
      lbound[x] = min(lbound[x],fbaFlux[x])
    }
    
    sybil::changeBounds(icdf838, 
                        checkReactId(icdf838, icdf838@react_id)[icdf838@react_id], 
                        lb = lbound, 
                        ub = ubound)
    
    FBAsoln = sybil::optimizeProb(icdf838, algorithm = "fba",solver = "glpkAPI")
  }
  
  # close the MATLAB server
  close(matlab)
  
  return(icdf838, FBAsoln, TFs)
}