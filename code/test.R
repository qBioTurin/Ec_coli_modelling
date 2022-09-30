
source("./Rfunction/readMATmod.R")
NAP08 = readMATmod("./models/NAP08/Clostridium_difficile_NAP08.mat")

R.matlab::Matlab$startServer()
matlab <- R.matlab::Matlab()
isOpen <- open(matlab)

R.matlab::evaluate(matlab, "cd /home/riccardo/cobratoolbox")
R.matlab::evaluate(matlab, "initCobraToolbox")
R.matlab::evaluate(matlab, "cd /home/riccardo/Documents/FBAandPN/models/NAP08")
R.matlab::evaluate(matlab, "load('Clostridium_difficile_NAP08.mat');")

R.matlab::evaluate(matlab,"model = buildRxnGeneMat(model);")
R.matlab::evaluate(matlab, ("[minFlux, maxFlux] = fluxVariability(model);"))

minFlux = R.matlab::getVariable(matlab, "minFlux")
maxFlux = R.matlab::getVariable(matlab, "maxFlux")
minFlux = minFlux[[1]]
maxFlux = maxFlux[[1]]

# threshold for FVA 
threshold = 10e-06

for (i in 1:length(minFlux))  {if(abs(i) < threshold) {minFlux[i] = 0}}
for (i in 1:length(maxFlux))  {if(abs(i) < threshold) {maxFlux[i] = 0}}

R.matlab::evaluate(matlab, "writecell(model.rxnNames,'rxnNames_NAP08.xlsx','Sheet',1);")
R.matlab::evaluate(matlab, "writematrix(model.rxnGeneMat,'rxnGeneMat_NAP08.xlsx','Sheet',1);")

R.matlab::evaluate(matlab, "[ReactionPosition, genePosition] = find(model.rxnGeneMat);")
ReactionPosition = R.matlab::getVariable(matlab, "ReactionPosition")
ReactionPosition = ReactionPosition[[1]]
genePosition = R.matlab::getVariable(matlab, "genePosition")
genePosition = genePosition[[1]]

NAP08@mod_desc = "C_diffile"
NAP08@met_single = sybil::singletonMetabolites(NAP08, 
                                                   retIds = FALSE)[["smet"]]
NAP08@react_single <- sybil::singletonMetabolites(NAP08, 
                                                      retIds = FALSE)[["sreact"]]
NAP08@met_de <- sybil::deadEndMetabolites(NAP08, retIds = FALSE)[["dem"]]
NAP08@react_de <- sybil::deadEndMetabolites(NAP08, retIds = FALSE)[["der"]]

rxnGeneMat <- readxl::read_excel("models/NAP08/rxnGeneMat_NAP08.xlsx", col_names = FALSE)

NAP08@rxnGeneMat = Matrix::Matrix(FALSE, nrow = NAP08@react_num, ncol = 
                                        length(NAP08@allGenes), sparse = TRUE)

# print reaction equation from reaction name
sybil::modelorg2tsv(NAP08, paste("BIGGdata", NAP08@mod_desc, sep = "_"), 
                    extMetFlag = "b", fielddelim = "\t", 
                    entrydelim = ", ", makeClosedNetwork = FALSE, 
                    onlyReactionList = FALSE, minimalSet = TRUE,
                    fpath = "./models/NAP08/")

save(NAP08, file = "./models/NAP08/NAP08.RData")

# close the MATLAB server
close(matlab)
