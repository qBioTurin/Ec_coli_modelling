
library(sybil)

source("./Rfunction/findReactEq.R")

# Problem: S*x = b -. size x = n and size of S = m X n
# With this function the file to pass to the GLPKsolver is generated:
# 1 row)    n_row n_col GLP_MAX (or GLP_MIN if the obj has to max or min)
# 2 row)    the coeff for defining the obj function 
#           (the num of coeff has to be == length of x)
# m rows)   definition of the S row bounds (i.e. b)
# n rows)   definition of the x bounds
# m*n rows) definition of the S coeffs: row_index col_index value

load(file = "Data/FBAmodel.RData")

# extracting stoichiometric matrix
S = as.matrix(Matlab.file@S)
ncol=length(S[1, ])
nrow=length(S[ ,1])

obj = c(Matlab.file@obj_coef)
b = as.matrix(Matlab.file@obj_coef)*0 # equilibrium condition

# reactions and reactant's nomenclature in accordance with VMH
ReactionsNames = unlist(Matlab.file@react_id)
# [e] = extracellular metabolites
# [c] = cytosolic metabolites
ReagentsNames = unlist(Matlab.file@met_id)

# objective reaction
print(ReactionsNames[obj!=0])

# The current set of EX reacts in the model can be accessed findExchReact()
# rSummaryExch = findExchReact(Matlab.file)
# The subset of uptake reactions can be retrieved by method uptReact:
# upt <- sybil::uptReact(rSummaryExch)

# The current set of reacts in the model can be accessed findReact()
# rSummary = findReact(Matlab.file)
# findReactEq(Matlab.file, reagent = "lcts[e]")

# print reaction equation from reaction name:

# modelorg2tsv(Matlab.file, "BIGGdata", extMetFlag = "b", fielddelim = "\t", 
#              entrydelim = ", ", makeClosedNetwork = FALSE, onlyReactionList = FALSE,
#              minimalSet = TRUE, fpath = "./Data")

BIGG = read.table(file = './Data/BIGGdata_react.tsv', sep = '\t', header = TRUE)

# NB: react.names' elements must be in alphabetical order
react.names = c("EX_glc(e)", "EX_lcts(e)", "EX_o2(e)", "LACZ", "LACZpp", "LCTStpp")

# default bounds:
# change.lb = as.matrix(lowbnd(Matlab.file))[react_id(Matlab.file) %in% react.names]
# change.ub = as.matrix(uppbnd(Matlab.file))[react_id(Matlab.file) %in% react.names]

# prohibit excretion of carbon sources:
# change.lb = c(0, 0, -1000, 0, 0, -1000)
# make lactose primary carbon source:
# change.lb = c(0, -1000, -1000, 0, 0, -1000)
# cobra example:
# change.lb = c(-18.5, -1000, -1000, 0, 0, -1000)

# experiment:
change.lb = c(-6.5, -3, -15, 0, 0, -1000)
change.ub = c(1000, 1000, 1000, 1000, 1000, 1000)

# The subset of uptake reactions can be retrieved by method uptReact:
# upt = uptReact(rSummaryExch)

Matlab.file <- changeBounds(Matlab.file, 
                            checkReactId(Matlab.file, react.names)[react.names], 
                            lb = change.lb, ub = change.ub)

optL = optimizeProb(Matlab.file, algorithm = "fba",solver = "glpkAPI")

print(optL@lp_obj)

dfFluxEq = data.frame(react.names, 
                 lowbnd=lowbnd(Matlab.file)[react_id(Matlab.file) %in% react.names], 
                 uppbnd=uppbnd(Matlab.file)[react_id(Matlab.file) %in% react.names],
                 eq = BIGG$equation[c(react_id(Matlab.file) %in% react.names)],
                 flx=fluxes(optL)[match(react.names, react_id(Matlab.file))])

# l and u bounds:
lb = as.matrix(lowbnd(Matlab.file))
ub = as.matrix(uppbnd(Matlab.file))

rb <- cbind(b,b)
cb <- cbind(lb,ub)

write.csv2(dfFluxEq, './Data/dfFluxEq.csv')
# dfFluxEq <- read_delim("Data/dfFluxEq.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
save(ReactionsNames, file = "Data/ReactionsNames.RData")
save(lb, file = "Data/lbValues.RData")
save(ub, file = "Data/ubValues.RData")

