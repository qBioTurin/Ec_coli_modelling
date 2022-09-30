
library(epimod)
library(sybil)
library(readr)
library(sybilSBML)
library(stringr)

source("./Rfunction/readMATmod.R")
source("./Rfunction/findReactEq.R")

system("cd /home/riccardo/Documents/AntibioticResistance/")
setwd("/home/riccardo/Documents/AntibioticResistance/")

# Pietro Liò model:
model.path = "data/icdf834.mat"
fbamodel <- readMATmod(model.path)
# optL = optimizeProb(fbamodel, algorithm = "fba", solver = "glpkAPI")

genesets = fbamodel@gpr # gene-reaction association rule for each reaction

##################################################################
## Performing flux Variability Analysis (FVA) for a given model ##
##################################################################

## Rationale

# 1. Antibiotic lethality may be initiated by altered cellular redox state. 
# The redox stress potentially rewires central metabolism, cellular respiration, 
# and iron resulting in cellular damage

# 2. FVA was used to show perturbed central metabolism and redox balance in
# the presence of antibiotics and reprogrammed metabolism as compensatory 
# mechanisms in resistant populations.

## Methods

# 1. To understand the effect of antibiotics on C. difficile cellular metabolism 
# the changes in feasible metabolic flux distributions in the presence of 
# metronizadole (WT+mtz) was delineated using FVA and analysed.

# 2. The metabolic reprogramming due to antibiotic selection pressures was 
# analysed by FVA using customized models to represent WT, ChlR and StrpR
# Differences in unique forced fluxes/rates between resistant and susceptible
# populations in central metabolic pathways indicate compensatory metabolic
# reprogramming.

# 3. One can thus determine the minimum and maximum flux value that each 
# reaction in the model can take up while satisfying all constraints on 
# the system for a specific objective

## optimizeProb() attribute algorithm = "fv"
## Class "sysBiolAlg_fv"

## Further arguments passed to optimizer. 
## Argument algorithm is set to "fv", further arguments passed to the constructor 
## for objects of class sysBiolAlg_fv, see there for details

fv <- fluxVar(icdf834)
plot(fv)

# fv <- fluxVar(iCN900)
# plot(fv)
