
# PRIME: Phenotype prediction tool

# Phenotype of Regulatory influences Integrated with Metabolism and Environment (PRIME) 
# is a computational framework to rationally elucidate conditional response phenotypes 
# using integrated networks of regulation and metabolism. PRIME algorithm integrates 
# Environment and Gene Regulatory Influence Network (EGRIN) inferred gene regulatory
# network (GRN) to the curated genome-scale metabolic network (MN) for TF essentiality
# and phenotype predictions mechanistic inferences associated with the predicted TF 
# essentially and their network vulnerabilities can be assessed by this approach. 
# PRIME has substantially improved the condition-specific integration of statistically 
# inferred GRN with MN by using gene expression data and environment-specific weights.
# Specifically, a new 𝜸 factor (Reaction flux influencer) to quantify the magnitude 
# of the influence of a particular TF on genes in a MN has been incorporated, 
# thereby constraining the reaction flux to predict condition-specific phenotypes.

setwd("./Documents/AntibioticResistance")

# Libraries

library(R.matlab)

# Methods

source("./Rfunction/PRIME.R")
source("./Rfunction/readMATmod.R")

# Coding 

icdf838 = readMATmod("models/icdf838.mat")
MonoColonizedCdiff = readMATmod("models/MonoColonizedCdiff.mat")

transcriptomes = 
  read.csv("./modules/Cdifficile630_transcriptional_compendia_log2_ratios.csv",
           header=T, row.names = 1)

transcriptomes = as.matrix(transcriptomes)

transcriptome.metadata = 
  readxl::read_excel("./modules/Transcriptional_compendium_metadata.xlsx")

modules = read.delim("modules/biclusters_details_import.txt")
tf.infer = readr::read_csv("modules/tf_modules_inferelator.csv", 
                           show_col_types = T)

beta.cutoff = 0.075
tf.infer = tf.infer[tf.infer$Beta >= beta.cutoff | 
                      tf.infer$Beta <= -beta.cutoff, ]

Regulator = tf.infer$TF
Magnitude = tf.infer$Beta
Target = icdf838@allGenes

GEO.Sample = c("GSM458943", "GSM458942", "GSM458941", "GSM458940", "GSM458939", 
               "GSM458938", "GSM458937", "GSM458936", "GSM458935", "GSM458934",
               "GSM458933", "GSM458932")

GeneExpressionIDs = 
  rownames(transcriptomes[, colnames(transcriptomes) %in% GEO.Sample])
GeneExpression = 
  transcriptomes[, colnames(transcriptomes) %in% GEO.Sample][, 1] 

c(RegulatedModel, F_KO, TFs) = PRIME(model, Regulator, 
                                      Target, Magnitude, 
                                      GeneExpressionIDs, 
                                      GeneExpression)
