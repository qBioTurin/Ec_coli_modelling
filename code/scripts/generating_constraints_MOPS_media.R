
###################################################

# FBA allows estimation of intracellular reaction rates using organism-specific GSMMs.
# This methodology employing complex media characterized by substrates uptake and products secretion.
# Metabolimic measurements obtained for a set of metabolites can being converted to exchange rates 
# and then used as constraints for FBA

# We are considering E. coli strain K-12 MG1655 cells growing aerobically in MOPS minimal medium 
# supplemted with 0.1% glucose as caborn supply.

# The E. coli Metabolome Database (ECMDB) is a database containing quantitative, 
# analytic or molecular-scale information about E. coli (strain K-12) metabolites and 
# their associated properties, pathways, reactions functions, sources, enzymes or transporters

# What does aerobic growth in MOPS minimal medium supplemented with 0.1% glucose mean?

# - *Aerobic growth: This refers to the growth of organisms in the presence of oxygen. 
# E. coli strain K-12 MG1655, use oxygen to efficiently generate energy through aerobic respiration.

# - *MOPS minimal medium*: This is a type of culture medium used for growing bacteria. 
# It’s called “minimal” because it contains only the essential nutrients that the bacteria need to grow. 
# MOPS is a buffering agent that helps maintain the pH of the medium.

# - *Supplemented with 0.1% glucose*: This means that glucose has been added to the medium
# as a carbon source, which the bacteria can use for energy and growth. 
# The concentration of glucose is 0.1%, or 1 gram per liter.

##### SETTING MODELLIG CONDITIONS #####
#
#    1) We could set 1g/L of glc_e initial concentration
#    2) We could define a dry mass weight of a E. coli cell cellular population
#    3) We could image to model the bacterial culture growing in bioreactor with a working volume of about 3 L
#
##### Setting initial concentration #####
#
# In a study on the growth of recombinant E.coli, the optimal concentration of glucose was found to be 1g/L
# DOI:10.4028/www.scientific.net/AMR.647.185
#
# (1) to convert 1g/L of glucose to mmol/L, we need to know the molar mass of glc:

conc_glc_e_gL = 1 # (g/L)
mm_glc = 180.16 # (g/mol)

# then can then apply the following relation:

conc_glc_e_mmolL = (conc_glc_e_gL/mm_glc)*1000 # (mmol/L)

##### Setting bacterial biomass #####
#
# To estimate the number of cells that can constitute 1 g of bacterial dry mass weight
# we can divide 1 g by the dry mass weight of a single E. coli cell.
# 
# (2) estimating the the dry mass weight of a E. coli cellular population in batch
# 
# from -> https://ecmdb.ca/e_coli_stats

gDW_ec = 3e-13 # (g/cell)
n_ec_cell = 1/gDW_ec # (cells)

# Given data
doubling_time <- 1.0583 # (h) 63.5 (minutes)

# Define the molecular weight of biomass (MW)
MW_biomass <- 1  # g mmol^(-1)

# Calculate specific growth rate (μ)
growth_rate <- log(2) / doubling_time

# Calculate biomass flux (mmol/gDW*h)
biomass_flux <- growth_rate * MW_biomass

# Print the result
# cat("Biomass flux (mmol/gDW*h):", biomass_flux, "\n")

##### Setting batch culture volume and glucose amount ##### 
#
# A saturated cell culture of E. coli contains about 1e+09 cell/mL

n_ec_cell_log_mid <- 0.5*1e+09 # (cells/mL)

# Minimum volume that can contain these cells
v_batch = n_ec_cell / n_ec_cell_log_mid # (mL)

# (3) setting glucose amount given a batch culture volume of about "v_batch"

iG <<- (conc_glc_e_mmolL/1000)*v_batch # (mmol)

#######################################################
#######################################################

# loading BiGG models metabolites annotation
#
# bigg_models_metabolites <- read.delim(file.path(wd, "input/gene_expression/bigg_models_metabolites.txt"))
# 
# bigg_met_iML1515 <- bigg_models_metabolites[bigg_models_metabolites$bigg_id %in% model@met_id, ]
# not_bigg_id <- as.character(model@met_id[!model@met_id %in% bigg_models_metabolites$bigg_id])

############################ Resources ###################################à

# (1) From -> https://mediadb.systemsbiology.net/defined_media/media/102/

# MOPS minimal medium (bergholz et al)
# 1 Source(s): Bergholz et al, 2007

MOPS_mediadb = read_delim(paste0(wd, "/input/gene_expression/MOPS_mediadb.txt"), 
                           delim = "\t", escape_double = FALSE, trim_ws = TRUE)

MOPS_mediadb$amount_batch = (MOPS_mediadb$amount_mM/1000)*v_batch # (mmol)
MOPS_mediadb$bigg_id = c("glc_D_e", "mg2_e", "ca2_e", NA, 
                         "k_e", NA, "na1_e", "nh4_e", "fe3_e", 
                         "cobalt2_e", "zn2_e", "mn2_e", "cu2_e",
                         "mobd_e", NA, NA)

MOPS_mediadb = na.omit(MOPS_mediadb)

MOPS_mediadb$exchanges = paste0("EX_", MOPS_mediadb$bigg_id)

# The units “mmol/L” and “mM” are actually equivalent

############################

# (2) From -> https://ecmdb.ca/concentrations?c=compounds.met_id&d=down

# (1) Sajed, T., Marcu, A., Ramirez, M., Pon, A., Guo, A., Knox, C., Wilson, M., Grant, J., Djoumbou, Y. and Wishart, D. (2015). ECMDB 2.0: A richer resource for understanding the biochemistry of E. coli. Nucleic Acids Res, p.gkv1060. 26481353 .
# (2) ECMDB: The E. coli Metabolome Database. Guo AC, Jewison T, Wilson M, Liu Y, Knox C, Djoumbou Y, Lo P, Mandal R, Krishnamurthy R, Wishart DS. Nucleic Acids Res. 2012 Jan;41(Database issue):D625-30. 23109553 .
