
# loading libraries
library(devtools)
# install_github("https://github.com/qBioTurin/epimod", ref="master", force = T)
library(epimod)
# downloadContainers(tag = "latest")
library(tidyr)        # load the dplyr package for data manipulation
library(ggplot2)      # load the ggplot2 package for data visualization
library(patchwork)    # load the patchwork package for combining plots
library(RColorBrewer) # load the RColorBrewer package for color palettes
library(parallel)
#   Data plots for selected GEO samples
library(GEOquery)
library(limma)
library(umap)
library(readr)
library(stringr)