
################################################################

# Gene expression data description
#
# Platform organism - Escherichia coli K-12
# Sample organism - Escherichia coli
# 
# The study investigated the effect of different carbon sources on E. coli global gene expression. 
# We grew MG1655 cells aerobically in MOPS minimal medium with either: 
# - glucose
# - glycerol
# - succinate
# - L-alanine
# - acetate
# - L-proline 
# as the carbon supply. 
# Samples were taken from each culture at mid log phase 
# and were RNA-stabilized using Qiagen RNAProtect Bacterial Reagent (Qiagen). 
# Total RNA was then isolated using MasterPure kits (Epicentre Technologies). 
# Purified RNA was reverse-transcribed to cDNA, labeled and hybridized to Affymetrix GeneChip.
# E. coli Antisense Genome Arrays as recommended in the technical manual (www.affymetrix.com).

# Samples (15)
# GSM37063 	glucose1
# GSM37064 	glucose2
# GSM37065 	glucose3
# GSM37066 	glucose4
# GSM37067 	glucose5
# GSM37068 	glycerol1
# GSM37069 	glycerol2
# GSM37070 	succinate1
# GSM37071 	succinate2
# GSM37072 	alanine1
# GSM37073 	alanine2
# GSM37074 	acetate1
# GSM37075 	acetate2
# GSM37076 	proline1
# GSM37077 	proline2

################################################################

model_type = "ecoli"
model_name = "iML1515"

# Load gene expression data from the GEO database for dataset GSE2037.
gset <- getGEO("GSE2037", GSEMatrix = TRUE, getGPL = FALSE)

# Check the number of datasets retrieved. If more than one, find the one with "GPL199" in its metadata.
if (length(gset) > 1) {
  idx <- grep("GPL199", attr(gset, "names"))
} else {
  idx <- 1  # If there's only one dataset, select it.
}

# Select the dataset based on the index obtained.
gset <- gset[[idx]]

# Extract the gene expression values from the selected dataset.
ex <- exprs(gset)

# Perform a log2 transformation under certain conditions.
# Calculate quantiles of the gene expression data.
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))

# Check conditions for log2 transformation.
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)

# If LogC is TRUE, apply the log2 transformation to the data.
if (LogC) {
  # Replace values less than or equal to 0 with NaN.
  ex[which(ex <= 0)] <- NaN
  # Perform the log2 transformation on the remaining values.
  ex <- log2(ex)
}

# box-and-whisker plot
par(mar = c(7, 4, 2, 1))
title <- paste ("GSE2037", "/", annotation(gset), sep ="")
boxplot(ex, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste ("GSE2037", "/", annotation(gset), " value distribution", sep = "")
plotDensities(ex, main=title, legend=F)

# mean-variance trend
ex <- na.omit(ex) # eliminate rows with NAs
plotSA(lmFit(ex), main="Mean variance trend, GSE2037")

# platofrm GPL199 - [Ecoli_ASv2] Affymetrix E. coli Antisense Genome Array
GPL199 <- read.delim(paste0(wd, "/input/gene_expression/GPL199.annot"))
rownames(GPL199) <- GPL199[, 1]
GPL199[, 1] <- NULL

ex = merge(ex, GPL199, by=0, all.x = TRUE)

# The identifier “bxxxx” is a gene identifier used in the KEGG database
# model link database -> "http://bigg.ucsd.edu/models/iML1515"

## Get a list of model reactions
# curl 'http://bigg.ucsd.edu/api/v2/models/iML1515/genes'

GRrules <- readRDS(paste0(wd, "/input/models/", model_type, "/", model_name,"/GRrules.rds"))
GenesGRrules <- readRDS(paste0(wd, "/input/models/", model_type, "/", model_name,"/genesfromGRrules.rds"))

# Print how many gene are envisaged
cat("num gene in expression dataset ...'", length(which(ex$Platform_ORF %in% GenesGRrules)), "'\n")

ex4model = ex[!is.na(match(ex$Platform_ORF, GenesGRrules)), ]
ex4model = ex4model[c(27, 2:15)]

# https://www.theseed.org/servers/TUTORIAL/SAPtutorial.html
# https://www.theseed.org/servers/
write.csv(ex4model, paste0(
  wd, "/input/gene_expression/GSE2037", "_", model_name, ".csv"), 
  row.names=FALSE)
