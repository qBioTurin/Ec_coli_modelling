# Load libraries
library(epimod)
library(epimodFBAfunctions)
library(tools)

wd = getwd()
print(wd)

# Set the values for the E. coli model
model_type = "ecoli"
model_name = "iML1515"
fba_fname = paste0(model_name, ".txt")

# Paths for network files
net_wd = paste0(wd, "/net")
net_fname = paste0(model_name, "_PN")
print(net_fname)

# Path for the matrix file for FBA
matfile = paste0(wd, '/input/', model_name, ".mat")

# Generate the FBA model
model = FBA4Greatmod.generation(fba_mat = matfile)

# Define output path for .rds files and move generated .rds files
output_rds_path <- paste0(wd, "/input/models/ecoli")
file_pattern <- "\\.rds$"  # Identify .rds files
files_to_move <- list.files(pattern = file_pattern)  # Find .rds files in the current directory
if (length(files_to_move) > 0) {
  sapply(files_to_move, function(f) {
    file.rename(from = f, to = file.path(output_rds_path, f))
  })
}

# Write the FBA file to the specified path
writeFBAfile(model, paste0('./result/', model_name, '.txt'))

# Model Generation
model.generation(
  net_fname = paste0(wd, "/net/", net_fname, ".PNPRO"),
  transitions_fname = paste0(wd, "/net/transitions.cpp"),
  fba_fname = paste0(wd, "/result/", fba_fname),
  debug = TRUE
)

# Move generated network files to the /net folder
system(paste0("mv ", wd, "/", net_fname, ".* ./net"))

# Model analysis 
model.analysis(
  parameters_fname = paste0(wd, "/OdeData/", 'CarbonAdmin.csv'),
  functions_fname = paste0(wd, '/OdeData/', "functions.R"),
  solver_fname = paste0(wd, "/net/", net_fname, ".solver"), 
  i_time = 0, f_time = 10, 
  s_time = 0.5, 
  fba_fname = paste0(wd, "/result/", fba_fname), 
  event_function = NULL,
  debug = TRUE
)
