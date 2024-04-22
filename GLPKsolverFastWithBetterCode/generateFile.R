wd = getwd()
print(wd)


source(paste0(wd, "/epimod_FBAfunction_fast/FBAgreatmodeClass.R"))
source(paste0(wd, "/epimod_FBAfunction_fast/class_generation.R"))
source(paste0(wd, "/epimod_FBAfunction_fast/readMat.R"))



model_type = "ecoli"
model_name = "iML1515"
fba_fname = paste0(model_name, ".txt")

matfile = paste0(wd, '/input/', model_name, ".mat")

model = FBA4Greatmod.generation(fba_mat = matfile);

# Define output path for .rds files and move generated .rds files
output_rds_path <- paste0(wd, "/input/models/ecoli")
file_pattern <- "\\.rds$"  # Identify .rds files
files_to_move <- list.files(pattern = file_pattern)  # Find .rds files in the current directory
if (length(files_to_move) > 0) {
  sapply(files_to_move, function(f) {
    file.rename(from = f, to = file.path(output_rds_path, f))
  })
}


writeFBAfile(model, paste0('./data/', model_name, '.txt'))

# ESEGUI DA TERMINALE IL FILE C++ PER I RISULTATI