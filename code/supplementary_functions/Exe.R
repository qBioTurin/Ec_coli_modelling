
Exe.exp = function(model_cat, 
                   model_name,
                   fba_fname,
                   atol, 
                   rtol,
                   f_time,
                   distance_measure,
                   reference_data,
                   Exper, 
                   carbon,
                   debug,
                   react,
                   replicates, 
                   param_target,
                   time.step,
                   net_fname,
                   wd,
                   supp_function.dir,
                   cores) {
  
  setwd(wd)
  
  system("rm -f dockerID *error.log *.log ExitStatusFile")
  
  # Read the CSV file into a data frame
  csv = read.csv(paste0(wd, "/input/csv/CarbonAdmin.csv"), header = F, quote = "")
  # Make changes to the data frame
  csv[2, ] = paste0("g; M; MatrixGeneration; frame='/home/docker/data/input/csv/carbon_regimes/", s, ".csv';")
  
  csv[1, ] = paste0("i; init; init.gen; iG = ", iG, "; iL = 0;")
  
  write.table(csv, paste0(wd, "/input/csv/CarbonAdmin.csv"), col.names = F, row.names = F, quote = F)
  
  model.generation(net_fname = paste0(wd, "/net/", net_fname, ".PNPRO"),
                   transitions_fname = paste0(wd, "/net/transitions.cpp"),
                   fba_fname = paste0(wd, "/input/compiled_models/", fba_fname))
  
  system(paste0("mv ", wd, "/", net_fname, ".* ./net"))
  
  switch(Exper, 
         Model_Analysis = {
           parallel_processors = 1
           n_config = 1
           parameters_fname <<- "input/csv/CarbonAdmin.csv"
         },
         Model_Sensitivity = {
           parallel_processors = cores
           n_config = parallel_processors*replicates
           parameters_fname <<- "input/csv/CarbonAdmin.csv"
         },
         FBA_Sensitivity = {
           parallel_processors = cores
           n_config = parallel_processors*replicates
           parameters_fname <<- "input/csv/CarbonAdmin.csv"
         }
  )
  
  if( Exper == "Model_Analysis" ) {
    
    model.analysis(solver_fname = paste0(wd, "/net/", net_fname, ".solver"),
                   i_time = 0, 
                   f_time = f_time, 
                   s_time = time.step,
                   parallel_processors = parallel_processors,
                   n_config = n_config,
                   FVA = F,
                   fba_fname = paste0(wd, "/input/compiled_models/", fba_fname),
                   atol = atol, rtol = rtol, 
                   debug = debug,
                   parameters_fname = paste0(wd, "/", parameters_fname),
                   functions_fname = paste0(supp_function.dir, "functions.R"),
                   event_function = NULL)
    
    system(paste0("rm -r ./results/", s, "_", Exper))
    system(paste0("mv ", net_fname, "_analysis* ./results/", s, "_", Exper))
    
  }
  
  if(Exper == "Model_Sensitivity") {
    
    model.sensitivity(
      solver_fname = paste0(wd, "/net/", net_fname, ".solver"),
      i_time = 0, f_time = f_time, s_time = time.step,
      parallel_processors = parallel_processors,
      n_config = n_config,
      fba_fname = paste0(wd, "/input/compiled_models/", fba_fname),
      atol = atol, 
      rtol = rtol, 
      debug = debug,
      target_value = "glc_e",
      distance_measure = distance_measure, 
      reference_data = reference_data, 
      parameters_fname = paste0(wd, "/", parameters_fname),
      functions_fname = paste0(wd, supp_function.dir, "functions.R"),
      event_function = NULL)
    
  }
  
}