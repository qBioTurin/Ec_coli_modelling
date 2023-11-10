
# setting the proper working directory
wd = getwd()
supp_function.dir = paste0(wd, "/code/supplementary_functions/")

source(paste0(wd, "/code/scripts/libraries.R"))
source(paste0(wd, "/code/scripts/plotting_functions.R"))
source(paste0(wd, "/code/scripts/generating_carbon regimes.R"))
source(paste0(wd, "/code/scripts/generating_constraints_MOPS_media.R"))

source(paste0(supp_function.dir, "FBAmodel.generation.R"))
# loading experiments executioner
source(paste0(supp_function.dir, "Exe.R"))

# setting working directory
setwd(wd)

# FBA model compiling phase
# setting FBA model tags
model_type = "ecoli"
model_name = "iML1515"

model_wd = paste0(wd, "/input/models/ecoli/", model_name)
model_R_file = paste0(model_wd, "/", model_name, ".RData")
GeneExpress.file = paste0(wd, "/input/gene_expression/GSE2037", "_", model_name, ".csv")
net_wd = paste0(wd, "/net")
net_fname = paste0(model_name, "_PN")

change_manual = c("BIOMASS_Ec_iML1515_core_75p37M", "FHL", MOPS_mediadb$exchanges)
lb_manual = c(-1e-12, -1e-12, -MOPS_mediadb$amount_batch)
ub_manual = c(biomass_flux, 1e-12, rep(1e-12, length(MOPS_mediadb$exchanges)))

# Measure execution time
time_taken_comp = system.time({
  # compiling the FBA models to R file format
  FBAmodel.gen(model_type = model_type,
               model_name = model_name,
               # you can choose if to model the environment using 
               # VMH or maintain template bounds
               # NutritionModeling = "VMH",
               NutritionModeling = "Template",
               diet_name = "EU_average",
               # (mmol/gWD*h)
               lb_off = NA, 
               # lb_off = -1e-08,
               write = T, 
               saveRData = T,
               # change_manual could a single reaction or a vector of reactions:
               # The genes encoding FHL are known to be active under anaerobic conditions, 
               # but this reaction is constrained to zero to avoid 
               # unrealistic aerobic hydrogen production.
               change_manual = change_manual,
               # change_manual and lb_manual must have the same length
               # (mmol/gDW*h)
               lb_manual = lb_manual,
               # (mmol/gWD*h)
               ub_manual = ub_manual,
               # if you have gene expression data to integrate
               gene.exp = T, 
               GeneExpress_file = GeneExpress.file,
               # glucose replicate 1
               gene.exp_sample = 1,
               wd = wd)
})

print(time_taken_comp)

# GreatMod unified multi-layer model generation & analysis
load(model_R_file)

all_react = readRDS(paste0(model_wd, "/param_all_react.rds"))
all_react = all_react[[1]]

system("rm -f dockerID *error.log *.log ExitStatusFile")

##### SETTING CARBON REGIMES #####
#
# 3 examples of glucose feeding time regimes or supplementation for a fed-batch:
# 
# 1) Constant_feeding:
#   The feed rate of glucose as growth-limiting substrate is constant throughout the culture.
# 
# 2) Linear_feeding:
#   The feed rate of glucose is increased linearly over time to match the growth rate of the cells
# 
# 3) Pulsed_feeding:
#   The glucose is added in repeated cycles of short duration, followed by a period of no feeding. 
#   e.g the same total amount of glucose was fed in repeated 300s (5 min) cycles, with the feed pump on 
#   for either 60 or 150 s during each cycle. 

carbon_reg = c("Constant_feeding", "Linear_feeding", 
                 "Pulsed_feeding_60", "Pulsed_feeding_150",
                 "blank")

# scale_color_manual(values = as.character(MoMAColors::moma.colors("Fritsch", 4)))
values = c("#0F8D7B", "#EADD17", "#532E6B", "#514A19", "darkgrey")

time_taken_analysis = system.time({
  # model analysis
  for (s in carbon_reg) {
    Exe.exp(model_cat = model_name,
            model_name = model_name,
            fba_fname = paste0(model_name, ".txt"),
            atol = 1e-06, 
            rtol = 1e-06,
            time.step = 1e-03, # (h) 
            f_time = 10, # (h)
            distance_measure = "ReferenceM",
            reference_data = paste0(wd, "/input/csv/ReferenceRanking.csv"),
            Exper = "Model_Analysis",
            carbon = s,
            debug = F,
            react = "EX_biomass_e",
            replicates = 1,
            param_target = "Tlcts",
            net_fname = net_fname,
            cores = detectCores(),
            wd = wd,
            supp_function.dir = supp_function.dir)
  }
})

print(time_taken_analysis)

# plotting results
# or pre-allocate for slightly more efficiency
exptrace = vector("list", length = length(carbon_reg))
expflux = vector("list", length = length(carbon_reg))

for (s in 1:length(carbon_reg)) {
  
  TracesPath = paste0(wd, "/results/", carbon_reg[s], 
                      "_Model_Analysis/", net_fname, "-analysis-1.trace")
  
  FluxesPath = paste0(wd, "/results/", carbon_reg[s], 
                      "_Model_Analysis/", net_fname, "-analysis-1-0.flux")
  
  trace = read.table(TracesPath, header = T)
  flux = read.table(FluxesPath, header = T)
  
  subtrace = trace %>% dplyr::select(Time, glc_e, lcts_e) %>% 
    tidyr::gather(key = "Places", value = "Marking", -Time)
  
  subflux = flux %>% 
    dplyr::select(Time, c("BIOMASS_Ec_iML1515_core_75p37M", "EX_glc_D_e", "EX_lcts_e", "LACZ")) %>% 
    tidyr::gather(key = "Reaction", value = "Flux", -Time)
  
  subtrace = cbind(subtrace, Scenario = rep(carbon_reg[s], length(subtrace$Marking)))
  subflux = cbind(subflux, Scenario = rep(carbon_reg[s], length(subflux$Flux)))
  
  exptrace[[s]] = subtrace
  expflux[[s]] = subflux
  
}

subtrace = do.call(rbind, exptrace)
subflux = do.call(rbind, expflux)

target_reacts = c("EX_glc_D_e", "EX_lcts_e", "LACZ", "BIOMASS_Ec_iML1515_core_75p37M")

# Filter the data based on target reactions
fluxes_filtered <- dplyr::filter(subflux, grepl(paste(target_reacts, collapse = "|"), Reaction))

subtrace_glc <- dplyr::filter(subtrace, grepl("glc_e", Places))
subtrace_lcts <- dplyr::filter(subtrace, grepl("lcts_e", Places))

# Create ggplots
# FBA solutions
F_glc <- create_flux_plot(data = fluxes_filtered, reaction_name = "EX_glc_D_e")
F_lcts <- create_flux_plot(fluxes_filtered, "EX_lcts_e")
f_lacz <- create_flux_plot(fluxes_filtered, "LACZ")
# A biomass flux of 0.001 (mmol/gDW*h) could be as a growth/no-growth cutoff.
F_bio <- create_flux_plot(data = fluxes_filtered, reaction_name = "BIOMASS_Ec_iML1515_core_75p37M")
# PN marking
G <- create_marking_plot(data = subtrace_glc, compound_name = "D-Glucose")
L <- create_marking_plot(data = subtrace_lcts, compound_name = "Lactose")

# Create the plots for carbon regimes
plot_cf <- plot_regimes(cf, title = carbon_reg[1], color  = values[1])
plot_lf <- plot_regimes(lf, title = carbon_reg[2], color  = values[2])
plot_pf_60 <- plot_regimes(pf_60[1:250, ], title = carbon_reg[3], color  = values[3])
plot_pf_150 <- plot_regimes(pf_150[1:250, ], title = carbon_reg[4], color = values[4])
plot_blank <- plot_regimes(blank, title = carbon_reg[5], color  = values[5])

# Define file paths
file_pic_1_pdf <- paste0(wd, "/results/plots/fluxes.pdf")
file_pic_1_png <- paste0(wd, "/results/plots/fluxes.png")
file_pic_2_pdf <- paste0(wd, "/results/plots/marking.pdf")
file_pic_2_png <- paste0(wd, "/results/plots/marking.png")

carbon_regimes_pdf <- paste0(wd, "/results/plots/carbon_regimes.pdf")
carbon_regimes_png <- paste0(wd, "/results/plots/carbon_regimes.png")

# Save ggplots to files
ggsave(file_pic_1_pdf, (F_glc | F_lcts | F_bio), width = 10, height = 3)
ggsave(file_pic_1_png, (F_glc | F_lcts | F_bio), width = 10, height = 3)
ggsave(file_pic_2_pdf, (G + L), width = 9, height = 3)
ggsave(file_pic_2_png, (G + L), width = 9, height = 3)

ggsave(carbon_regimes_pdf, (plot_cf / plot_lf / plot_pf_60 / plot_pf_150 / plot_blank), 
       width = 2.5, height = 8)
ggsave(carbon_regimes_png, (plot_cf / plot_lf / plot_pf_60 / plot_pf_150 / plot_blank), 
       width = 2.5, height = 8)
