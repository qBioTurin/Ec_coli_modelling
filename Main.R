
library(epimod)
library(dplyr)
library(ggplot2)
library(patchwork)

setwd("~/Documents/ClostridiumDiff_FBAandPN/Ec_coli_modelling")

################################################################################
########## Let's compile the FBA model from the RData storing the CRN ##########
################################################################################

# source("./Rfunction/FBAfileGeneration.R")
# RDataPresent = F
# NutritionModeling = "VMH"
# ShowExBounds = T
# write = T

# model.name = "Ec_K12"
# path = "./Input/Models/Ec_K12/Ec_K12.mat"

# "EU_average"
# lbNOTdiet = -0.1
# uppbndBiomass = 5
# lowbndBiomass = lbNOTdiet
# Ub = 10

################################################################################
########################### epimod Model Generation ############################
################################################################################

model.generation(net_fname = "./Net/Ec_corePN.PNPRO",
                 transitions_fname = "./Net/transitions.cpp", 
                 fba_fname = "FBAmodelK12")

system("mv Ec_corePN.* ./Net")

################################################################################
############################ epimod Model Analysis #############################
################################################################################

model.analysis(solver_fname = "./Net/Ec_corePN.solver",
               i_time = 0, f_time = 42, s_time = 0.1, 
               fba_fname = "./Input/CompiledModels/FBAmodelK12", 
               atol = 1e-08, rtol = 1e-08, debug = T, 
               parameters_fname = "Input/csv/CarbonAdmin.csv",
               functions_fname = "Rfunction/Functions.R")

file.remove("dockerID")
file.remove("ExitStatusFile")

################################################################################
############################### Visualizations #################################
################################################################################

grDevices::pdf(file = "./Results/CarbonSourcesUtilization.pdf",
               width = 12, height = 8)

conditions = c("AlternedPulses", "PairedPulses", "StepwisePulses")

for (con in conditions) {
  
  TracesPath = paste("Results/Ec_corePN_", con, 
                     "/Ec_corePN-analysis-1.trace", sep = "")
  
  trace = read.table(TracesPath, header = F)
  
  ReactionsNames = trace[1, ]
  ReactionsNames = gsub("\\(", replacement = "_", ReactionsNames)
  ReactionsNames = gsub("\\)", replacement = "", ReactionsNames)
  
  trace = read.table(TracesPath, header = T)
  
  colnames(trace) = ReactionsNames
  
  subflux = trace %>% 
    dplyr::select(Time, c("EX_biomass_e", "biomass525", "EX_glc_D_e", 
                          "EX_lcts_e", "LACZ", "LCTSt")) %>% 
    tidyr::gather(key = "Reaction", value = "Flux", -Time)
  
  subtrace = trace %>% dplyr::select(Time, glc_e, lcts_e) %>% 
    tidyr::gather(key = "Places", value = "Marking", -Time)
  
  assign(paste("subtrace", con, sep = ""),
         cbind(subtrace, Scenario = rep(con, length(subtrace$Marking))))
  
  assign(paste("subflux", con, sep = ""), 
         cbind(subflux, Scenario = rep(con, length(subflux$Flux))))
  
}

subtrace = rbind(subtraceAlternedPulses, subtracePairedPulses, subtraceStepwisePulses)
subflux = rbind(subfluxAlternedPulses, subfluxPairedPulses, subfluxStepwisePulses)

fluxes = filter(subflux, grepl("EX_glc_D_e|EX_lcts_e|LACZ|LCTSt", Reaction))

for (i in 1:length(fluxes$Reaction)) {
  if(fluxes$Reaction[i] == "EX_glc_D_e") {
    fluxes$Reaction[i] = "EX_glc_D(e)"}
}

for (i in 1:length(fluxes$Reaction)) {
  if(fluxes$Reaction[i] == "EX_lcts_e") {
    fluxes$Reaction[i] = "EX_lcts(e)"}
}

f_lacz = ggplot(filter(fluxes, Reaction == "LACZ"),
                  aes(Time, Flux, colour = Scenario)) +
  geom_line(size = 0.80) + theme_gray(base_size = 14) +
  ggtitle("LACZ speed") +
  xlab("Time (hour)") + labs(color  = "", linetype = "") +
  scale_color_manual(values = c('darkred', 'darkcyan', "purple")) +
  theme(plot.title = element_text(color="black", size = 10, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none") +
  labs(y = expression("Reaction Fluxes (mmol)"))

f_LCTSt = ggplot(filter(fluxes, Reaction == "LCTSt"),
               aes(Time, Flux, colour = Scenario)) +
  geom_line(size = 0.80) + theme_gray(base_size = 14) +
  ggtitle("LCTSt speed") +
  xlab("Time (hour)") + labs(color  = "", linetype = "") +
  scale_color_manual(values = c('darkred', 'darkcyan', "purple")) +
  theme(plot.title = element_text(color="black", size = 10, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none") +
  labs(y = expression("Reaction Fluxes (mmol)"))

F_glc = ggplot(filter(fluxes, Reaction == "EX_glc_D(e)"),
                  aes(Time, Flux, colour = Scenario)) +
  geom_line(size = 0.80) + theme_gray(base_size = 14) +
  ggtitle("EX_glc_D(e)") +
  xlab("Time (hour)") + labs(color  = "", linetype = "") +
  scale_color_manual(values = c('darkred', 'darkcyan', "purple")) +
  theme(plot.title = element_text(color="black", size = 10, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none") +
  labs(y = expression("Reaction Fluxes (mmol)"))

F_lcts = ggplot(filter(fluxes, Reaction == "EX_lcts(e)"),
              aes(Time, Flux, colour = Scenario)) +
  geom_line(size = 0.80) + theme_gray(base_size = 14) +
  ggtitle("EX_lcts(e)") +
  xlab("Time (hour)") + labs(color  = "", linetype = "") +
  scale_color_manual(values = c('darkred', 'darkcyan', "purple")) +
  theme(plot.title = element_text(color="black", size = 10, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none") +
  labs(y = expression("Reaction Fluxes (mmol)"))

D_Glucose = filter(subtrace, grepl("glc_e", Places))
Lactose = filter(subtrace, grepl("lcts_e", Places))

G = ggplot(D_Glucose, aes(Time, Marking, colour = Scenario)) + 
  geom_line(size = 0.80) + theme_bw() +
  ggtitle("D-Glucose") +
  xlab("Time (hour)") + 
  ylab("extracellular Glucose (mmol)") +
  labs(colour = "Condition") +
  scale_color_manual(values = c('darkred', 'darkcyan', "purple")) +
  theme(plot.title = element_text(color="black", size = 10, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.14, 0.8),
        legend.background = 
          element_rect(size = 0.20, fill = "white", color = "black"))

L = ggplot(Lactose, aes(Time, Marking, colour = Scenario)) + 
  geom_line(size = 0.80) + theme_bw() +
  ggtitle("Lactose") +
  xlab("Time (hour)") + 
  ylab("extracellular Lactose (mmol)") +
  labs(colour = "Condition") +
  scale_color_manual(values = c('darkred', 'darkcyan', "purple")) +
  theme(plot.title = element_text(color="black", size = 10, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.14, 0.8),
        legend.background = 
          element_rect(size = 0.20, fill = "white", color = "black"))

(F_lcts | F_glc | f_LCTSt | f_lacz) / (G + L)
  
dev.off()