library(parallel)
library(dplyr)
library(ggplot2)
library(patchwork)
library(epimod)
library(grid)
library(gtable)
library(stringr)
library(ggplotify)

changecolorsFacet = function(pt,colors){
  # Convert the plot to a gtable object
  g <- ggplotGrob(pt)
  
  # Locate the strip grobs for the x-axis (top strips)
  strip_idx <- c(which(grepl("strip-t", g$layout$name)),which(grepl("strip-r", g$layout$name)))
  
  # Customize the strip backgrounds
  for (i in seq_along(strip_idx)) {

    strip_grob <- g$grobs[[strip_idx[i]]]
    label_grob <- strip_grob$grobs[[1]]$children[[which(grepl("text", names(strip_grob$grobs[[1]]$children)))]]
    facet_label <- label_grob$children[[1]]$label
    
    # Set the fill color based on the mapping
    fill_color <- colors[facet_label]
    
    # Replace the background rect
    g$grobs[[strip_idx[i]]]$grobs[[1]]$children[[1]]$gp$fill <- fill_color
  }
  
  return(as.ggplot(g))
  
}

ModelAnalysisPlot=function(TracesPath, FluxPath, FluxVec,new_eps_value,tag) {
  
  trace = read.table(TracesPath, header = T)
  
  subtrace = trace %>% tidyr::gather(key = "Places", value = "Marking", -Time)
  
  
  flux = read.csv(FluxPath, header=FALSE, sep=" ")
  colnames(flux) = flux[1, ]
  flux = flux[-1, ]
  
  if(is.null(FluxVec)){
    
    FluxVec = colnames(flux)
    flux %>% dplyr::select(FluxVec, Time) %>% tidyr::gather(key = "Reaction", value = "Flux", -Time)
    
    subflux = flux %>% select(FluxVec, Time) %>%
      tidyr::gather(key = "Reaction", value = "Flux", -Time) %>%
      ggplot2::ggplot() + ggplot2::geom_line(ggplot2::aes(x=Time, y=Flux, group = Reaction)) +
      ggplot2::facet_wrap(~Reaction, scales = "free_y") + ggplot2::scale_fill_brewer("Accent")
    
    print(subflux, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
    
  } else {
    
    subflux = flux[,c(FluxVec, "Time")] %>%
      tidyr::gather(key = "Reaction", value = "Flux", -Time)
    
    subflux$Time = as.numeric(as.character(subflux$Time))
    subflux$Flux = as.numeric(as.character(subflux$Flux))
  }
  
  subflux$new_eps_value = subtrace$new_eps_value = new_eps_value
  subflux$tag = subtrace$tag = tag
  
  return(list(subflux = subflux,subtrace = subtrace))
}

set_closest_directory <- function(target_name) {
  wd = getwd()
  if( basename(wd) != target_name){
    # Get all directories relative to the current working directory
    all_dirs <- list.dirs(path = ".", full.names = TRUE, recursive = TRUE)
    
    # Filter directories with the specific name
    matching_dirs <- all_dirs[basename(all_dirs) == target_name]
    
    # Check if there are matches
    if (length(matching_dirs) == 0) {
      stop("No folder with the name '", target_name, "' found in the current working directory or its subdirectories.")
    }
    
    # Find the directory with the shortest path (closest to the current working directory)
    closest_dir <- matching_dirs[which.min(nchar(matching_dirs))]
    
    # Set the working directory
    setwd(closest_dir)
  }
  # Return the new working directory
  return(getwd())
}

wd <- set_closest_directory("Ec_coli_modelling")
setwd(wd)

system("rm -f dockerID *error.log *.log ExitStatusFile" )

if(!dir.exists("./results")) dir.create("./results")
if(!dir.exists("./results/TimingEval")) dir.create("./results/TimingEval")

######## Variables definition:
model_type = "ecoli"
model_name = "iML1515"
carbon_reg = c("Constant_feeding", "Linear_feeding", 
               "Pulsed_feeding_60", "Pulsed_feeding_150",
               "blank")
fba_fname = paste0(model_name, ".txt")
distance_measure = "ReferenceM"
reference_data = paste0(wd,"/input/csv/ReferenceRanking.csv")
Exper = "Model_Analysis"
react = "EX_biomass_e"
param_target = "Tlcts"
net_fname =  paste0(model_name, "_PN")
cores = detectCores()
supp_function.dir = "/code/supplementary_functions/"

model.generation(net_fname = paste0(wd, "/net/", net_fname, ".PNPRO"),
                 transitions_fname = paste0(wd, "/net/transitions.cpp"),
                 fba_fname = paste0(wd, "/input/compiled_models/", fba_fname))
system(paste0("mv ", net_fname, ".* ./net"))

### Color association ####

epstimes = c("1e-6", "1e-5",  "1e-4","1e-3", "1e-2")
colors_new_eps_value <- c("#30123BFF", "#BB1656FF", "#5801A4FF" ,"#EDD03AFF" , "#D23105FF")
names(colors_new_eps_value) <- epstimes

place = c("glc_e","lcts_e")
colors_place <- c("#862781FF", "#FEB67CFF")
names(colors_place) <- place

scenarios = c("Constant_feeding","Linear_feeding","Pulsed_feeding_60","Pulsed_feeding_150", "blank"  )
colors_scenarios <- c("#0B0405FF" , "#3B5698FF", "#359BAAFF", "#49C1ADFF", "#96DDB5FF")
names(colors_scenarios) <- scenarios

####
paramsgrid = expand.grid(epstimes, carbon_reg)

MultipleAnalysis = lapply(seq_along(paramsgrid[, 1]),
                          function(c, paramsConfig, paramsgrid) {
                            # browser()
                            new_eps_value = paramsgrid[c, "Var1"]
                            tag = paramsgrid[c, "Var2"]
                            
                            # Read the CSV file into a data frame
                            csv = read.csv(paste0(wd, "/input/csv/CarbonAdmin_eps.csv"), header = F, quote = "")
                            # Make changes to the data frame
                            # csv[1, ] = paste0("i; init; init.gen; iG = ", iG, "; iL = 0;")
                            csv[2, ] = paste0("g; M; MatrixGeneration; frame='/home/docker/data/input/csv/carbon_regimes/", tag, ".csv';")
                            csv[3, ] = paste0("g; eps; ", new_eps_value)
                            write.table(csv, paste0(wd, "/input/csv/CarbonAdmin_eps.csv"), col.names = F, row.names = F, quote = F)

                            execution_start <- Sys.time()
                            model.analysis(solver_fname = paste0(wd, "/net/", net_fname, ".solver"),
                                           i_time = 0, 
                                           f_time = 10, 
                                           s_time = 0.5,
                                           FVA = F,debug = T,
                                           fba_fname = paste0(wd, "/input/compiled_models/", fba_fname),
                                           parameters_fname = paste0(wd, "/input/csv/CarbonAdmin_eps.csv"),
                                           functions_fname = paste0(wd, "/",supp_function.dir, "functions.R"))
                            execution_end <- Sys.time()
                            
                            execution_time = as.numeric(difftime(execution_end, execution_start, units = "secs"))
                            
                            lines <- readLines(list.files(path = wd, pattern = "\\.log$", full.names = TRUE))
                            file.remove(list.files(path = wd, pattern = "\\.log$", full.names = TRUE))
                            
                            extracted_lines <- grep("Total memory used|Total time required", lines, value = TRUE)
                            
                            n_FBA = length(grep("FBA call", lines, value = TRUE))
                            
                            times <- as.numeric(sub(".*Total time required: ([0-9]+)s.*", "\\1", extracted_lines[grepl("time required", extracted_lines)]))
                            memory <- as.numeric(sub(".*Total memory used: ([0-9]+)KB.*", "\\1", extracted_lines[grepl("memory used", extracted_lines)]))
                            
                            result_df <- data.frame(Time_s = sum(times),
                                                    Memory_KB = sum(memory),
                                                    GlobalExecution_time = execution_time,
                                                    eps = new_eps_value,
                                                    FBA_calls = n_FBA,
                                                    Scenario = tag)
                            
                            resFolder = paste0("./results/TimingEval/Ec_coli", "_", tag, "_", new_eps_value)
                            
                            if(dir.exists(resFolder)) system(paste("rm -r ", resFolder))
                            
                            system(paste0("mv ", net_fname, "_analysis* ", resFolder))
                            
                            traces = ModelAnalysisPlot(
                              paste0(resFolder,"/iML1515_PN-analysis-1.trace"),
                              paste0(resFolder,"/iML1515_PN-analysis-1-0.flux"), 
                              FluxVec = c("BIOMASS_Ec_iML1515_core_75p37M", "EX_glc_D_e", "EX_lcts_e", "LACZ"),
                              new_eps_value,
                              tag
                            )
                            
                            debug_plot = ggplot(
                              traces[["subtrace"]]) +
                              geom_line(aes(x = Time, y = Marking,
                                            linetype = tag),
                                        linewidth = 1) +
                              facet_wrap(~Places, scales = "free_y") +
                              theme_minimal() +
                              theme(legend.position = "bottom")
                            
                            # Save the debug plot
                            ggsave(filename = paste0(wd, "/results/plots/", new_eps_value,
                                                    "_", tag, "_debug_plot.pdf"), 
                                   plot = debug_plot,
                                   width = 12, height = 8)
                            
                            return(list(traces = traces, result_df = result_df))
                          },
                          paramsConfig = resParams$Config, paramsgrid
)

saveRDS(MultipleAnalysis, file = "trajectories.rds")
MultipleAnalysis = readRDS(file = "trajectories.rds")

timing = do.call(rbind, lapply(MultipleAnalysis, "[[", 2))
flux = do.call(rbind, lapply(lapply(MultipleAnalysis,"[[", 1),"[[", 1))
trajectories = do.call(rbind, lapply(lapply(MultipleAnalysis, "[[",1), "[[", 2))
trajectories$new_eps_value = factor(trajectories$new_eps_value, levels = epstimes)
flux$new_eps_value = factor(flux$new_eps_value, levels = epstimes)

###

plot_varying_eps = ggplot(trajectories ) +
  geom_line(aes(x = Time, y = (Marking),  col = new_eps_value)) +
  scale_color_manual(values = colors_new_eps_value[unique(trajectories$new_eps_value)])+
  facet_grid(Places~tag,scales = "free") +
  theme_bw() +
  theme(
    plot.subtitle = element_text(size = 10, face = "bold", color = "#2a475e"),
    plot.title.position = "plot", 
    axis.text = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 15, face = "bold"),
    legend.key.size = unit(0.4, "cm"),
    strip.text.x = element_text(size = 10, face = "bold", colour = "white"),
    strip.text.y = element_text(size = 10, face = "bold", colour = "white"),
    legend.position = "top") +
  labs(x = "Time (h)", y = "Quantity", col = "Threshold") +
  geom_vline(xintercept = 8)

plot_varying_eps

# Calculate percentage difference

df_diff <- trajectories %>%
  mutate(Marking = round(Marking,digits = 12))%>%
  mutate(new_eps_value = as.character(new_eps_value) ) %>% 
  mutate(new_eps_value = str_replace(new_eps_value, "1e-6", "baseline")) %>%
  group_by(Time, tag) %>%
  tidyr::spread(key = new_eps_value, value = Marking) %>%
  tidyr::gather(-Time,-Places,-baseline, -tag, key = "new_eps_value", value = "Marking") %>%
  mutate(perc_diff = if_else(baseline!= 0 , 100 * (Marking - baseline) / baseline , 100 * Marking))

pt_perc = ggplot(df_diff)+
  geom_bar(aes(x = Time, y = perc_diff, group = tag, fill = Places),
           stat = "identity", position = "dodge")+
  scale_fill_manual(values = colors_place[unique(trajectories$Places)])+
  facet_grid(tag~new_eps_value,scales = "free")+
  theme_bw()+
  theme(
    plot.subtitle = element_text(size = 10, face = "bold", color = "#2a475e"),
    plot.title.position = "plot", 
    axis.text = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 15, face = "bold"),
    legend.key.size = unit(0.4, "cm"),
    strip.text.x = element_text(size = 10, face = "bold", colour = "white"),
    strip.text.y = element_text(size = 10, face = "bold", colour = "white"),
    legend.position = "none")+
  labs(x = "Time (h)",y = "% difference with 1e-6")
pt_perc = changecolorsFacet(pt_perc,c(colors_new_eps_value,colors_scenarios) )

pt = ggplot(trajectories) +
  geom_line(aes(x = Time, y = Marking,  col = new_eps_value ))+
  scale_color_manual(values = colors_new_eps_value[unique(trajectories$new_eps_value)])+
  facet_grid(Places~tag,scales = "free")+
  theme_bw()+
  theme(
    plot.subtitle = element_text(size = 10, face = "bold", color = "#2a475e"),
    plot.title.position = "plot", 
    axis.text = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 15, face = "bold"),
    legend.key.size = unit(0.4, "cm"),
    strip.text.x = element_text(size = 10, face = "bold", colour = "white"),
    strip.text.y = element_text(size = 10, face = "bold", colour = "white"),
    legend.position = "top")+
  labs(x = "Time (h)",y = "Quantity", col = "Threshold")
pt = changecolorsFacet(pt,c(colors_place,colors_scenarios) )

Fig2D = (pt|pt_perc) & theme(legend.position = "bottom")
Fig2D

## ggsave(plot = Fig2D, filename = "Figures/Fig2D.pdf", width = 16, height = 18)

#################

range_time <- range(timing$GlobalExecution_time)
range_calls <- range(timing$FBA_calls)
scale_factor <- diff(range_time) / diff(range_calls)
timing$FBA_calls_scaled <- (timing$FBA_calls - min(range_calls)) * scale_factor + min(range_time)

timing$ID = 1
timing$eps = factor(timing$eps,levels = epstimes)

p <- ggplot(timing, aes(x = eps, group = ID) ) +
  geom_line(aes(y = GlobalExecution_time, color = "Global Time (s)"), linewidth = 1, linetype = "dashed") +  # Time_s line
  geom_line(aes(y = FBA_calls_scaled, color = "Number of FBA"), linewidth = 1) +  # Scaled Memory_KB line
  scale_y_continuous(
    name = "Time (s)",
    sec.axis = sec_axis(~ (. - min(range_time)) / scale_factor + min(range_calls), 
                        name = "Number of FBA")
  ) +
  scale_color_manual(values = c("Number of FBA" = "red","Global Time (s)" = "blue")) +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 10, face = "bold", color = "#2a475e"),
    plot.title.position = "plot", 
    axis.text = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 15, face = "bold"),
    legend.key.size = unit(0.4, "cm"),
    strip.text.x = element_text(size = 10, face = "bold", colour = "white"),
    strip.text.y = element_text(size = 10, face = "bold", colour = "black"),
    strip.background.y = element_rect( fill = "white"),
    legend.position = "nonw",
    axis.title.y.left = element_text(color = "blue"),  # Secondary axis styling
    axis.text.y.left = element_text(color = "blue"),
    axis.title.y.right = element_text(color = "red"),  # Secondary axis styling
    axis.text.y.right = element_text(color = "red"),
    legend.title = element_blank()
  ) +
  labs(x = "Threshold (eps)", color = "")+
  facet_grid(~Scenario)


pl = (Fig2D + plot_layout(guides = "collect")) /p + plot_layout(heights = c(1,0.4)) + plot_annotation(tag_levels = 'A')

pl

ggsave(plot = pl, filename = "results/plots/PerformanceEcoli.pdf", width = 21, height = 12)

saveRDS(pl,file = "results/plots/PerformanceEcoli.RDs")
