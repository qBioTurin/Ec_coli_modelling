library(dplyr)
library(ggplot2)
library(cowplot)

ModelAnalysisPlot=function(Traces_path, Flx_path) {
  
  par(mfrow=c(2,1))
  
  trace = read.delim(Traces_path, header = TRUE, sep = "")
  
  data = data.frame(var1 = trace$glc,
                    var2 = trace$lcts,
                    Time = trace$Time)
  
  pCS.glc = ggplot(data, aes(date)) + 
    geom_line(data=data,
              aes(x=Time, y=var1), 
              color = "blue", size = 1)+
    theme(legend.position = c(0.95, 0.95),
          legend.justification = c("right", "top")) +
    xlab("Time (h)") +
    ylab("glc concentration (mmol/L)") +
    ylim(0, 6)
  
  pCS.lcts = ggplot(data, aes(date)) +
    geom_line(data=data,
              aes(x=Time, y=var2), 
              color = "red", size = 1) +
    theme(legend.position = c(0.95, 0.95),
          legend.justification = c("right", "top")) +
    xlab("Time (h)") +
    ylab("lcts concentration (mmol/L)")
  
  pCS = plot_grid(pCS.glc, pCS.lcts, labels=c("A", "C"), ncol = 1, nrow = 2)
  
  fluxes = read.delim(Flx_path, header = F, sep = ";")
  
  data.flx = data.frame(var1 = fluxes$V2,
                        var2 = fluxes$V3,
                        var3 = fluxes$V4,
                        Time = fluxes$V1)
  
  data.flx = unique(data.flx)
  
  pFlx.glc = ggplot(data.flx, aes(date)) + 
    geom_line(data = data.flx,
              aes(x=Time, y=var1), 
              color = "blue", size = 1) +
    theme(legend.position = c(0.95, 0.95),
          legend.justification = c("right", "top")) +
    xlab("Time [FBA call]") +
    ylab("flux EX_glc_D(e) [(mmol/gDW)/h]")+
    ylim(-3, 3)
  
  pFlx.lcts = ggplot(data.flx, aes(date))+
    geom_line(data = data.flx,
              aes(x=Time, y=var2), 
              color = "red", size = 1) +
    theme(legend.position = c(0.95, 0.95),
          legend.justification = c("right", "top")) +
    xlab("Time [FBA call]") +
    ylab("flux EX_lcts(e) [(mmol/gDW)/h]")+
    ylim(-3, 3)
  
  pFlx = plot_grid(pFlx.glc, pFlx.lcts, labels=c("B", "D"), ncol = 1, nrow = 2)
  
  pComp = plot_grid(pCS, pFlx, labels=c("", ""), ncol = 2, nrow = 1)
  
  new_var = c()
  
  # One unit of OD600 corresponds to a cell dry weight of. Value, 0.3 g/L
  e.coli_gDW_OD600 = 0.3 # [g]
  
  new_var[1] = e.coli_gDW_OD600 # iniCellDensity [g/l]
  # MaxsDensity = 5 [g/l]
  
  for (i in 1:length(data.flx$var3)) {
    new_var[i+1] =  new_var[i] * exp((data.flx$var3[i+1]*0.1)*((5-new_var[i])/5))
  }
  
  data.flx$var3  = head(new_var, -1)
    
  pBiomass = ggplot(data.flx, aes(date))+
    geom_line(data = data.flx,
              aes(x=Time, y=var3),  
              color = "green", size = 1) +
    theme(legend.position = c(0.95, 0.95),
          legend.justification = c("right", "top")) +
    xlab("Time (h)") +
    ylab("Cell Density (g/l)") +
    ylim(0, 6)
  
  plot_grid(pComp)
  # plot_grid(pComp, pBiomass, labels=c("", "E"), ncol = 2, nrow = 1)
  
}
