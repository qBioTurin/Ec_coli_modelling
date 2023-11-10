
# Define a function to create a ggplot
create_flux_plot <- function(data, reaction_name) {
  
  data_colors <- data.frame(
    Scenario = carbon_reg,
    Color = values
  )

  ggplot(dplyr::filter(data, grepl(reaction_name, Reaction)), 
         aes(Time, Flux, colour = factor(Scenario, levels = data_colors$Scenario))) +
    geom_line(size = 1.25) +
    ggtitle(ifelse(
      reaction_name == "BIOMASS_Ec_iML1515_core_75p37M", "BIOMASS_iML1515", 
      reaction_name)) +
    xlab("Time (hour)") +
    ylab(expression("Reaction Fluxes (mmol/gDW*h)")) +
    scale_color_manual(values = data_colors$Color) +
        theme(plot.title = element_text(color = "black", size = 13, face = "bold"),
              axis.title = element_text(size = 10),
              axis.text = element_text(size = 10),
              legend.position = "none",
              panel.background = element_rect(
                fill = "#ededed", colour = "#ededed", 
                linewidth = 0.5, linetype = "solid"),
              plot.background = element_rect(fill = "white"))
}

# Define a function to create a ggplot for extracellular compounds
create_marking_plot <- function(data, compound_name) {
  
  data_colors <- data.frame(
    Scenario = carbon_reg,
    Color = values
  )
  
  ggplot(data, aes(
    Time, Marking, colour = factor(Scenario, levels = data_colors$Scenario))) +
    geom_line(size = 1) +
    ggtitle(compound_name) +
    xlab("Time (hour)") +
    ylab(paste("Extracellular", compound_name, "(mmol)")) +
    labs(color = "carbon regimes") +
    scale_color_manual(values = data_colors$Color) +
        theme(plot.title = element_text(color = "black", size = 17, face = "bold"),
              axis.title = element_text(size = 12),
              axis.text = element_text(size = 10),
              legend.position = ifelse(compound_name == "D-Glucose", "none", "right"),
              legend.direction = "vertical",
              legend.box = "vertical",
              legend.title = element_text(size = 12, face = "bold"),
              legend.text = element_text(size = 10),
              panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "grey"))
}

# Define the function
plot_regimes <- function(df, title, color) {
  p <- ggplot(df, aes(x=time_h)) +
    geom_line(aes(y=glc, linetype="glc"), color=color) +
    geom_line(aes(y=lcts, linetype="lcts"), color=color) +
    labs(title=title, x="Time (h)", y = "Value (mmol)") +
    theme_minimal() +
    theme(
      plot.title = element_text(color = "black", size = 9, face = "bold")
    )
  
  # Add legend only for the first plot
  if (title == "cf") {
    p <- p + labs(linetype="Variable")
  } else {
    p <- p + theme(legend.position="none")
  }
  
  return(p)
}
