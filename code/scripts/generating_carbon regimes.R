
# the time resolution is approximately 5 [s] for step -> 0.001 [h]
t_0 = 0; t_f = 15;
points = 82
# 82 time points per interval
n_steps = points * floor(t_f*1000 / points)
time_window = seq(t_0, t_f, length.out = n_steps)

# generating blank
blank = data.frame(time_h = time_window, 
                   glc = rep(0, length(time_window)),
                   lcts = rep(0, length(time_window)), 
                   row.names = NULL, check.names = FALSE)

cf = data.frame(time_h = time_window, 
                glc = rep(0.5, length(time_window)),
                lcts = rep(0.2, length(time_window)), 
                row.names = NULL, check.names = FALSE)

lf = data.frame(time_h = time_window,
                glc = 0.5*time_window,
                lcts = 0.05*time_window, 
                row.names = NULL, check.names = FALSE)

# 300s ~= 0.08 hours -> 30s ~= 0.008 hours -> 60s ~= 0.016 hours -> 150s ~= 0.04 hours
pf_60 = data.frame(time_h = time_window,
                   glc = rep(c(rep(5, length(seq(0, 0.016, by = 0.001))), 
                               rep(0, length(seq(0.016, 0.08, by = 0.001)))), n_steps/points),
                   lcts = rep(c(rep(2, length(seq(0, 0.008, by = 0.001))), 
                                rep(0, length(seq(0.008, 0.08, by = 0.001)))), n_steps/points))

pf_150 = data.frame(time_h = time_window,
                    glc = rep(c(rep(5, length(seq(0, 0.04, by = 0.001))), 
                                rep(0, length(seq(0.04, 0.08, by = 0.001)))), n_steps/points),
                    lcts = rep(c(rep(2, length(seq(0, 0.02, by = 0.001))), 
                                 rep(0, length(seq(0.02, 0.08, by = 0.001)))), n_steps/points))

data = list(cf, lf, pf_60, pf_150, blank)

mat = matrix(0, nrow = dim(cf)[1], ncol = dim(cf)[2])
mock = paste0(wd, "/input/csv/carbon_regimes/mock.csv")
write.table(mat, file = mock, col.names = F, row.names = F, quote = F)

# Create a for loop to write each row of the data frame to a CSV file
for(df in 1:length(data)) {
  row.names(data[[df]]) = NULL
  colnames(data[[df]]) = NULL
  csv = readr::read_table(mock, col_names = FALSE)
  for (i in 1:nrow(data[[df]])) {
    csv[i, ] <- data[[df]][i, ]
  }
  filename = paste0(wd, "/input/csv/carbon_regimes/", carbon_reg[df], ".csv")
  write.table(csv, file = filename, col.names = F, row.names = F, quote = F, sep = ", ")
}

file.remove(mock)
