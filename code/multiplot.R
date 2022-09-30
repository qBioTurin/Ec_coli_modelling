# Multiplot creates a layout structure for multiple plots in one
# window. The function accepts ggplot objects via "..." or as a
# list of objects. Use input "cols" to define the number of
# columns in the plot layout. Alternatively, input "layout" to
# define the row/col matrix of plot elements. If "layout" is
# present, "cols" is ignored. For example, given a layout defined
# by matrix(c(1,2,3,3), nrow=2, byrow=TRUE), then plot #1 goes in
# the upper left, #2 goes in the upper right, and #3 will cross
# the entire bottom row, utilizing both columns.

multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
   # launch the grid graphical system
   require(grid)

   # make a list from the "..." arguments and/or plotlist
   plots <- c(list(...), plotlist)
   numPlots = length(plots)

   # if layout = NULL, then use 'cols' to determine layout
   if (is.null(layout)) {
      # make the panel using ncol; nrow is calculated from ncol
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)), ncol = cols, nrow = ceiling(numPlots/cols))
   }

   if (numPlots==1) {
      print(plots[[1]])
   } else {
      # set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
      # put each plot, in the correct location
      for (i in 1:numPlots) {
         # get the i,j matrix position of the subplot
         matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
         print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
      }
   }
}
