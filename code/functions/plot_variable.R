plot_variable <- function(gam, covariate, bounds, variable, ylabel, yvalues, bounds_min, bounds_max){
  plot(gam[[2]],
       pages = 0,
       select = covariate, 
       shade = T,
       shade.col = "#eb8055b2",
       xlim = c(bounds_min, bounds_max),
       ylim = bounds,
       xlab = variable,
       ylab = ylabel,
       yaxt = yvalues,
       xaxs = "i",
       seWithMean = T,
       scale = 0,
       cex.axis = 6.5,
       cex.lab = 6.5,
       family = "serif",
       lwd = 2.5)
}
