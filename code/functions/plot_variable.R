plot_variable <- function(gam, covariate, bounds, variable, ylabel, yvalues){
  plot(gam[[2]],
       pages = 0,
       select = covariate, 
       shade = T,
       shade.col = "lightsalmon2",
       ylim = bounds,
       xlab = variable,
       ylab = ylabel,
       yaxt = yvalues,
       seWithMean = T,
       scale = 0,
       cex.axis = 6,
       cex.lab = 6,
       family = "serif",
       lwd = 2.5)
}

plot_variable_pk <- function(gam, covariate, bounds, variable, ylabel, yvalues){
  plot(gam[[2]],
       pages = 0,
       select = covariate, 
       shade = T,
       shade.col = "lightsalmon2",
       yaxp = c(-4, 8, 3),
       ylim = bounds,
       xlab = variable,
       ylab = ylabel,
       yaxt = yvalues,
       seWithMean = T,
       scale = 0,
       cex.axis = 6,
       cex.lab = 6,
       family = "serif",
       lwd = 2.5)
}
