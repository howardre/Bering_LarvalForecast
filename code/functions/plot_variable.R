plot_variable <- function(gam, covariate, bounds, variable, ylabel, yvalues){
  plot(gam[[2]],
       pages = 0,
       select = covariate, # 1 = year/PDO/NPGO, 2 = lat/lon, 3 = depth, 4 = julian, 5 = temp
       shade = T,
       shade.col = "lemonchiffon3",
       ylim = bounds,
       xlab = variable,
       ylab = ylabel,
       yaxt = yvalues,
       seWithMean = T,
       scale = 0,
       cex.axis = 6,
       cex.lab = 6,
       family = "serif",
       lwd = 2)
}