plot_variable <- function(gam, covariate, bounds, variable, ylabel, yvalues){
  par(
    mar = c(6.4, 7.2, .5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0))
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
       cex.axis = 3,
       cex.lab = 3,
       family = "serif",
       lwd = 2)
}