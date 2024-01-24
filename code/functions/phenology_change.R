phenology_change <- function(grids, bounds_min, bounds_max, ylab){
  plot(grids[[2]]$doy,
       grids[[2]]$pred,
       main = '',
       type = 'l',
       ylim = range(c(grids[[2]]$pred_up,
                      grids[[2]]$pred2_up,
                      grids[[2]]$pred_lw,
                      grids[[2]]$pred2_lw), na.rm = T),
       xlim = c(min(bounds_min, na.rm = T), 
                max(bounds_max, na.rm = T)),
       col = '#000000',
       lwd = 2.5,
       xlab = 'Day of the year',
       ylab = ylab,
       cex.lab = 6.5,
       cex.axis = 6.5)
  polygon(c(grids[[2]]$doy, rev(grids[[2]]$doy)),
          c(grids[[2]]$pred_lw, rev(grids[[2]]$pred_up)),
          col = alpha('#000000', 0.2), 
          lty = 0)
  polygon(c(grids[[2]]$doy, rev(grids[[2]]$doy)),
          c(grids[[2]]$pred2_lw, rev(grids[[2]]$pred2_up)),
          col = alpha('#1C0D51', 0.6),
          lty = 0)
  lines(grids[[2]]$doy,
        grids[[2]]$pred2,
        col = '#000000',
        lwd = 2.5)
}
