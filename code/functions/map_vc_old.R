map_vc <- function(data, gam, grids, title){
  nlat = 70
  nlon = 90
  latd = seq(min(data$lat), max(data$lat), length.out = nlat)
  lond = seq(min(data$lon), max(data$lon), length.out = nlon)
  myvis_gam(gam,
            view = c('lon', 'lat'),
            too.far = 0.04,
            plot.type = 'contour',
            contour.col = contour_col,
            color = "jet",
            type = 'link',
            xlim = c(-176.5, -156.5),
            ylim = c(52, 62),
            xlab = "",
            ylab = "",
            main = "",
            axes = FALSE)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "mintcream")
  par(new = TRUE)
  myvis_gam(gam,
            view = c('lon', 'lat'),
            too.far = 0.04,
            plot.type = 'contour',
            contour.col = contour_col,
            color = "jet" ,
            type = 'link',
            xlim = c(-176.5, -156.5),
            ylim = c(52, 62),
            family = "serif",
            xlab = "Longitude",
            ylab = "Latitude",
            main = title,
            cex.main = 1.5,
            cex.lab = 1.5,
            cex.axis = 1.5)
  maps::map('worldHires',
            add = T,
            col = 'wheat4',
            fill = T)
  image.plot(legend.only = T,
             col = jet.colors(100),
             legend.shrink = 0.2,
             smallplot = c(.72, .75, .24, .38),
             legend.cex = 0.7,
             axis.args = list(cex.axis = 0.9),
             legend.width = 0.8,
             legend.mar = 6,
             zlim = c(min(gam$linear.predictors), 
                      max(gam$linear.predictors)),
             legend.args = list("Occurrence",
                                side = 2, cex = 0.8))
  my_color = colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D",
                                    "#F4A582", "#FDDBC7", "#F7F7F7",
                                    "#D1E5F0", "#92C5DE", "#4393C3",
                                    "#2166AC", "#053061")))
  color_levels = 100
  max_absolute_value = max(abs(c(min(grids[[1]]$diff, na.rm = T), 
                                 max(grids[[1]]$diff, na.rm = T)))) 
  color_sequence = seq(-max_absolute_value, max_absolute_value, 
                       length.out = color_levels + 1)
  n_in_class = hist(grids[[1]]$diff, breaks = color_sequence, plot = F)$counts > 0
  col_to_include = min(which(n_in_class == T)):max(which(n_in_class == T))
  breaks_to_include = min(which(n_in_class == T)):(max(which(n_in_class == T)) + 1)
  image(lond,
        latd,
        t(matrix(grids[[1]]$diff,
                 nrow = length(latd),
                 ncol = length(lond),
                 byrow = T)),
        xlim = c(-176.5, -156.5),
        ylim = c(52, 62),
        axes = FALSE,
        xlab = "",
        ylab = "")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "mintcream")
  par(new = TRUE)
  image(lond,
        latd,
        t(matrix(grids[[1]]$diff,
                 nrow = length(latd),
                 ncol = length(lond),
                 byrow = T)),
        col = my_color(n = color_levels)[col_to_include],
        breaks = color_sequence[breaks_to_include],
        ylab = "Longitude",
        xlab = "Latitude",
        xlim = c(-176.5, -156.5),
        ylim = c(52, 62),
        main = 'Change in distribution',
        cex.main = 1.5,
        cex.lab = 1.5,
        cex.axis = 1.5)
  maps::map("worldHires",
            fill = T,
            col = "wheat4",
            add = T)
  image.plot(legend.only = T,
             col = my_color(n = color_levels)[col_to_include],
             legend.shrink = 0.2,
             smallplot = c(.73, .76, .22, .36),
             legend.cex = 0.7,
             axis.args = list(cex.axis = 0.9),
             legend.width = 0.8,
             legend.mar = 6,
             zlim = c(min(grids[[1]]$diff, na.rm = T), 
                      max(grids[[1]]$diff, na.rm = T)),
             legend.args = list("Predicted \n Change",
                                side = 2, cex = 0.8))
  plot(grids[[2]]$doy,
       grids[[2]]$pred,
       main = 'Change in phenology',
       type = 'l',
       ylim = range(c(grids[[2]]$pred_up,
                      grids[[2]]$pred2_up,
                      grids[[2]]$pred_lw,
                      grids[[2]]$pred2_lw), na.rm = T),
       xlim = c(min(grids[[2]]$doy, na.rm = T), 
                max(grids[[2]]$doy, na.rm = T)),
       col = '#000000',
       lwd = 2,
       xlab = 'Day of the year',
       ylab = 'Egg density ln(n/10m2)',
       cex.lab = 1.5,
       cex.axis = 1.5,
       cex.main = 1.5)
  polygon(c(grids[[2]]$doy, rev(grids[[2]]$doy)),
          c(grids[[2]]$pred_lw, rev(grids[[2]]$pred_up)),
          col = alpha('#000000', 0.2), 
          lty = 0)
  lines(grids[[2]]$doy,
        grids[[2]]$pred2,
        col = '#E69F00',
        lwd = 2)
  polygon(c(grids[[2]]$doy, rev(grids[[2]]$doy)),
          c(grids[[2]]$pred2_lw, rev(grids[[2]]$pred2_up)),
          col = alpha('#E69F00', 0.2),
          lty = 0)
}
