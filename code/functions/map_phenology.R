map_phenology <- function(data, grids, ratio_base_geog,
                          vc_ratio_space, ratio_base_pheno,
                          vc_ratio_time){
  nlat = 70
  nlon = 95
  latd = seq(min(grids[[1]]$lat), max(grids[[1]]$lat), length.out = nlat)
  lond = seq(min(grids[[1]]$lon), max(grids[[1]]$lon), length.out = nlon)
  my_color = colorRampPalette(c("#1C0D51", "#4C408E", "#7E77B0",
                                "#AFABCB", "#DAD9E5", "#F9F9F9",
                                "#FFDAB7", "#FFB377","#E18811",
                                "#AC6000", "#743700"))
  color_levels = 100
  max_absolute_value = max(abs(c(min(grids[[1]]$pred, na.rm = T), max(grids[[1]]$pred, na.rm = T)))) 
  color_sequence = seq(-max_absolute_value, max_absolute_value, 
                       length.out = color_levels + 1)
  n_in_class = hist(grids[[1]]$pred, breaks = color_sequence, plot = F)$counts > 0
  col_to_include = min(which(n_in_class == T)):max(which(n_in_class == T))
  breaks_to_include = min(which(n_in_class == T)):(max(which(n_in_class == T)) + 1)
  image(lond,
        latd,
        t(matrix(grids[[1]]$pred,
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
        t(matrix(grids[[1]]$pred,
                 nrow = length(latd),
                 ncol = length(lond),
                 byrow = T)),
        col = my_color(n = color_levels)[col_to_include], 
        breaks = color_sequence[breaks_to_include],
        ylab = "Latitude",
        xlab = "Longitude",
        xlim = c(-176.5, -156.5),
        ylim = c(52, 62),
        main = "Distribution",
        cex.main = 1.2,
        cex.lab = 1.1,
        cex.axis = 1.1)
  symbols(data$lon[data$larvalcatchper10m2 > 0],
          data$lat[data$larvalcatchper10m2 > 0],
          circles = log(data$larvalcatchper10m2 + 1)[data$larvalcatchper10m2 > 0],
          inches = 0.1,
          bg = alpha('grey', 0.1),
          fg = alpha('black', 0.05),
          add = T)
  points(data$lon[data$larvalcatchper10m2 == 0], 
         data$lat[data$larvalcatchper10m2 == 0], 
         pch = '')
  maps::map("worldHires",
            fill = T,
            col = "wheat4",
            add = T)
  image.plot(legend.only = T,
             col = my_color(n = color_levels)[col_to_include],
             legend.shrink = 0.2,
             smallplot = c(.765, .79, .25, .37),
             legend.cex = 0.4,
             axis.args = list(cex.axis = 0.8),
             legend.width = 0.5,
             legend.mar = 6,
             zlim = c(min(grids[[1]]$pred, na.rm = T), max(grids[[1]]$pred, na.rm = T)),
             legend.args = list("Predicted \n Change",
                                side = 2, cex = 0.7))
  plot(grids[[2]]$doy,
       grids[[2]]$pred,
       main = 'Phenology',
       type = 'l',
       ylab = 'Egg density ln(n/10m2)',
       xlab = 'Day of the year',
       cex.lab = 1.1,
       cex.axis = 1.1,
       cex.main = 1.2,
       xlim = c(min(grids[[2]]$doy, na.rm = T), max(grids[[2]]$doy, na.rm = T)),
       ylim = range(c(grids[[2]]$pred_up, grids[[2]]$pred_lw)),
       col = 'blue',
       lwd = 2)
  polygon(c(grids[[2]]$doy, rev(grids[[2]]$doy)),
          c(grids[[2]]$pred_lw, rev(grids[[2]]$pred_up)),
          col = alpha('gray', 0.6),
          lty = 0)
  barplot(matrix(c(ratio_base_geog,
                   vc_ratio_space,
                   ratio_base_pheno,
                   vc_ratio_time) *  100,
                 ncol = 2,
                 nrow = 2),
          ylab = 'Delta MSE (%)',
          names.arg = c('Variable distribution', 'Variable phenologies'),
          ylim = c(0, 60),
          main = 'Delta MSE',
          cex.lab = 1.1,
          cex.main = 1.2,
          cex.axis = 1.1,
          cex.names = 1.1,
          space = c(0, 1),
          beside = T,
          col = c('azure4', 'azure3'))
  box()
  legend("topright",
         legend = c('D-MSE|Year', 'D-MSE|SST'),
         bty = 'n',
         col = c('azure4', 'azure3'),
         pch = 15,
         pt.cex = 2.5,
         cex = 1.1)
}