# Function to make maps
grid_predict <- function(grid, title){
  nlat = 40
  nlon = 60
  latd = seq(min(grid$lat), max(grid$lat), length.out = nlat)
  lond = seq(min(grid$lon), max(grid$lon), length.out = nlon)
  my_color = colorRampPalette(rev(c("#FFFFCC", "#FBF2A8", "#F9E585",
                                    "#F5D363", "#EFBA55", "#EAA352",
                                    "#E68C51", "#E0754F", "#D75C4D",
                                    "#BB4A48", "#994240", "#763931", 
                                    "#542D20", "#352311", "#191900")))
  image(lond,
        latd,
        t(matrix(grid$pred_scaled,
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
        t(matrix(grid$pred_scaled,
                 nrow = length(latd),
                 ncol = length(lond),
                 byrow = T)),
        col = my_color(100), 
        ylab = "Latitude",
        xlab = "Longitude",
        xlim = c(-176.5, -156.5),
        ylim = c(52, 62),
        zlim = c(min(grid$pred_scaled, na.rm = T), 
                 max(grid$pred_scaled, na.rm = T)),
        main = title,
        cex.main = 1.2,
        cex.lab = 1.1,
        cex.axis = 1.1)
  maps::map("worldHires",
            fill = T,
            col = "wheat4",
            add = T)
  image.plot(legend.only = T,
             col = my_color(100),
             legend.shrink = 0.2,
             smallplot = c(.79, .82, .20, .37),
             legend.cex = 0.8,
             axis.args = list(cex.axis = 0.8),
             legend.width = 0.5,
             legend.mar = 6,
             zlim = c(min(grid$pred_scaled, na.rm = T), 
                      max(grid$pred_scaled, na.rm = T)),
             legend.args = list("Avg. Predicted \n Occurrence",
                                side = 2, cex = 1))
}