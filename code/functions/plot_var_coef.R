# Plot variable coefficient term
plot_var_coef <- function(gam, data, predictions){
  jet.colors <- colorRampPalette(c(sequential_hcl(15, palette = "Mint")))
  par(mar = c(6.4, 7.2, .5, 0.6) + 0.1,
      oma = c(1, 1, 1, 1),
      mgp = c(5, 2, 0))
  myvis_gam(gam[[2]],
            view = c('lon', 'lat'),
            too.far = 0.07,
            plot.type = 'contour',
            contour.col = contour_col,
            color = "jet" ,
            type = 'response',
            xlim = c(-176.5, -156.5),
            ylim = c(52, 62),
            family = "serif",
            xlab = "Longitude",
            ylab = "Latitude",
            main = " ",
            cex.lab = 2.5,
            cex.axis =  2.5)
  symbols(data$lon[predictions[[2]]],
          data$lat[predictions[[2]]],
          circle = predictions[[3]][predictions[[2]]],
          inches = 0.12,
          add = T,
          bg = alpha('darkred', 0.4),
          fg = alpha('black', 0.08))
  symbols(data$lon[predictions[[1]]],
          data$lat[predictions[[1]]],
          circle = (-1) * predictions[[3]][predictions[[1]]],
          inches = 0.12,
          add = T,
          bg = alpha('navy', 0.4),
          fg = alpha('black', 0.08))
  maps::map("worldHires",
            fill = T,
            col = "wheat4",
            add = T)
  image.plot(legend.only = T,
             col = jet.colors(100),
             legend.shrink = 0.2,
             smallplot = c(.8, .83, .15, .35),
             legend.cex = 1.3,
             axis.args = list(cex.axis = 1.8,
                              family = "serif"),
             legend.width = 0.8,
             legend.mar = 6,
             zlim = c(min(gam[[2]]$fitted.values),
                      max(gam[[2]]$fitted.values)),
             legend.args = list("Catch",
                                side = 2, cex = 2,
                                family = "serif"))
}

plot_var_coef2 <- function(gam, data, predictions){
  jet.colors <- colorRampPalette(c(sequential_hcl(15, palette = "Mint")))
  par(mar = c(6.4, 7.2, .5, 0.6) + 0.1,
      oma = c(1, 1, 1, 1),
      mgp = c(5, 2, 0))
  myvis_gam(gam[[2]],
            view = c('lon', 'lat'),
            too.far = 0.07,
            plot.type = 'contour',
            contour.col = contour_col,
            color = "jet" ,
            type = 'response',
            xlim = c(-176.5, -156.5),
            ylim = c(52, 62),
            family = "serif",
            xlab = "Longitude",
            ylab = "Latitude",
            main = " ",
            cex.lab = 2.5,
            cex.axis =  2.5)
  symbols(data$lon[predictions[[2]]],
          data$lat[predictions[[2]]],
          circle = predictions[[3]][predictions[[2]]],
          inches = 0.12,
          add = T,
          bg = alpha('darkred', 0.4),
          fg = alpha('black', 0.08))
  # symbols(data$lon[predictions[[1]]],
  #         data$lat[predictions[[1]]],
  #         circle = (-1) * predictions[[3]][predictions[[1]]],
  #         inches = 0.12,
  #         add = T,
  #         bg = alpha('navy', 0.4),
  #         fg = alpha('black', 0.08))
  maps::map("worldHires",
            fill = T,
            col = "wheat4",
            add = T)
  image.plot(legend.only = T,
             col = jet.colors(100),
             legend.shrink = 0.2,
             smallplot = c(.8, .83, .15, .35),
             legend.cex = 1.3,
             axis.args = list(cex.axis = 1.8,
                              family = "serif"),
             legend.width = 0.8,
             legend.mar = 6,
             zlim = c(min(gam[[2]]$fitted.values),
                      max(gam[[2]]$fitted.values)),
             legend.args = list("Catch",
                                side = 2, cex = 2,
                                family = "serif"))
}


plot_var_coef3 <- function(gam, data, predictions){
  jet.colors <- colorRampPalette(c(sequential_hcl(15, palette = "Mint")))
  par(mar = c(6.4, 7.2, .5, 0.6) + 0.1,
      oma = c(1, 1, 1, 1),
      mgp = c(5, 2, 0))
  myvis_gam(gam[[2]],
            view = c('lon', 'lat'),
            too.far = 0.07,
            plot.type = 'contour',
            contour.col = contour_col,
            color = "jet" ,
            type = 'response',
            xlim = c(-176.5, -156.5),
            ylim = c(52, 62),
            family = "serif",
            xlab = "Longitude",
            ylab = "Latitude",
            main = " ",
            cex.lab = 2.5,
            cex.axis =  2.5)
  # symbols(data$lon[predictions[[2]]],
  #         data$lat[predictions[[2]]],
  #         circle = predictions[[3]][predictions[[2]]],
  #         inches = 0.12,
  #         add = T,
  #         bg = alpha('darkred', 0.4),
  #         fg = alpha('black', 0.08))
  symbols(data$lon[predictions[[1]]],
          data$lat[predictions[[1]]],
          circle = (-1) * predictions[[3]][predictions[[1]]],
          inches = 0.12,
          add = T,
          bg = alpha('navy', 0.4),
          fg = alpha('black', 0.08))
  maps::map("worldHires",
            fill = T,
            col = "wheat4",
            add = T)
  image.plot(legend.only = T,
             col = jet.colors(100),
             legend.shrink = 0.2,
             smallplot = c(.8, .83, .15, .35),
             legend.cex = 1.3,
             axis.args = list(cex.axis = 1.8,
                              family = "serif"),
             legend.width = 0.8,
             legend.mar = 6,
             zlim = c(min(gam[[2]]$fitted.values),
                      max(gam[[2]]$fitted.values)),
             legend.args = list("Catch",
                                side = 2, cex = 2,
                                family = "serif"))
}

