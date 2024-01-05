threepanel_base_mse <- function(data, base_model) {
  nlat = 70
  nlon = 95
  latd = seq(min(data$lat), max(data$lat), length.out = nlat)
  lond = seq(min(data$lon), max(data$lon), length.out = nlon)
  grid <- expand.grid(lond, latd)
  names(grid) <- c('lon', 'lat')
  
  # Geography grid
  grid$dist <- NA
  for (k in 1:nrow(grid)) {
    dist <-
      distance_function(grid$lat[k],
                        grid$lon[k],
                        data$lat,
                        data$lon)
    grid$dist[k] <- min(dist)
  }
  grid$year <- 2014
  grid$doy <- median(data$doy)
  grid$pred <- predict(base_model, newdata = grid)
  grid$pred[grid$dist > 30000] <- NA
  
  # Phenology grid
  grid2 <-  data.frame('lon' = rep(-170, 100),
                       'lat' = rep(57, 100),
                       'doy' = seq(min(data$doy),
                                   max(data$doy),
                                   length = 100),
                       'year' = rep(2014, 100))
  grid2$pred <- predict(base_model, newdata = grid2)
  grid2$se <- predict(base_model, newdata = grid2, se = T)[[2]]
  grid2$pred_up <- grid2$pred + 1.96 * grid2$se
  grid2$pred_lw <- grid2$pred - 1.96 * grid2$se
  return(list(grid, grid2))
}