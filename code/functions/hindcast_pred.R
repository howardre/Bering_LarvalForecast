# Get predictions per year
hindcast_pred <- function(data, the_year, the_month, doy, gam, roms){
  nlat = 40
  nlon = 60
  latd = seq(min(data$lat), max(data$lat), length.out = nlat)
  lond = seq(min(data$lon), max(data$lon), length.out = nlon)
  
  grid_extent <- expand.grid(lond, latd)
  names(grid_extent) <- c('lon', 'lat')
  
  grid_extent$dist <- NA
  
  for (k in 1:nrow(grid_extent)) {
    dist <- distance_function(grid_extent$lat[k],
                              grid_extent$lon[k],
                              data$lat,
                              data$lon)
    grid_extent$dist[k] <- min(dist)
  }
  
  grid_extent$year <- the_year
  grid_extent$doy <- rep(doy, length(grid_extent))
  grid_extent$month <- the_month
  grid_extent$mean_temp <- roms$mean[roms$year == 2014]
  
  grid_extent[, c(8, 9)] <- as.data.frame(RANN::nn2(data[, c('lat', 'lon')],
                                                    grid_extent[, c('lat', 'lon')],
                                                    k = 1))
  grid_extent$roms_temperature <-data[c(grid_extent$nn.idx), 10] 
  grid_extent$roms_salinity <- data[c(grid_extent$nn.idx), 11] 
  grid_extent <- grid_extent[-c(8, 9)] 
  grid_extent$pred <- exp(predict(gam[[2]],
                                  newdata = grid_extent,
                                  type = "link"))
  grid_extent$pred[grid_extent$dist > 30000] <- NA
  grid_extent$roms_temperature[grid_extent$dist > 30000] <- NA
  return(grid_extent)
}

# Predict all the year
hindcast_loop <- function(range, data, doy, month, gam, roms, temp_range){
  grids <- list()
  for(j in range) {
    grid <- hindcast_pred(data, j, month, doy, gam, roms)
    grids[[paste("year", j, sep = "")]] <- grid
  }
  df <- data.frame(lat = grids[[1]]$lat,
                   lon = grids[[1]]$lon,
                   temperature = rowMeans(do.call(cbind, lapply(grids, "[", "roms_temperature"))),
                   avg_pred = rowMeans(do.call(cbind, lapply(grids, "[", "pred"))))
  df$pred_scaled <- rescale(df$avg_pred)
  return(df)
}
