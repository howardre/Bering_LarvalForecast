# Get predictions per year
hindcast_pred <- function(data, the_year, the_month, doy, gam, roms, temp_range){
  nlat = 70
  nlon = 95
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
  grid_extent$mean_temp <- roms$mean[roms$year == the_year]
  
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
  grid_extent$pred_scaled <- rescale(grid_extent$pred)
  grid_extent$temperature <- ifelse(between(grid_extent$roms_temperature, 
                                            min(temp_range$x),
                                            max(temp_range$x)), 1, 0)
  COG_year <- COG_pred(grid_extent)
  COG_temp <- COG_temp(grid_extent)
  grid_extent <- subset(grid_extent, select = -c(pred_scaled, temperature))
  final_list <- list(grid_extent, COG_year, COG_temp)
  return(final_list)
}


hindcast_cells <- function(range, data, doy, month, gam, roms, temp_range){
  preds <- hindcast_loop(range, data, doy, month, gam, roms, temp_range)
  new_list <- lapply(preds, "[[", 1)
  COG_abun <- lapply(preds, "[[", 2)
  COG_temp <- lapply(preds, "[[", 3)
  df <- data.frame(lat = new_list[[1]]$lat,
                   lon = new_list[[1]]$lon,
                   avg_pred = rowMeans(do.call(cbind, lapply(new_list, "[", "pred"))))
  df$pred_scaled <- rescale(df$avg_pred)
  final_list <- list(df, COG_abun, COG_temp)
  return(final_list)
}

# Predict all the year
hindcast_loop <- function(range, data, doy, month, gam, roms, temp_range){
  grids <- list()
  for(j in range) {
    grid <- hindcast_pred(data, j, month, doy, gam, roms, temp_range)
    grids[[paste("year", j, sep = "")]] <- grid
  }
  return(grids)
}

# Calculate COG for predictions
COG_pred <- function(data){
  temp1 <- data %>% dplyr::mutate(pred = pred * lat) %>%
    dplyr::summarize(value = sum(pred, na.rm = TRUE)) %>%
    dplyr::select(value)
  
  temp2 <- data %>%
    dplyr::summarize(value = sum(pred, na.rm = TRUE)) %>%
    dplyr::select(value)
  
  COG_x <- temp1 / temp2
  
  temp3 <- data %>% dplyr::mutate(pred = pred * lon) %>%
    dplyr::summarize(value1 = sum(pred, na.rm = TRUE)) %>%
    dplyr::select(value1)
  
  temp4 <- data %>%
    dplyr::summarize(value1 = sum(pred, na.rm = TRUE)) %>%
    dplyr::select(value1)
  
  COG_y <- temp3 / temp4
  
  COG <- as.data.frame(cbind(COG_x, COG_y))
  return(COG)
}

# Calculate COG for temperature
COG_temp <- function(data){
  temp1 <- data %>% dplyr::mutate(temp = temperature * lat) %>%
    dplyr::summarize(value = sum(temp, na.rm = TRUE)) %>%
    dplyr::select(value)
  
  temp2 <- data %>%
    dplyr::summarize(value = sum(temperature, na.rm = TRUE)) %>%
    dplyr::select(value)
  
  COG_x <- temp1 / temp2
  
  temp3 <- data %>% dplyr::mutate(temp = temperature * lon) %>%
    dplyr::summarize(value1 = sum(temp, na.rm = TRUE)) %>%
    dplyr::select(value1)
  
  temp4 <- data %>%
    dplyr::summarize(value1 = sum(temperature, na.rm = TRUE)) %>%
    dplyr::select(value1)
  
  COG_y <- temp3 / temp4
  
  COG <- as.data.frame(cbind(COG_x, COG_y))
  return(COG)
}