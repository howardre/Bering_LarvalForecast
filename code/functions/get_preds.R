# New method to use bias corrected ROMS output
get_preds <- function(data, the_year, doy,
                      the_month, proj, temp_output, 
                      salt_output, formula, temp_range){
  nlat = 70
  nlon = 95
  latd = seq(min(data$lat), max(data$lat), length.out = nlat)
  lond = seq(min(data$lon), max(data$lon), length.out = nlon)
  
  grid_extent <- expand.grid(lond, latd)
  names(grid_extent) <- c('lon', 'lat')
  
  # Calculate distance of each grid point to closest 'positive observation'
  grid_extent$dist <- NA
  
  for (k in 1:nrow(grid_extent)) {
    dist <- distance_function(grid_extent$lat[k],
                              grid_extent$lon[k],
                              data$lat,
                              data$lon)
    grid_extent$dist[k] <- min(dist)
  }
  
  # Assign a within sample year and doy to the grid data
  grid_extent$year <- the_year
  grid_extent$doy <- rep(doy, length(grid_extent))
  grid_extent$month <- the_month
  
  temp_output <- temp_output %>% 
    mutate(lat = lat,
           lon = case_when(lon >= 180 ~ lon - 360,
                           lon < 180 ~ lon))
  salt_output <- salt_output %>% 
    mutate(lat = lat,
           lon = case_when(lon >= 180 ~ lon - 360,
                           lon < 180 ~ lon))
  
  
  # Use RANN package to match nearest temperature value based on lat, lon, month, year
  bc_temps <- temp_output %>% filter(month == the_month & year == the_year & projection == proj)
  bc_salts <- salt_output %>% filter(month == the_month & year == the_year & projection == proj)
  
  grid_extent[, c(7, 8)] <- as.data.frame(RANN::nn2(bc_temps[, c('lat', 'lon')],
                                                    grid_extent[, c('lat', 'lon')],
                                                    k = 1))
  grid_extent$roms_temperature <- bc_temps[c(grid_extent$nn.idx), 10] # Match nearest temp
  grid_extent$roms_salinity <- bc_salts[c(grid_extent$nn.idx), 10] # Match nearest salinity
  grid_extent <- grid_extent[-c(6:8)] # remove extra columns before predicting
  
  # Calculate mean temperature
  temp_filtered <- temp_output %>% filter(lon >= -170 & lon <= -165, 
                                          lat >= 56 & lat <= 58,
                                          month >= 2 & month <= 4,
                                          year == the_year)
  mean <- mean(temp_filtered$bc, na.rm = T)
  grid_extent$mean_temp <- mean
  
  # Parameterized model
  gam <- formula
  
  # Predict on forecasted output
  grid_extent$pred <- exp(predict(gam,
                                  newdata = grid_extent,
                                  type = "link",
                                  exclude = "s(year)"))
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

# New function to loop through years
pred_loop <- function(range, data, doy, month, 
                      proj, temp_output, salt_output,
                      the_formula, temp_range){
  grids <- list()
  for(j in range) {
    grid <- get_preds(data, j, doy, month,
                      proj, temp_output, salt_output,
                      the_formula, temp_range)
    grids[[paste("year", j, sep = "")]] <- grid
  }
  return(grids)
}

# sandbox for COG
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