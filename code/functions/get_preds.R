# New method to use bias corrected ROMS output
get_preds <- function(data, the_year, doy,
                      the_month, proj, temp_output, 
                      salt_output, formula){
  nlat = 40
  nlon = 60
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
  grid_extent$roms_salinity <- bc_salts[c(grid_extent$nn.idx), 10] # Match nearest temp
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
  return(grid_extent)
}


# New function to loop through years
pred_loop <- function(range, data, doy, month, 
                      proj, temp_output, salt_output,
                      the_formula){
  grids <- list()
  for(j in range) {
    grid <- get_preds(data, j, doy, month,
                      proj, temp_output, salt_output,
                      the_formula)
    grids[[paste("year", j, sep = "")]] <- grid
  }
  return(grids)
}
