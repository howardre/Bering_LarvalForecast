# Create grids with VC GAMs
vc_grids <- function(data, gam_list){
  nlat = 70
  nlon = 90
  latd = seq(min(data$lat), max(data$lat), length.out = nlat)
  lond = seq(min(data$lon), max(data$lon), length.out = nlon)
  spatial_grid <- expand.grid(lond, latd)
  names(spatial_grid) <- c('lon', 'lat')
  spatial_grid$dist <- NA
  for (k in 1:nrow(spatial_grid)) {
    dist <-  distance_function(spatial_grid$lat[k],
                               spatial_grid$lon[k],
                               data$lat,
                               data$lon)
    spatial_grid$dist[k] <- min(dist)
  }
  spatial_grid$year <- 2014
  spatial_grid$doy <- median(data$doy, na.rm = T)
  spatial_grid$mean_temp <- mean(data$mean_temp, na.rm = T)
  spatial_grid$roms_temperature <- mean(data$roms_temperature, na.rm = T)
  spatial_grid$roms_salinity <- mean(data$roms_salinity, na.rm = T)
  spatial_grid$pred <- predict(gam_list[[2]], newdata = spatial_grid)
  spatial_grid$se <- predict(gam_list[[2]], newdata = spatial_grid, se = T)[[2]]
  spatial_grid$pred_up <- spatial_grid$pred + 1.96 * spatial_grid$se
  spatial_grid$pred_lw <- spatial_grid$pred - 1.96 * spatial_grid$se
  spatial_grid$pred[spatial_grid$dist > 30000] <- NA
  spatial_grid$mean_temp <- mean(data$mean_temp, na.rm = T) + 1
  spatial_grid$pred2 <- predict(gam_list[[2]], newdata = spatial_grid)
  spatial_grid$se2 <- predict(gam_list[[2]], newdata = spatial_grid, se = T)[[2]]
  spatial_grid$pred2_up <- spatial_grid$pred2 + 1.96 * spatial_grid$se2
  spatial_grid$pred2_lw <- spatial_grid$pred2 - 1.96 * spatial_grid$se2
  spatial_grid$diff <- spatial_grid$pred2 - spatial_grid$pred
  spatial_grid$sig_pos <- c(spatial_grid$pred2_lw > spatial_grid$pred_up)
  spatial_grid$sig_neg <- c(spatial_grid$pred2_up < spatial_grid$pred_lw)
  spatial_grid$pos_diff <- spatial_grid$diff * spatial_grid$sig_pos
  spatial_grid$neg_diff <- spatial_grid$diff * spatial_grid$sig_neg
  max_slope <- max(spatial_grid$diff, na.rm = T)
  
  phenology_grid <-  data.frame('lon' = rep(-170, 100),
                                'lat' = rep(57, 100),
                                'doy' = seq(min(data$doy, na.rm = T),
                                            max(data$doy, na.rm = T),
                                            length = 100),
                                mean_temp = rep(mean(data$mean_temp, na.rm = T), 100),
                                'year' = rep(2014, 100),
                                'roms_temperature' = rep(mean(data$roms_temperature, na.rm = T), 100),
                                'roms_salinity' = rep(mean(data$roms_salinity, na.rm = T), 100))
  phenology_grid$pred <- predict(gam_list[[1]][[2]], newdata = phenology_grid)
  phenology_grid$se <- predict(gam_list[[1]][[2]], newdata = phenology_grid, se = T)[[2]]
  phenology_grid$pred_up <- phenology_grid$pred + 1.96 * phenology_grid$se
  phenology_grid$pred_lw <- phenology_grid$pred - 1.96 * phenology_grid$se
  phenology_grid$mean_temp <- rep(mean(data$mean_temp, na.rm = T), 100) + 1
  phenology_grid$pred2 <- predict(gam_list[[1]][[2]], newdata = phenology_grid)
  phenology_grid$se2 <- predict(gam_list[[1]][[2]], newdata = phenology_grid, se = T)[[2]]
  phenology_grid$pred2_up <- phenology_grid$pred2 + 1.96 * phenology_grid$se2
  phenology_grid$pred2_lw <- phenology_grid$pred2 - 1.96 * phenology_grid$se2
  return(list(spatial_grid, phenology_grid))
}
