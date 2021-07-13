# Title: Forecasting
# Date: complete 5/25/2021

### Libraries, functions, and data ----
library(maps)
library(maptools)
library(marmap)
library(raster)
library(ncdf4)
library(spacetime)
library(fields)
library(here)
library(tidyverse)
library(lubridate)
library(date)
library(rgdal)
source(here('code/functions', 'distance_function.R'))

# Using avg surface temperatures & salinity
bering_model_temp <- nc_open(here('data/temperature_netcdf', 
                                   'B10K-K20_CORECFS_2020-2024_average_temp_surface5m.nc'))
bering_model_salt <- nc_open(here('data/salinity_netcdf', 
                                   'B10K-K20_CORECFS_2020-2024_average_salt_surface5m.nc'))

# fish data
yfs_egg <- readRDS(here('data', 'yfs_egg.rds'))
yfs_larvae <- readRDS(here('data', 'yfs_larvae.rds'))
akp_egg <- readRDS(here('data', 'akp_egg.rds'))
akp_larvae <- readRDS(here('data', 'akp_larvae.rds'))
fhs_egg <- readRDS(here('data', 'fhs_egg.rds'))
fhs_larvae <- readRDS(here('data', 'fhs_larvae.rds'))
pk_egg <- readRDS(here('data', 'pk_egg.rds'))
pk_larvae <- readRDS(here('data', 'pk_larvae.rds'))

# Function to extract data from .nc files
nc_extract <- function(file, variable, varid_name){
  lon <- ncvar_get(file, varid = 'lon_rho')
  lon1 <- ifelse(lon >= 180, lon -360, lon)
  lat <- ncvar_get(file, varid = 'lat_rho')
  time <- ncvar_get(file, varid = 'ocean_time')
  time1 <- as.Date(time / (60 * 60 * 24), origin = "1900-01-01 00:00:00")
  fillvalue_t <- ncatt_get(file, varid_name, "_FillValue")
  variable <- ncvar_get(file, varid = varid_name)
  variable[variable == fillvalue_t$value] <- NA
  return_list <- list("lon1" = lon1, "lat" = lat, "time1" = time1, "variable" = variable)
} 

# create lists of the four variables, one for salinity and temperature each
temp_output <- nc_extract(bering_model_temp, temp, 'temp')
salt_output <- nc_extract(bering_model_salt, salt, 'salt')

# Close netcdfs
nc_close(bering_model_temp)
nc_close(bering_model_salt)

### Predict future distributions ----
#### Pollock ----
# Two-part binomial and gaussian
# binomial
pk_egg$presence <- 1 * (pk_egg$count > 0)
pk_egg_gam1 <- gam(presence ~ 
                     s(doy) +
                     s(lon, lat),
                   data = pk_egg,
                   family = "binomial")


# gaussian
pk_egg_gam2 <- gam(log(larvalcatchper10m2 + 1) ~ 
                     s(doy, k = 4) +
                     s(lon, lat) +
                     s(roms_salinity, k = 4) +
                     s(roms_temperature, k = 4),
                   data = pk_egg[pk_egg$larvalcatchper10m2 > 0, ])

varid_match <- function(data, model_output1, model_output2){
  data$roms_date <- NA
  data$roms_temperature <- NA
  data$roms_salinity <- NA
  for (i in 1:nrow(data)) {
    idx_time <- order(abs(model_output1[[3]] - data$date[i]))[1]
    data$roms_date[i] <- model_output1[[3]][idx_time]
    idx_grid <- order(distance_function(
      data$lat[i],
      data$lon[i],
      c(model_output1[[2]]),
      c(model_output1[[1]])
    ))[1]
    data$roms_temperature[i] <- c(model_output1[[4]][, , idx_time])[idx_grid]
    data$roms_salinity[i] <- c(model_output2[[4]][, , idx_time])[idx_grid]
  }
  return(data)
}

# Prediction grid
nlat = 80
nlon = 120
latd = seq(min(pk_egg$lat), max(pk_egg$lat), length.out = nlat)
lond = seq(min(pk_egg$lon), max(pk_egg$lon), length.out = nlon)
grid_extent_eggpk <- expand.grid(lond, latd)
names(grid_extent_eggpk) <- c('lon', 'lat')

# Calculate distance of each grid point to closest 'positive observation'
grid_extent_eggpk$dist <- NA
for (k in 1:nrow(grid_extent_eggpk)) {
  dist <- distance_function(grid_extent_eggpk$lat[k],
                            grid_extent_eggpk$lon[k],
                            pk_egg$lat,
                            pk_egg$lon)
  grid_extent_eggpk$dist[k] <- min(dist)
}

# Assign a within sample year and doy to the grid data
grid_extent_eggpk$year <- 2020
grid_extent_eggpk$date <- rep(as.Date("2022-02-18"), length(grid_extent_eggpk))
grid_extent_eggpk$doy <- rep(49, length(grid_extent_eggpk))

# Attach ROMs forecast
grid_extent_eggpk <- varid_match(grid_extent_eggpk, temp_output, salt_output)

grid_extent_eggpk$pred1 <- predict(pk_egg_gam1, newdata = grid_extent_eggpk)
grid_extent_eggpk$pred2 <- predict(pk_egg_gam2, newdata = grid_extent_eggpk)
grid_extent_eggpk$pred <- grid_extent_eggpk$pred1 * grid_extent_eggpk$pred2
grid_extent_eggpk$pred[grid_extent_eggpk$dist > 30000] <- NA

# Plot
windows(width = 6, height = 5)
par(mfrow = c(1, 1), 
    mai = c(0.8, 0.9, 0.5, 0.5))
image(lond,
      latd,
      t(matrix(grid_extent_eggpk$pred,
               nrow = length(latd),
               ncol = length(lond),
               byrow = T)),
      col = viridis(100, option = "F", direction = 1),
      ylab = "Latitude",
      xlab = "Longitude",
      xlim = c(-176.5, -156.5),
      ylim = c(52, 62),
      main = 'Forecasted Distribution',
      cex.main = 1.2,
      cex.lab = 1.1,
      cex.axis = 1.1)
maps::map("worldHires",
    fill = T,
    col = "wheat4",
    add = T)