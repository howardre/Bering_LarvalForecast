# Title: Forecasting
# Date: complete 7/07/2021

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
library(RColorBrewer)
library(mgcv)
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

# Match ROMS output function
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

# Function to make maps
grid_predict <- function(grid, title){
  nlat = 80
  nlon = 120
  latd = seq(min(grid$lat), max(grid$lat), length.out = nlat)
  lond = seq(min(grid$lon), max(grid$lon), length.out = nlon)
  my_color = colorRampPalette(c("#1C0D51", "#4C408E", "#7E77B0",
                                "#AFABCB", "#DAD9E5", "#F9F9F9",
                                "#FFDAB7", "#FFB377","#E18811",
                                "#AC6000", "#743700"))
  color_levels = 100
  max_absolute_value = max(abs(c(min(grid$pred, na.rm = T), max(grid$pred, na.rm = T)))) 
  color_sequence = seq(-max_absolute_value, max_absolute_value, 
                       length.out = color_levels + 1)
  n_in_class = hist(grid$pred, breaks = color_sequence, plot = F)$counts > 0
  col_to_include = min(which(n_in_class == T)):max(which(n_in_class == T))
  breaks_to_include = min(which(n_in_class == T)):(max(which(n_in_class == T)) + 1)
  image(lond,
        latd,
        t(matrix(grid$pred,
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
        t(matrix(grid$pred,
                 nrow = length(latd),
                 ncol = length(lond),
                 byrow = T)),
        col = my_color(n = color_levels)[col_to_include], 
        breaks = color_sequence[breaks_to_include],
        ylab = "Latitude",
        xlab = "Longitude",
        xlim = c(-176.5, -156.5),
        ylim = c(52, 62),
        main = title,
        cex.main = 1.2,
        cex.lab = 1.1,
        cex.axis = 1.1)
  maps::map("worldHires",
            fill = T,
            col = "wheat4",
            add = T)
  image.plot(legend.only = T,
             col = my_color(n = color_levels)[col_to_include],
             legend.shrink = 0.2,
             smallplot = c(.79, .82, .20, .37),
             legend.cex = 0.8,
             axis.args = list(cex.axis = 0.8),
             legend.width = 0.5,
             legend.mar = 6,
             zlim = c(min(grid$pred, na.rm = T), max(grid$pred, na.rm = T)),
             legend.args = list("Predicted \n Change",
                                side = 2, cex = 1))
}

# create lists of the four variables, one for salinity and temperature each
temp_output <- nc_extract(bering_model_temp, temp, 'temp')
salt_output <- nc_extract(bering_model_salt, salt, 'salt')

# Close netcdfs
nc_close(bering_model_temp)
nc_close(bering_model_salt)

### Predict future distributions ----
#### Pollock ----
##### Eggs ----
# Two-part binomial and gaussian
# binomial
pk_egg$presence <- 1 * (pk_egg$count > 0)
pk_egg_gam1 <- gam(presence ~ 
                     s(doy, k = 4)
                   data = pk_egg,
                   family = "binomial")
plot(pk_egg_gam1)


# gaussian
pk_egg_gam2 <- gam(log(larvalcatchper10m2 + 1) ~ 
                     s(doy, k = 4) +
                     s(lon, lat) +
                     s(roms_salinity, k = 4) +
                     s(roms_temperature, k = 4),
                   data = pk_egg[pk_egg$larvalcatchper10m2 > 0, ])
par(mfrow = c(2, 2))
plot(pk_egg_gam2)


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
grid_extent_eggpk$year <- 2022
grid_extent_eggpk$date <- rep(as.Date("2022-02-18"), length(grid_extent_eggpk))
grid_extent_eggpk$doy <- rep(49, length(grid_extent_eggpk))

# Attach ROMs forecast
grid_extent_eggpk <- varid_match(grid_extent_eggpk, temp_output, salt_output)

# Predict on forecasted output
grid_extent_eggpk$pred1 <- predict(pk_egg_gam1, newdata = grid_extent_eggpk)
grid_extent_eggpk$pred2 <- predict(pk_egg_gam2, newdata = grid_extent_eggpk)
grid_extent_eggpk$pred <- grid_extent_eggpk$pred1 * grid_extent_eggpk$pred2
grid_extent_eggpk$pred[grid_extent_eggpk$dist > 30000] <- NA

# Plot
windows(width = 6, height = 5, family = "serif")
grid_predict(grid_extent_eggpk, "Forecasted Distribution")

##### Larvae ----
# Two-part binomial and gaussian
# binomial
pk_larvae$presence <- 1 * (pk_larvae$count > 0)
pk_larvae_gam1 <- gam(presence ~ 
                     s(doy, k = 4) +
                     s(lon, lat),
                   data = pk_larvae,
                   family = "binomial")

# gaussian
pk_larvae_gam2 <- gam(log(larvalcatchper10m2 + 1) ~ 
                     s(doy, k = 4) +
                     s(lon, lat) +
                     s(roms_salinity, k = 4) +
                     s(roms_temperature, k = 4),
                   data = pk_larvae[pk_larvae$larvalcatchper10m2 > 0, ])

# Prediction grid
nlat = 80
nlon = 120
latd = seq(min(pk_larvae$lat), max(pk_larvae$lat), length.out = nlat)
lond = seq(min(pk_larvae$lon), max(pk_larvae$lon), length.out = nlon)
grid_extent_larvaepk <- expand.grid(lond, latd)
names(grid_extent_larvaepk) <- c('lon', 'lat')

# Calculate distance of each grid point to closest 'positive observation'
grid_extent_larvaepk$dist <- NA
for (k in 1:nrow(grid_extent_larvaepk)) {
  dist <- distance_function(grid_extent_larvaepk$lat[k],
                            grid_extent_larvaepk$lon[k],
                            pk_larvae$lat,
                            pk_larvae$lon)
  grid_extent_larvaepk$dist[k] <- min(dist)
}

# Assign a within sample year and doy to the grid data
grid_extent_larvaepk$year <- 2022
grid_extent_larvaepk$date <- rep(as.Date("2022-02-18"), length(grid_extent_larvaepk))
grid_extent_larvaepk$doy <- rep(49, length(grid_extent_larvaepk))

# Attach ROMs forecast
grid_extent_larvaepk <- varid_match(grid_extent_larvaepk, temp_output, salt_output)

# Predict on forecasted output
grid_extent_larvaepk$pred1 <- predict(pk_larvae_gam1, newdata = grid_extent_larvaepk)
grid_extent_larvaepk$pred2 <- predict(pk_larvae_gam2, newdata = grid_extent_larvaepk)
grid_extent_larvaepk$pred <- grid_extent_larvaepk$pred1 * grid_extent_larvaepk$pred2
grid_extent_larvaepk$pred[grid_extent_larvaepk$dist > 30000] <- NA

# Plot
windows(width = 6, height = 5, family = "serif")
grid_predict(grid_extent_larvaepk, "Forecasted Distribution")

# Combined plot
windows(width = 12, height = 5, family = "serif")
par(mfrow = c(1, 2), 
    mai = c(0.9, 0.9, 0.5, 0.5))
grid_predict(grid_extent_eggpk, "Egg Distribution")
grid_predict(grid_extent_larvaepk, "Larval Distribution")
dev.copy(jpeg, here('results/pollock_forecast', 'pollock_forecast.jpg'), 
         height = 5, width = 12, res = 200, units = 'in', family = "serif")
dev.off()

#### Flathead Sole ----
##### Eggs ----
# Two-part binomial and gaussian
# binomial
fhs_egg$presence <- 1 * (fhs_egg$count > 0)
fhs_egg_gam1 <- gam(presence ~ 
                     s(doy, k = 4) +
                     s(lon, lat),
                   data = fhs_egg,
                   family = "binomial")


# gaussian
fhs_egg_gam2 <- gam(log(larvalcatchper10m2 + 1) ~ 
                     s(doy, k = 4) +
                     s(lon, lat) +
                     s(roms_salinity, k = 4) +
                     s(roms_temperature, k = 4),
                   data = fhs_egg[fhs_egg$larvalcatchper10m2 > 0, ])

# Prediction grid
nlat = 80
nlon = 120
latd = seq(min(fhs_egg$lat), max(fhs_egg$lat), length.out = nlat)
lond = seq(min(fhs_egg$lon), max(fhs_egg$lon), length.out = nlon)
grid_extent_eggfhs <- expand.grid(lond, latd)
names(grid_extent_eggfhs) <- c('lon', 'lat')

# Calculate distance of each grid point to closest 'positive observation'
grid_extent_eggfhs$dist <- NA
for (k in 1:nrow(grid_extent_eggfhs)) {
  dist <- distance_function(grid_extent_eggfhs$lat[k],
                            grid_extent_eggfhs$lon[k],
                            fhs_egg$lat,
                            fhs_egg$lon)
  grid_extent_eggfhs$dist[k] <- min(dist)
}

# Assign a within sample year and doy to the grid data
grid_extent_eggfhs$year <- 2022
grid_extent_eggfhs$date <- rep(as.Date("2022-02-18"), length(grid_extent_eggfhs))
grid_extent_eggfhs$doy <- rep(49, length(grid_extent_eggfhs))

# Attach ROMs forecast
grid_extent_eggfhs <- varid_match(grid_extent_eggfhs, temp_output, salt_output)

# Predict on forecasted output
grid_extent_eggfhs$pred1 <- predict(fhs_egg_gam1, newdata = grid_extent_eggfhs)
grid_extent_eggfhs$pred2 <- predict(fhs_egg_gam2, newdata = grid_extent_eggfhs)
grid_extent_eggfhs$pred <- grid_extent_eggfhs$pred1 * grid_extent_eggfhs$pred2
grid_extent_eggfhs$pred[grid_extent_eggfhs$dist > 30000] <- NA

# Plot
windows(width = 6, height = 5, family = "serif")
grid_predict(grid_extent_eggfhs, "Forecasted Distribution")

##### Larvae ----
# Two-part binomial and gaussian
# binomial
fhs_larvae$presence <- 1 * (fhs_larvae$count > 0)
fhs_larvae_gam1 <- gam(presence ~ 
                        s(doy, k = 4) +
                        s(lon, lat),
                      data = fhs_larvae,
                      family = "binomial")

# gaussian
fhs_larvae_gam2 <- gam(log(larvalcatchper10m2 + 1) ~ 
                        s(doy, k = 4) +
                        s(lon, lat) +
                        s(roms_salinity, k = 4) +
                        s(roms_temperature, k = 4),
                      data = fhs_larvae[fhs_larvae$larvalcatchper10m2 > 0, ])

# Prediction grid
nlat = 80
nlon = 120
latd = seq(min(fhs_larvae$lat), max(fhs_larvae$lat), length.out = nlat)
lond = seq(min(fhs_larvae$lon), max(fhs_larvae$lon), length.out = nlon)
grid_extent_larvaefhs <- expand.grid(lond, latd)
names(grid_extent_larvaefhs) <- c('lon', 'lat')

# Calculate distance of each grid point to closest 'positive observation'
grid_extent_larvaefhs$dist <- NA
for (k in 1:nrow(grid_extent_larvaefhs)) {
  dist <- distance_function(grid_extent_larvaefhs$lat[k],
                            grid_extent_larvaefhs$lon[k],
                            fhs_larvae$lat,
                            fhs_larvae$lon)
  grid_extent_larvaefhs$dist[k] <- min(dist)
}

# Assign a within sample year and doy to the grid data
grid_extent_larvaefhs$year <- 2022
grid_extent_larvaefhs$date <- rep(as.Date("2022-02-18"), length(grid_extent_larvaefhs))
grid_extent_larvaefhs$doy <- rep(49, length(grid_extent_larvaefhs))

# Attach ROMs forecast
grid_extent_larvaefhs <- varid_match(grid_extent_larvaefhs, temp_output, salt_output)

# Predict on forecasted output
grid_extent_larvaefhs$pred1 <- predict(fhs_larvae_gam1, newdata = grid_extent_larvaefhs)
grid_extent_larvaefhs$pred2 <- predict(fhs_larvae_gam2, newdata = grid_extent_larvaefhs)
grid_extent_larvaefhs$pred <- grid_extent_larvaefhs$pred1 * grid_extent_larvaefhs$pred2
grid_extent_larvaefhs$pred[grid_extent_larvaefhs$dist > 30000] <- NA

# Plot
windows(width = 6, height = 5, family = "serif")
grid_predict(grid_extent_larvaefhs, "Forecasted Distribution")

# Combined plot
windows(width = 12, height = 5, family = "serif")
par(mfrow = c(1, 2), 
    mai = c(0.9, 0.9, 0.5, 0.5))
grid_predict(grid_extent_eggfhs, "Egg Distribution")
grid_predict(grid_extent_larvaefhs, "Larval Distribution")
dev.copy(jpeg, here('results/flathead_forecast', 'flathead_forecast.jpg'), 
         width = 12, height = 5, res = 200, units = 'in', family = "serif")
dev.off()

#### Alaska Plaice ----
##### Eggs ----
# Two-part binomial and gaussian
# binomial
akp_egg$presence <- 1 * (akp_egg$count > 0)
akp_egg_gam1 <- gam(presence ~ 
                      s(doy, k = 4) +
                      s(lon, lat),
                    data = akp_egg,
                    family = "binomial")


# gaussian
akp_egg_gam2 <- gam(log(larvalcatchper10m2 + 1) ~ 
                      s(doy, k = 4) +
                      s(lon, lat) +
                      s(roms_salinity, k = 4) +
                      s(roms_temperature, k = 4),
                    data = akp_egg[akp_egg$larvalcatchper10m2 > 0, ])

# Prediction grid
nlat = 80
nlon = 120
latd = seq(min(akp_egg$lat), max(akp_egg$lat), length.out = nlat)
lond = seq(min(akp_egg$lon), max(akp_egg$lon), length.out = nlon)
grid_extent_eggakp <- expand.grid(lond, latd)
names(grid_extent_eggakp) <- c('lon', 'lat')

# Calculate distance of each grid point to closest 'positive observation'
grid_extent_eggakp$dist <- NA
for (k in 1:nrow(grid_extent_eggakp)) {
  dist <- distance_function(grid_extent_eggakp$lat[k],
                            grid_extent_eggakp$lon[k],
                            akp_egg$lat,
                            akp_egg$lon)
  grid_extent_eggakp$dist[k] <- min(dist)
}

# Assign a within sample year and doy to the grid data
grid_extent_eggakp$year <- 2022
grid_extent_eggakp$date <- rep(as.Date("2022-02-18"), length(grid_extent_eggakp))
grid_extent_eggakp$doy <- rep(49, length(grid_extent_eggakp))

# Attach ROMs forecast
grid_extent_eggakp <- varid_match(grid_extent_eggakp, temp_output, salt_output)

# Predict on forecasted output
grid_extent_eggakp$pred1 <- predict(akp_egg_gam1, newdata = grid_extent_eggakp)
grid_extent_eggakp$pred2 <- predict(akp_egg_gam2, newdata = grid_extent_eggakp)
grid_extent_eggakp$pred <- grid_extent_eggakp$pred1 * grid_extent_eggakp$pred2
grid_extent_eggakp$pred[grid_extent_eggakp$dist > 30000] <- NA

# Plot
windows(width = 6, height = 5, family = "serif")
grid_predict(grid_extent_eggakp, "Forecasted Distribution")

##### Larvae ----
# Two-part binomial and gaussian
# binomial
akp_larvae$presence <- 1 * (akp_larvae$count > 0)
akp_larvae_gam1 <- gam(presence ~ 
                         s(doy, k = 4) +
                         s(lon, lat),
                       data = akp_larvae,
                       family = "binomial")

# gaussian
akp_larvae_gam2 <- gam(log(larvalcatchper10m2 + 1) ~ 
                         s(doy, k = 4) +
                         s(lon, lat) +
                         s(roms_salinity, k = 4) +
                         s(roms_temperature, k = 4),
                       data = akp_larvae[akp_larvae$larvalcatchper10m2 > 0, ])

# Prediction grid
nlat = 80
nlon = 120
latd = seq(min(akp_larvae$lat), max(akp_larvae$lat), length.out = nlat)
lond = seq(min(akp_larvae$lon), max(akp_larvae$lon), length.out = nlon)
grid_extent_larvaeakp <- expand.grid(lond, latd)
names(grid_extent_larvaeakp) <- c('lon', 'lat')

# Calculate distance of each grid point to closest 'positive observation'
grid_extent_larvaeakp$dist <- NA
for (k in 1:nrow(grid_extent_larvaeakp)) {
  dist <- distance_function(grid_extent_larvaeakp$lat[k],
                            grid_extent_larvaeakp$lon[k],
                            akp_larvae$lat,
                            akp_larvae$lon)
  grid_extent_larvaeakp$dist[k] <- min(dist)
}

# Assign a within sample year and doy to the grid data
grid_extent_larvaeakp$year <- 2022
grid_extent_larvaeakp$date <- rep(as.Date("2022-02-18"), length(grid_extent_larvaeakp))
grid_extent_larvaeakp$doy <- rep(49, length(grid_extent_larvaeakp))

# Attach ROMs forecast
grid_extent_larvaeakp <- varid_match(grid_extent_larvaeakp, temp_output, salt_output)

# Predict on forecasted output
grid_extent_larvaeakp$pred1 <- predict(akp_larvae_gam1, newdata = grid_extent_larvaeakp)
grid_extent_larvaeakp$pred2 <- predict(akp_larvae_gam2, newdata = grid_extent_larvaeakp)
grid_extent_larvaeakp$pred <- grid_extent_larvaeakp$pred1 * grid_extent_larvaeakp$pred2
grid_extent_larvaeakp$pred[grid_extent_larvaeakp$dist > 30000] <- NA

# Plot
windows(width = 6, height = 5, family = "serif")
grid_predict(grid_extent_larvaeakp, "Forecasted Distribution")

# Combined plot
windows(width = 12, height = 5, family = "serif")
par(mfrow = c(1, 2), 
    mai = c(0.9, 0.9, 0.5, 0.5))
grid_predict(grid_extent_eggakp, "Egg Distribution")
grid_predict(grid_extent_larvaeakp, "Larval Distribution")
dev.copy(jpeg, here('results/plaice_forecast', 'plaice_forecast.jpg'), 
         width = 12, height = 5, res = 200, units = 'in', family = "serif")
dev.off()

#### Yellowfin Sole ----
##### Eggs ----
# Two-part binomial and gaussian
# binomial
yfs_egg$presence <- 1 * (yfs_egg$count > 0)
yfs_egg_gam1 <- gam(presence ~ 
                      s(doy, k = 4) +
                      s(lon, lat),
                    data = yfs_egg,
                    family = "binomial")


# gaussian
yfs_egg_gam2 <- gam(log(larvalcatchper10m2 + 1) ~ 
                      s(doy, k = 4) +
                      s(lon, lat) +
                      s(roms_salinity, k = 4) +
                      s(roms_temperature, k = 4),
                    data = yfs_egg[yfs_egg$larvalcatchper10m2 > 0, ])

# Prediction grid
nlat = 80
nlon = 120
latd = seq(min(yfs_egg$lat), max(yfs_egg$lat), length.out = nlat)
lond = seq(min(yfs_egg$lon), max(yfs_egg$lon), length.out = nlon)
grid_extent_eggyfs <- expand.grid(lond, latd)
names(grid_extent_eggyfs) <- c('lon', 'lat')

# Calculate distance of each grid point to closest 'positive observation'
grid_extent_eggyfs$dist <- NA
for (k in 1:nrow(grid_extent_eggyfs)) {
  dist <- distance_function(grid_extent_eggyfs$lat[k],
                            grid_extent_eggyfs$lon[k],
                            yfs_egg$lat,
                            yfs_egg$lon)
  grid_extent_eggyfs$dist[k] <- min(dist)
}

# Assign a within sample year and doy to the grid data
grid_extent_eggyfs$year <- 2022
grid_extent_eggyfs$date <- rep(as.Date("2022-02-18"), length(grid_extent_eggyfs))
grid_extent_eggyfs$doy <- rep(49, length(grid_extent_eggyfs))

# Attach ROMs forecast
grid_extent_eggyfs <- varid_match(grid_extent_eggyfs, temp_output, salt_output)

# Predict on forecasted output
grid_extent_eggyfs$pred1 <- predict(yfs_egg_gam1, newdata = grid_extent_eggyfs)
grid_extent_eggyfs$pred2 <- predict(yfs_egg_gam2, newdata = grid_extent_eggyfs)
grid_extent_eggyfs$pred <- grid_extent_eggyfs$pred1 * grid_extent_eggyfs$pred2
grid_extent_eggyfs$pred[grid_extent_eggyfs$dist > 30000] <- NA

# Plot
windows(width = 6, height = 5, family = "serif")
grid_predict(grid_extent_eggyfs, "Forecasted Distribution")

##### Larvae ----
# Two-part binomial and gaussian
# binomial
yfs_larvae$presence <- 1 * (yfs_larvae$count > 0)
yfs_larvae_gam1 <- gam(presence ~ 
                         s(doy, k = 4) +
                         s(lon, lat),
                       data = yfs_larvae,
                       family = "binomial")

# gaussian
yfs_larvae_gam2 <- gam(log(larvalcatchper10m2 + 1) ~ 
                         s(doy, k = 4) +
                         s(lon, lat) +
                         s(roms_salinity, k = 4) +
                         s(roms_temperature, k = 4),
                       data = yfs_larvae[yfs_larvae$larvalcatchper10m2 > 0, ])

# Prediction grid
nlat = 80
nlon = 120
latd = seq(min(yfs_larvae$lat), max(yfs_larvae$lat), length.out = nlat)
lond = seq(min(yfs_larvae$lon), max(yfs_larvae$lon), length.out = nlon)
grid_extent_larvaeyfs <- expand.grid(lond, latd)
names(grid_extent_larvaeyfs) <- c('lon', 'lat')

# Calculate distance of each grid point to closest 'positive observation'
grid_extent_larvaeyfs$dist <- NA
for (k in 1:nrow(grid_extent_larvaeyfs)) {
  dist <- distance_function(grid_extent_larvaeyfs$lat[k],
                            grid_extent_larvaeyfs$lon[k],
                            yfs_larvae$lat,
                            yfs_larvae$lon)
  grid_extent_larvaeyfs$dist[k] <- min(dist)
}

# Assign a within sample year and doy to the grid data
grid_extent_larvaeyfs$year <- 2022
grid_extent_larvaeyfs$date <- rep(as.Date("2022-02-18"), length(grid_extent_larvaeyfs))
grid_extent_larvaeyfs$doy <- rep(49, length(grid_extent_larvaeyfs))

# Attach ROMs forecast
grid_extent_larvaeyfs <- varid_match(grid_extent_larvaeyfs, temp_output, salt_output)

# Predict on forecasted output
grid_extent_larvaeyfs$pred1 <- predict(yfs_larvae_gam1, newdata = grid_extent_larvaeyfs)
grid_extent_larvaeyfs$pred2 <- predict(yfs_larvae_gam2, newdata = grid_extent_larvaeyfs)
grid_extent_larvaeyfs$pred <- grid_extent_larvaeyfs$pred1 * grid_extent_larvaeyfs$pred2
grid_extent_larvaeyfs$pred[grid_extent_larvaeyfs$dist > 30000] <- NA

# Plot
windows(width = 6, height = 5, family = "serif")
grid_predict(grid_extent_larvaeyfs, "Forecasted Distribution")

# Combined plot
windows(width = 12, height = 5, family = "serif")
par(mfrow = c(1, 2), 
    mai = c(0.9, 0.9, 0.5, 0.5))
grid_predict(grid_extent_eggyfs, "Egg Distribution")
grid_predict(grid_extent_larvaeyfs, "Larval Distribution")
dev.copy(jpeg, here('results/yellowfin_forecast', 'yellowfin_forecast.jpg'), 
         width = 12, height = 5, res = 200, units = 'in', family = "serif")
dev.off()