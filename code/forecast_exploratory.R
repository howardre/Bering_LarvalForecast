# Title: Forecasting
# Date: complete 7/07/2021

### Libraries, functions, and data ----
library(maps)
library(maptools)
library(mapdata)
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

# Fish data
# Load ROMS temperature means
roms_temps <- readRDS(here('data', 'roms_temps.rds'))

# Load fish data
pk_egg <- as.data.frame(filter((readRDS(here('data', 'pk_egg.rds'))),
                               lat >= 52 & lat <= 62,
                               lon >= -176.5 & lon <= -156.5))
pk_egg$mean_temp <- roms_temps$mean[match(pk_egg$year, roms_temps$year)]

pk_larvae <- as.data.frame(filter(readRDS(here('data', 'pk_larvae.rds')),
                                  lat >= 52 & lat <= 62,
                                  lon >= -176.5 & lon <= -156.5))
pk_larvae$mean_temp <- roms_temps$mean[match(pk_larvae$year, roms_temps$year)]

fhs_egg <- as.data.frame(filter(readRDS(here('data', 'fhs_egg.rds')),
                                lat >= 52 & lat <= 62,
                                lon >= -176.5 & lon <= -156.5))
fhs_egg$mean_temp <- roms_temps$mean[match(fhs_egg$year, roms_temps$year)]

fhs_larvae <- as.data.frame(filter(readRDS(here('data', 'fhs_larvae.rds')),
                                   lat >= 52 & lat <= 62,
                                   lon >= -176.5 & lon <= -156.5))
fhs_larvae$mean_temp <- roms_temps$mean[match(fhs_larvae$year, roms_temps$year)]

yfs_egg <- as.data.frame(filter(readRDS(here('data', 'yfs_egg.rds')),
                                lat >= 52 & lat <= 62,
                                lon >= -176.5 & lon <= -156.5))
yfs_egg$mean_temp <- roms_temps$mean[match(yfs_egg$year, roms_temps$year)]

yfs_larvae <- as.data.frame(filter(readRDS(here('data', 'yfs_larvae.rds')),
                                   lat >= 52 & lat <= 62,
                                   lon >= -176.5 & lon <= -156.5))
yfs_larvae$mean_temp <- roms_temps$mean[match(yfs_larvae$year, roms_temps$year)]

akp_egg <- as.data.frame(filter(readRDS(here('data', 'akp_egg.rds')),
                                lat >= 52 & lat <= 62,
                                lon >= -176.5 & lon <= -156.5))
akp_egg$mean_temp <- roms_temps$mean[match(akp_egg$year, roms_temps$year)]

akp_larvae <- as.data.frame(filter(readRDS(here('data', 'akp_larvae.rds')),
                                   lat >= 52 & lat <= 62,
                                   lon >= -176.5 & lon <= -156.5))
akp_larvae$mean_temp <- roms_temps$mean[match(akp_larvae$year, roms_temps$year)]


# Match ROMS output function
varid_match <- function(data, model_output1, model_output2){
  data$roms_date <- NA
  data$roms_temperature <- NA
  data$roms_salinity <- NA
  for (i in 1:nrow(data)) {
    idx_time <- order(abs(model_output1[[3]] - data$date[i]))[1]
    data$roms_date[i] <- model_output1[[3]][idx_time]
    idx_grid <- order(distance_function(data$lat[i],
                                        data$lon[i],
                                        c(model_output1[[2]]),
                                        c(model_output1[[1]])))[1]
    data$roms_temperature[i] <- c(model_output1[[4]][, , , idx_time])[idx_grid]
    data$roms_salinity[i] <- c(model_output2[[4]][, , , idx_time])[idx_grid]
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

# Function to extract netcdf
years_list <- list('2015-2019', '2020-2024', '2025-2029', '2030-2034',
                '2035-2039', '2040-2044', '2045-2049', '2050-2054',
                '2055-2059', '2060-2064', '2065-2069', '2070-2074',
                '2075-2079', '2080-2084', '2085-2089', '2090-2094',
                '2095-2099')

get_temp_filepath <- function(years){
  filepath = paste('D:/B10K-K20P19_CMIP6_gfdl_ssp126/Level1/B10K-K20P19_CMIP6_gfdl_ssp126_', years, '_average_temp.nc', sep = '')
  return(filepath)
}

get_salt_filepath <- function(years){
  filepath = paste('D:/B10K-K20P19_CMIP6_gfdl_ssp126/Level1/B10K-K20P19_CMIP6_gfdl_ssp126_', years, '_average_salt.nc', sep = '')
  return(filepath)
}

# Function to get predictions
get_preds <- function(y, years_list, data, year, 
                      date, doy, start_date, end_date){
  # Prediction grid
  nlat = 80
  nlon = 120
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
  
  # Extract from netcdf
  bering_model_temp = nc_open(get_temp_filepath(years_list[y]))
  bering_model_salt = nc_open(get_salt_filepath(years_list[y]))
  
  temp_output <- nc_extract(bering_model_temp, temp, 'temp')
  salt_output <- nc_extract(bering_model_salt, salt, 'salt')
  
  nc_close(bering_model_temp)
  nc_close(bering_model_salt)
  
  # Assign a within sample year and doy to the grid data
  grid_extent$year <- year
  grid_extent$date <- rep(as.Date(date),
                          length(grid_extent))
  grid_extent$doy <- rep(doy, length(grid_extent))
  
  # Attach ROMs forecast
  grid_extent <- varid_match(grid_extent, temp_output, salt_output)
  
  # Calculate mean temperature
  time_index <- temp_output[[3]] >= start_date & temp_output[[3]] <= end_date
  temp_array <- temp_output[[4]][, , , time_index]
  
  # Select out box on the shelf
  temp_data <- as.data.frame(cbind(lon = as.vector(temp_output[[1]]), 
                                     lat = as.vector(temp_output[[2]]), 
                                     temp = as.vector(temp_array)))
  temp_filtered <- temp_data %>% filter(lon >= -170 & lon <= -165, lat >= 56 & lat <= 58)
  mean <- mean(temp_filtered$temp, na.rm = T)
  
  grid_extent$mean_temp <- mean
  
  gam <- gam(larvalcatchper10m2 + 1 ~ s(year) +
               s(doy, k = 8) +
               s(lon, lat) +
               s(roms_temperature, k = 6) +
               s(roms_salinity, k = 6) +
               s(lat, lon, by = mean_temp, k = 6),
             data = data,
             family = tw(link = 'log'),
             method = 'REML')
  
  # Predict on forecasted output
  grid_extent$pred <- predict(gam, newdata = grid_extent)
  grid_extent$pred[grid_extent$dist > 30000] <- NA
  
  return(grid_extent)
  
}



### Predict future distributions ----
#### Pollock ----
##### Eggs ----
grids_pkegg1 <- list()
for(k in 2015:2019){
grid <- get_preds(1, years_list, pk_egg, k, "2015-05-10", 
                  130, "2015-02-01", "2015-04-30")
grids_pkegg1[[paste("year", k)]] <- grid
} # only produces one year, probably not enough memory to finish



# Plot
windows(width = 6, height = 5, family = "serif")
grid_predict(grid_pkegg1, "Forecasted Distribution")

