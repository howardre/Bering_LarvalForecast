### Libraries, functions, and data ----
library(maps)
library(maptools)
library(mapdata)
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
# Load ROMS temperature means and forecast
roms_temps <- readRDS(here('data', 'roms_temps.rds'))
temps_cesm_ssp126 <- readRDS(here('data', 'temps_cesm_ssp126.rds'))
salts_cesm_ssp126 <- readRDS(here('data', 'salts_cesm_ssp126.rds'))

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
    idx_time <- order(abs(model_output1[[1]][[3]] - data$date[i]))[1]
    data$roms_date[i] <- model_output1[[1]][[3]][idx_time]
    idx_grid <- order(distance_function(data$lat[i],
                                        data$lon[i],
                                        c(model_output1[[1]][[2]]),
                                        c(model_output1[[1]][[1]])))[1]
    data$roms_temperature[i] <- c(model_output1[[1]][[4]][, , idx_time])[idx_grid]
    data$roms_salinity[i] <- c(model_output2[[1]][[4]][, , idx_time])[idx_grid]
  }
  return(data)
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


# Function to get predictions
get_preds <- function(data, year, date, doy, 
                      start_date, end_date,
                      temp_output, salt_output){
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
  
  # Assign a within sample year and doy to the grid data
  grid_extent$year <- year
  grid_extent$date <- rep(as.Date(date),
                          length(grid_extent))
  grid_extent$doy <- rep(doy, length(grid_extent))
  
  # Attach ROMs forecast
  grid_extent <- varid_match(grid_extent, temp_output, salt_output)
  
  # Calculate mean temperature
  time_index <- temp_output[[1]][[3]] >= start_date & temp_output[[1]][[3]] <= end_date
  temp_array <- temp_output[[1]][[4]][, , time_index]
  
  # Select out box on the shelf
  temp_data <- as.data.frame(cbind(lon = as.vector(temp_output[[1]][[1]]), 
                                   lat = as.vector(temp_output[[1]][[2]]), 
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



# Troubleshoot function
nlat = 80
nlon = 120
latd = seq(min(pk_egg$lat), max(pk_egg$lat), length.out = nlat)
lond = seq(min(pk_egg$lon), max(pk_egg$lon), length.out = nlon)
grid_extent <- expand.grid(lond, latd)
names(grid_extent) <- c('lon', 'lat')

# Calculate distance of each grid point to closest 'positive observation'
grid_extent$dist <- NA
for (k in 1:nrow(grid_extent)) {
  dist <- distance_function(grid_extent$lat[k],
                            grid_extent$lon[k],
                            pk_egg$lat,
                            pk_egg$lon)
  grid_extent$dist[k] <- min(dist)
}

# Assign a within sample year and doy to the grid data
grid_extent$year <- 2015
grid_extent$date <- rep(as.Date("2015-05-10"),
                        length(grid_extent))
grid_extent$doy <- rep(130, length(grid_extent))

# Attach ROMs forecast
grid_extent <- varid_match(grid_extent, temps_cesm_ssp126, salts_cesm_ssp126)

# Calculate mean temperature
time_index <- temps_cesm_ssp126[[1]][[3]] >= "2015-03-01" & temps_cesm_ssp126[[1]][[3]] <= "2015-04-30"
temp_array <- temps_cesm_ssp126[[1]][[4]][, , time_index]

# Select out box on the shelf
temp_pk_egg <- as.pk_egg.frame(cbind(lon = as.vector(temp_output[[1]][[1]]), 
                                 lat = as.vector(temp_output[[1]][[2]]), 
                                 temp = as.vector(temp_array)))
temp_filtered <- temp_pk_egg %>% filter(lon >= -170 & lon <= -165, lat >= 56 & lat <= 58)
mean <- mean(temp_filtered$temp, na.rm = T)

grid_extent$mean_temp <- mean

gam <- gam(larvalcatchper10m2 + 1 ~ s(year) +
             s(doy, k = 8) +
             s(lon, lat) +
             s(roms_temperature, k = 6) +
             s(roms_salinity, k = 6) +
             s(lat, lon, by = mean_temp, k = 6),
           pk_egg = pk_egg,
           family = tw(link = 'log'),
           method = 'REML')

# Predict on forecasted output
grid_extent$pred <- predict(gam, newpk_egg = grid_extent)
grid_extent$pred[grid_extent$dist > 30000] <- NA




grids_pkegg2 <- list()
for(j in 2015:2019){
  date1 <- paste(j, "-05-10", sep = "")
  date2 <- paste(j, "-03-01", sep = "")
  date3 <- paste(j, "-04-30", sep = "")
  grid <- get_preds(pk_egg, j, date1, 130, 
                    date2, date3, temps_cesm_ssp126,
                    salts_cesm_ssp126)
  grids_pkegg2[[paste("year", j, sep = "")]] <- grid
}
