# Title: Import and cleaning of species data

# Load appropriate libraries and functions
library(mgcv)
library(maps)
library(mapdata)
library(marmap)
library(raster)
library(ncdf4)
library(spacetime)
library(fields)
library(date)
library(colorRamps) 
library(itsadug)   
library(RColorBrewer)
library(here)
library(tidyverse)
library(lubridate)
source(here('code/functions', 'distance_function.R'))
source(here('code/functions', 'vis_gam_COLORS.R'))

##### Depth data obtained from NGDC grid extract tool for ETOPO1 ----
str_name <- (here('data', 'bering_bathy.tiff'))
bering_bathy <- as.bathy(raster(str_name) * -1)
bathy_lat <- as.numeric(colnames(bering_bathy))
bathy_lon <- as.numeric(rownames(bering_bathy))
bathy_ylim = range(bathy_lat)
bathy_xlim = range(bathy_lon)
bering_bathy[bering_bathy <= -1] <- NA

# Bering 10K model output: avg surface temperatures
bering_model_temp1 <- nc_open(here('data/temperature_netcdf', 
                                  'B10K-K20_CORECFS_1985-1989_average_temp_surface5m.nc'))
bering_model_temp2 <- nc_open(here('data/temperature_netcdf', 
                                   'B10K-K20_CORECFS_1990-1994_average_temp_surface5m.nc'))
bering_model_temp3 <- nc_open(here('data/temperature_netcdf', 
                                   'B10K-K20_CORECFS_1995-1999_average_temp_surface5m.nc'))
bering_model_temp4 <- nc_open(here('data/temperature_netcdf', 
                                   'B10K-K20_CORECFS_2000-2004_average_temp_surface5m.nc'))
bering_model_temp5 <- nc_open(here('data/temperature_netcdf', 
                                   'B10K-K20_CORECFS_2005-2009_average_temp_surface5m.nc'))
bering_model_temp6 <- nc_open(here('data/temperature_netcdf', 
                                   'B10K-K20_CORECFS_2010-2014_average_temp_surface5m.nc'))
bering_model_temp7 <- nc_open(here('data/temperature_netcdf', 
                                   'B10K-K20_CORECFS_2015-2019_average_temp_surface5m.nc'))
bering_model_salt1 <- nc_open(here('data/salinity_netcdf', 
                                   'B10K-K20_CORECFS_1985-1989_average_salt_surface5m.nc'))
bering_model_salt2 <- nc_open(here('data/salinity_netcdf', 
                                   'B10K-K20_CORECFS_1990-1994_average_salt_surface5m.nc'))
bering_model_salt3 <- nc_open(here('data/salinity_netcdf', 
                                   'B10K-K20_CORECFS_1995-1999_average_salt_surface5m.nc'))
bering_model_salt4 <- nc_open(here('data/salinity_netcdf', 
                                   'B10K-K20_CORECFS_2000-2004_average_salt_surface5m.nc'))
bering_model_salt5 <- nc_open(here('data/salinity_netcdf', 
                                   'B10K-K20_CORECFS_2005-2009_average_salt_surface5m.nc'))
bering_model_salt6 <- nc_open(here('data/salinity_netcdf', 
                                   'B10K-K20_CORECFS_2010-2014_average_salt_surface5m.nc'))
bering_model_salt7 <- nc_open(here('data/salinity_netcdf', 
                                   'B10K-K20_CORECFS_2015-2019_average_salt_surface5m.nc'))
# check to see contents
print(bering_model_temp)

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
} # need to get time to return properly

# list of the four variables 
temp_output1 <- nc_extract(bering_model_temp1, temp, 'temp')
temp_output2 <- nc_extract(bering_model_temp2, temp, 'temp')
temp_output3 <- nc_extract(bering_model_temp3, temp, 'temp')
temp_output4 <- nc_extract(bering_model_temp4, temp, 'temp')
temp_output5 <- nc_extract(bering_model_temp5, temp, 'temp')
temp_output6 <- nc_extract(bering_model_temp6, temp, 'temp')
temp_output7 <- nc_extract(bering_model_temp7, temp, 'temp')

salt_output1 <- nc_extract(bering_model_sal1, salt, 'salt')
salt_output2 <- nc_extract(bering_model_sal2, salt, 'salt')
salt_output3 <- nc_extract(bering_model_sal3, salt, 'salt')
salt_output4 <- nc_extract(bering_model_sal4, salt, 'salt')
salt_output5 <- nc_extract(bering_model_sal5, salt, 'salt')
salt_output6 <- nc_extract(bering_model_sal6, salt, 'salt')
salt_output7 <- nc_extract(bering_model_sal7, salt, 'salt')

nc_close(bering_model_temp1)
nc_close(bering_model_temp2)
nc_close(bering_model_temp3)
nc_close(bering_model_temp4)
nc_close(bering_model_temp5)
nc_close(bering_model_temp6)
nc_close(bering_model_temp7)

nc_close(bering_model_salt1)
nc_close(bering_model_salt2)
nc_close(bering_model_salt3)
nc_close(bering_model_salt4)
nc_close(bering_model_salt5)
nc_close(bering_model_salt6)
nc_close(bering_model_salt7)

# Use if want to see dimensions of .nc variables
# dim(temp)
# range(temp, na.rm = T)
# dim(lat)
# range(lat, na.rm = T)
# dim(lon)
# range(lon, na.rm = T)
# dim(time1)
# range(time1, na.rm = T)
# nc_close(bering_model)

#### Egg and larval data ----
# From Steve email (8/02/19): The files are for 60 cm bongo catches and include both net 1 and 2 in case net 1 was a fail. For the case where catches for both net 1 and 2 are included at the same station, I would suggest using net 1 because that is the net typically used for quantitative catch.
ys_egg_raw <-
  read_csv(here('data/species_data', 'YFSole_Egg_Catch.csv'))
ys_larvae_raw <-
  read_csv(here('data/species_data', 'YFSole_Larvae_Catch.csv'))

# change to lowercase
names(ys_egg_raw) <- tolower(names(ys_egg_raw))
names(ys_larvae_raw) <- tolower(names(ys_larvae_raw))
colnames(ys_egg_raw)[23:25] <- c("month", "year", "day")
colnames(ys_larvae_raw)[23:25] <- c("month", "year", "day")

# See email from Steve about removing station with ID == 1SS02 81 1 60BON 2
ys_larvae_raw <- ys_larvae_raw[ys_larvae_raw$haul_id != '1SS02 81 1 60BON 2', ]
ys_egg_raw <- ys_egg_raw[ys_egg_raw$haul_id != '1SS02 81 1 60BON 2', ]


ys_larvae_raw$doy <- as.numeric(mdy.date(ys_larvae_raw$month, 
                                         ys_larvae_raw$day, 1960))
ys_egg_raw$doy <- as.numeric(mdy.date(ys_egg_raw$month, 
                                      ys_egg_raw$day, 1960))
ys_egg_raw$id <- paste(ys_egg_raw$cruise,
                       ys_egg_raw$lat,
                       ys_egg_raw$lon,
                       ys_egg_raw$gmt_date_time,
                       ys_egg_raw$mesh,
                       sep = "_")
ys_larvae_raw$id <- paste(ys_larvae_raw$cruise,
                          ys_larvae_raw$lat,
                          ys_larvae_raw$lon,
                          ys_larvae_raw$gmt_date_time,
                          ys_larvae_raw$mesh,
                          sep = "_")

# Check attributes of egg and larval data
windows()
par(mfrow = c(2, 4))
plot(table(ys_egg_raw$year[ys_egg_raw$larvalcatchper10m2 > 0]),
     ylab = 'Frequency',
     xlab = 'Year',
     main = 'Eggs')
plot(table(ys_egg_raw$month[ys_egg_raw$larvalcatchper10m2 > 0]),
     ylab = 'Frequency',
     xlab = 'Month',
     main = 'Eggs')
plot(table(ys_egg_raw$gear_name[ys_egg_raw$larvalcatchper10m2 > 0]),
     ylab = 'Frequency',
     xlab = 'Gear',
     main = 'Eggs')
hist(ys_egg_raw$lat[ys_egg_raw$larvalcatchper10m2 > 0],
     ylab = 'Frequency',
     xlab = 'Lat',
     main = 'Eggs')
plot(table(ys_larvae_raw$year[ys_larvae_raw$larvalcatchper10m2 > 0]),
     ylab = 'Frequency',
     xlab = 'Year',
     main = 'Larvae')
plot(table(ys_larvae_raw$month[ys_larvae_raw$larvalcatchper10m2 > 0]),
     ylab = 'Frequency',
     xlab = 'Month',
     main = 'Larvae')
plot(table(ys_larvae_raw$gear_name[ys_larvae_raw$larvalcatchper10m2 > 0]),
     ylab = 'Frequency',
     xlab = 'Gear',
     main = 'Larvae')
hist(ys_larvae_raw$lat[ys_larvae_raw$larvalcatchper10m2 > 0],
     ylab = 'Frequency',
     xlab = 'Lat',
     main = 'Larvae')

table(ys_egg_raw$haul_performance)
table(ys_larvae_raw$haul_performance)

# Egg data: check number of bongo nets per each station
table(ys_egg_raw$net) # Most sampled net is 1, but there are instances of net 2 (and 3!) being sampled.
ys_egg_raw[ys_egg_raw$net == 3, ] # All net 3 samples are from a 153 um mesh size. They need to be discarded.
ys_egg_1_2 <- ys_egg_raw[ys_egg_raw$net < 3, ] # Let us now check the net 2 samples
table(ys_egg_1_2$net) # There are 659 stations with net = 2 recorded
tmp <- table(ys_egg_1_2$id)
table(tmp) # Of these 84 have been sampled twice, which means that we only need to retain one of the two nets from these 84 stations
id_2 <- names(tmp)[tmp == 2] # ID of stations with 2 net records
ys_egg_1 <- ys_egg_1_2[ys_egg_1_2$net == 1, ] # These are data for records with net 1 only (n=5495)
dim(ys_egg_1)
ys_egg_2 <- ys_egg_1_2[ys_egg_1_2$id %in% id_2 &
               ys_egg_1_2$net == 1, ] # These are the data in which two nets were recorded, but I am only retaining net 1. We now need to remove data that were recorded with 2 nets, and then re-attach data with only net 1 (n=84)
dim(ys_egg_2)

ys_egg_1net <- ys_egg_1_2[!ys_egg_1_2$id %in% id_2, ] # These are data where only 1 net was recorded, either 1 or 2 (n=5495+660-84*2=5987)
dim(ys_egg_1net)
ys_egg <- rbind(ys_egg_1net, ys_egg_2)
dim(ys_egg) # Original n = 6182, of which n = 27 were net 3 (discarded), and n = 84*2 had two nets recorded. We should end up with n = 6181-27-84= 6070 data points

# Larval catch data: check number of bongo nets per each station
table(ys_larvae_raw$net) # Most sampled net is 1, but there are instances of net 2 (n=659) being sampled.
ys_larvae_raw[ys_larvae_raw$net == 3, ] # All net 3 samples are from a 153 um mesh size. They need to be discarded.
ys_larvae_1_2 <- ys_larvae_raw[ys_larvae_raw$net < 3, ] # Let us now check the net 2 samples
table(ys_larvae_1_2$net) # There are 659 stations with net = 2 recorded
tmp <- table(ys_larvae_1_2$id)
table(tmp) # Of these 84 have been sampled twice, which means that we only need to retain one of the two nets from these 84 stations
id_2 <- names(tmp)[tmp == 2] #ID of stations with 2 net records
ys_larvae_2 <- ys_larvae_1_2[ys_larvae_1_2$id %in% id_2 &
                  ys_larvae_1_2$net == 1, ] # These are the data in which two nets were recorded, but I am only retaining net 1 (n = 84).
dim(ys_larvae_2)
# We now need to remove data that were recorded with 2 nets (n=84), and then re-attach data with only net 1

ys_larvae_1net <- ys_larvae_1_2[!ys_larvae_1_2$id %in% id_2, ]#These are data where only 1 net was recorded, either 1 or 2 (n=net1+net2-net1&2=5495+660-84*2=5987)
dim(ys_larvae_1net)
ys_larvae <- rbind(ys_larvae_1net, ys_larvae_2)
dim(ys_larvae)#Original n = 6182, of which n = 27 were net 3 (discarded), and n = 84*2 had two nets recorded. We should end up with n = 6181-27-84= 6070 data points

# As per Lauren email (9/12/19) only retain station where primary_net == Y
table(ys_egg$primary_net)
table(ys_larvae$primary_net)
ys_egg <- ys_egg[ys_egg$primary_net == 'Y', ]
ys_larvae <- ys_larvae[ys_larvae$primary_net == 'Y', ]
table(ys_egg$primary_net)
table(ys_larvae$primary_net)

##### Trim egg and larval data -----
# Year: 1988 forwad
# Month: all
# Latitude: all
ys_egg <- ys_egg[ys_egg$year > 1987, ]
ys_larvae <- ys_larvae[ys_larvae$year > 1987, ]

# For each station, add distance to closest positive catch
ys_egg$dist <- NA
ys_larvae$dist <- NA
for (i in 1:nrow(ys_egg)) {
  ys_egg$dist[i] <- min(distance_function(ys_egg$lat[i], 
                                          ys_egg$lon[i], 
                                          ys_egg$lat[ys_egg$larvalcatchper10m2 > 0], 
                                          ys_egg$lon[ys_egg$larvalcatchper10m2 > 0]))
  ys_larvae$dist[i] <- min(distance_function(ys_larvae$lat[i],
                                             ys_larvae$lon[i],
                                             ys_larvae$lat[ys_larvae$larvalcatchper10m2 > 0],
                                             ys_larvae$lon[ys_larvae$larvalcatchper10m2 > 0]))
}
ys_egg <- mutate(ys_egg, date = make_date(year, month, day))
ys_larvae <- mutate(ys_larvae, date = make_date(year, month, day))

# Inspect new data
windows()
par(mfrow = c(2, 4))
plot(table(ys_egg$year[ys_egg$larvalcatchper10m2 > 0]),
     ylab = 'Frequency',
     xlab = 'Year',
     main = 'Eggs')
plot(table(ys_egg$month[ys_egg$larvalcatchper10m2 > 0]),
     ylab = 'Frequency',
     xlab = 'Month',
     main = 'Eggs')
plot(table(ys_egg$gear_name[ys_egg$larvalcatchper10m2 > 0]),
     ylab = 'Frequency',
     xlab = 'Gear',
     main = 'Eggs')
hist(ys_egg$lat[ys_egg$larvalcatchper10m2 > 0],
     ylab = 'Frequency',
     xlab = 'Lat',
     main = 'Eggs')
plot(table(ys_larvae$year[ys_larvae$larvalcatchper10m2 > 0]),
     ylab = 'Frequency',
     xlab = 'Year',
     main = 'Larvae')
plot(table(ys_larvae$month[ys_larvae$larvalcatchper10m2 > 0]),
     ylab = 'Frequency',
     xlab = 'Month',
     main = 'Larvae')
plot(table(ys_larvae$gear_name[ys_larvae$larvalcatchper10m2 > 0]),
     ylab = 'Frequency',
     xlab = 'Gear',
     main = 'Larvae')
hist(ys_larvae$lat[ys_larvae$larvalcatchper10m2 > 0],
     ylab = 'Frequency',
     xlab = 'Lat',
     main = 'Larvae')


# Select data for constraint analyses, including stations that are <30km away from the closest positive catch
subset_egg <- ys_egg[ys_egg$dist < 30000, 
                     c('larvalcatchper10m2', 'year', 'lat', 'lon', 'doy', 'date')]
subset_larvae <- ys_larvae[ys_larvae$dist < 30000, 
                           c('larvalcatchper10m2', 'year', 'lat', 'lon', 'doy', 'date')]
names(subset_egg) <- c('larvalcatchper10m2', 'year', 'lat', 'lon', 'doy', 'date')
names(subset_larvae) <- c('larvalcatchper10m2', 'year', 'lat', 'lon', 'doy', 'date')
table(subset_egg$year)
table(subset_larvae$year)

# Add Bering10K model temperatures and salinities
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

complete_egg <- varid_match(subset_egg, temp_output, salt_output)

# outside of function
subset_egg$roms_date <- NA
subset_egg$roms_temperature <- NA
subset_egg$roms_salinity <- NA
for (i in 1:nrow(subset_egg)) {
  idx_time <- order(abs(temp_output[[3]] - subset_egg$date[2]))[1]
  subset_egg$roms_date[2] <- temp_output[[3]][idx_time]
  idx_grid <- order(distance_function(
    subset_egg$lat[2],
    subset_egg$lon[2],
    c(temp_output[[2]]),
    c(temp_output[[1]])
  ))[1]
  subset_egg$roms_temperature[2] <- c(temp_output[[4]][, , idx_time])[idx_grid]
  subset_egg$roms_salinity[2] <- c(salt_output[[4]][, , idx_time])[idx_grid]}
