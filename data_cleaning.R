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

salt_output1 <- nc_extract(bering_model_salt1, salt, 'salt')
salt_output2 <- nc_extract(bering_model_salt2, salt, 'salt')
salt_output3 <- nc_extract(bering_model_salt3, salt, 'salt')
salt_output4 <- nc_extract(bering_model_salt4, salt, 'salt')
salt_output5 <- nc_extract(bering_model_salt5, salt, 'salt')
salt_output6 <- nc_extract(bering_model_salt6, salt, 'salt')
salt_output7 <- nc_extract(bering_model_salt7, salt, 'salt')

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
yfs_egg_raw <-
  read_csv(here('data/species_data', 'YFSole_Egg_Catch.csv'))
yfs_larvae_raw <-
  read_csv(here('data/species_data', 'YFSole_Larvae_Catch.csv'))
akp_egg_raw <-
  read_csv(here('data/species_data', 'BS_PlaiceEggCatch.csv'))
akp_larvae_raw <-
  read_csv(here('data/species_data', 'BS_PlaiceLarvaeCatch.csv'))
fhs_egg_raw <-
  read_csv(here('data/species_data', 'BS_FlatheadEggCatch.csv'))
fhs_larvae_raw <-
  read_csv(here('data/species_data', 'BS_FlatheadLarvaeCatch.csv'))
pk_egg_raw <-
  read_csv(here('data/species_data', 'BS_PollockEggCatch.csv'))
pk_larvae_raw <-
  read_csv(here('data/species_data', 'Bs_PollockLarvaeCatch.csv'))

### Functions ----
format_data <- function(data){
  names(data) <- tolower(names(data)) # change to lowercase
  colnames(data)[23:25] <- c("month", "year", "day")
  # See email from Steve about removing station with ID == 1SS02 81 1 60BON 2
  data <- data[data$haul_id != '1SS02 81 1 60BON 2', ]
  data$doy <- as.numeric(mdy.date(data$month, data$day, 1960))
  data$id <- paste(data$cruise,
                   data$lat,
                   data$lon,
                   data$gmt_date_time,
                   data$mesh,
                   sep = "_")
  return(data)
}
data_check <- function(egg, larvae){
  par(mfrow = c(2, 4))
  plot(table(egg$year[egg$larvalcatchper10m2 > 0]),
       ylab = 'Frequency',
       xlab = 'Year',
       main = 'Eggs')
  plot(table(egg$month[egg$larvalcatchper10m2 > 0]),
       ylab = 'Frequency',
       xlab = 'Month',
       main = 'Eggs')
  plot(table(egg$gear_name[egg$larvalcatchper10m2 > 0]),
       ylab = 'Frequency',
       xlab = 'Gear',
       main = 'Eggs')
  hist(egg$lat[egg$larvalcatchper10m2 > 0],
       ylab = 'Frequency',
       xlab = 'Lat',
       main = 'Eggs')
  plot(table(larvae$year[larvae$larvalcatchper10m2 > 0]),
       ylab = 'Frequency',
       xlab = 'Year',
       main = 'Larvae')
  plot(table(larvae$month[larvae$larvalcatchper10m2 > 0]),
       ylab = 'Frequency',
       xlab = 'Month',
       main = 'Larvae')
  plot(table(larvae$gear_name[larvae$larvalcatchper10m2 > 0]),
       ylab = 'Frequency',
       xlab = 'Gear',
       main = 'Larvae')
  hist(larvae$lat[larvae$larvalcatchper10m2 > 0],
       ylab = 'Frequency',
       xlab = 'Lat',
       main = 'Larvae')
}
clean_data <- function(data_formatted){
  # Most sampled net is 1, but there are instances of net 2 (and 3!) being sampled.
  data_formatted[data_formatted$net == 3, ] # All net 3 samples are from a 153 um mesh size. They need to be discarded.
  data_two <- data_formatted[data_formatted$net < 3, ] # Let us now check the net 2 samples
  # There are 659 stations with net = 2 recorded
  data_tmp <- table(data_two$id)
  # Of these 84 have been sampled twice, which means that we only need to retain one of the two nets from these 84 stations
  data_net2 <- names(data_tmp)[data_tmp == 2] # ID of stations with 2 net records
  data_net1 <- data_formatted[data_formatted$net == 1, ] # These are data for records with net 1 only (n=5495)
  data_both <- data_two[data_two$id %in% id_2 &
                                data_two$net == 1, ] # These are the data in which two nets were recorded, but I am only retaining net 1. We now need to remove data that were recorded with 2 nets, and then re-attach data with only net 1 (n=84)
  data_one <- data_two[!data_two$id %in% id_2, ] # These are data where only 1 net was recorded, either 1 or 2 (n=5495+660-84*2=5987)
  data_clean <- rbind(data_one, data_both)
  data_clean <- data_clean[data_clean$primary_net == 'Y', ]
  return(data_clean)
}
trim_data <- function(data_clean){
  data_clean <- data_clean[data_clean$year > 1987, ]
  # For each station, add distance to closest positive catch
  data_clean$dist <- NA
  for (i in 1:nrow(data_clean)) {
    data_clean$dist[i] <- min(distance_function(data_clean$lat[i], 
                                             data_clean$lon[i], 
                                             data_clean$lat[data_clean$larvalcatchper10m2 > 0], 
                                             data_clean$lon[data_clean$larvalcatchper10m2 > 0]))
  }
  data_final <- return(mutate(data_clean, date = make_date(year, month, day)))
}
final_data <- function(data_trim){
  data_subset <- data_trim[data_trim$dist < 30000, 
                            c('larvalcatchper10m2', 'year', 'lat', 'lon', 'doy', 'date')]
  names(data_subset) <- c('larvalcatchper10m2', 'year', 'lat', 'lon', 'doy', 'date')
  return(data_subset)
}

### Yellowfin Sole ----
# change to lowercase
yfs_egg_formatted <- format_data(yfs_egg_raw)
yfs_larvae_formatted <- format_data(yfs_larvae_raw)

# Check attributes of egg and larval data
data_check(yfs_egg_formatted, yfs_larvae_formatted)

table(yfs_egg_formatted$haul_performance)
table(yfs_larvae_formatted$haul_performance)

# Egg data: check number of bongo nets per each station
# These appear to be the same for all species, therefore use clean_data function for all other species
# table(yfs_egg_formatted$net) 
# yfs_egg_formatted[yfs_egg_formatted$net == 3, ] 
# yfs_egg_two <- yfs_egg_formatted[yfs_egg_formatted$net < 3, ] 
# yfs_egg_tmp <- table(yfs_egg_two$id)
# yfs_egg_net2 <- names(yfs_egg_tmp)[yfs_egg_tmp == 2] 
# yfs_egg_net1 <- yfs_egg_formatted[yfs_egg_formatted$net == 1, ] 
# dim(yfs_egg_net1)
# yfs_egg_both <- yfs_egg_two[yfs_egg_two$id %in% id_2 &
#                yfs_egg_two$net == 1, ] 
# dim(yfs_egg_both)
# yfs_egg_one <- yfs_egg_two[!yfs_egg_two$id %in% id_2, ] 
# dim(yfs_egg_one)
# yfs_egg <- rbind(yfs_egg_one, yfs_egg_both)
# dim(yfs_egg)
# 
# # Larval catch data: check number of bongo nets per each station
# # These appear to be the same for all species, use clean_data
# table(yfs_larvae_raw$net) 
# yfs_larvae_raw[yfs_larvae_raw$net == 3, ] 
# yfs_larvae_1_2 <- yfs_larvae_raw[yfs_larvae_raw$net < 3, ] 
# table(yfs_larvae_1_2$net) 
# tmp <- table(yfs_larvae_1_2$id)
# table(tmp) 
# id_2 <- names(tmp)[tmp == 2]
# yfs_larvae_2 <- yfs_larvae_1_2[yfs_larvae_1_2$id %in% id_2 &
#                   yfs_larvae_1_2$net == 1, ] 
# dim(yfs_larvae_2)
# yfs_larvae_1net <- yfs_larvae_1_2[!yfs_larvae_1_2$id %in% id_2, ]
# dim(yfs_larvae_1net)
# yfs_larvae <- rbind(yfs_larvae_1net, yfs_larvae_2)
# dim(yfs_larvae)

# As per Lauren email (9/12/19) only retain station where primary_net == Y
# table(yfs_egg$primary_net)
# table(yfs_larvae$primary_net)
# yfs_egg <- yfs_egg[yfs_egg$primary_net == 'Y', ]
# yfs_larvae <- yfs_larvae[yfs_larvae$primary_net == 'Y', ]

# Clean up data to remove unnecessary hauls
yfs_egg_clean <- clean_data(yfs_egg_formatted)
yfs_larvae_clean <- clean_data(yfs_larvae_formatted)

table(yfs_egg_clean$primary_net)
table(yfs_larvae_clean$primary_net)

##### Trim egg and larval data
# Year: 1988 forward
# Month: all
# Latitude: all
yfs_egg_trim <- trim_data(yfs_egg_clean)
yfs_larvae_trim <- trim_data(yfs_larvae_clean)

# Inspect new data
data_check(yfs_egg_trim, yfs_larvae_trim)

# Select data for constrained analyses, including stations that are <30km away from the closest positive catch
yfs_subset_egg <- final_data(yfs_egg_trim)
yfs_subset_larvae <- final_data(yfs_larvae_trim)
table(yfs_subset_egg$year)
table(yfs_subset_larvae$year)

#### Split egg and larvae data into separate datasets by date
yfs_subset_egg1 <- filter(yfs_subset_egg, year < 1990)
yfs_complete_larvae1 <- filter(yfs_subset_egg, year < 1990)
yfs_subset_egg2 <- filter(yfs_subset_egg, year > 1989 & year < 1995)
yfs_complete_larvae2 <- filter(yfs_subset_larvae, year > 1989 & year < 1995)
yfs_subset_egg3 <- filter(yfs_subset_egg, year > 1994 & year < 2000)
yfs_complete_larvae3 <- filter(yfs_subset_larvae, year > 1994 & year < 2000)
yfs_subset_egg4 <- filter(yfs_subset_egg, year > 1999 & year < 2005)
yfs_complete_larvae4 <- filter(yfs_subset_larvae, year > 1999 & year < 2005)
yfs_subset_egg5 <- filter(yfs_subset_egg, year > 2004 & year < 2010)
yfs_complete_larvae5 <- filter(yfs_subset_larvae, year > 2004 & year < 2010)
yfs_subset_egg6 <- filter(yfs_subset_egg, year > 2009 & year < 2015)
yfs_complete_larvae6 <- filter(yfs_subset_larvae, year > 2009 & year < 2015)
yfs_subset_egg7 <- filter(yfs_subset_egg, year > 2014 & year < 2020)
yfs_complete_larvae7 <- filter(yfs_subset_larvae, year > 2014 & year < 2020)

#### Add Bering10K model temperatures and salinities
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

yfs_complete_egg1 <- varid_match(yfs_subset_egg1, temp_output1, salt_output1)
yfs_complete_egg2 <- varid_match(yfs_subset_egg2, temp_output2, salt_output2)
yfs_complete_egg3 <- varid_match(yfs_subset_egg3, temp_output3, salt_output3)
yfs_complete_egg4 <- varid_match(yfs_subset_egg4, temp_output4, salt_output4)
yfs_complete_egg5 <- varid_match(yfs_subset_egg5, temp_output5, salt_output5)
yfs_complete_egg6 <- varid_match(yfs_subset_egg6, temp_output6, salt_output6)
yfs_complete_egg7 <- varid_match(yfs_subset_egg7, temp_output7, salt_output7)

yfs_complete_larvae1 <- varid_match(yfs_complete_larvae1, temp_output1, salt_output1)
yfs_complete_larvae2 <- varid_match(yfs_complete_larvae2, temp_output2, salt_output2)
yfs_complete_larvae3 <- varid_match(yfs_complete_larvae3, temp_output3, salt_output3)
yfs_complete_larvae4 <- varid_match(yfs_complete_larvae4, temp_output4, salt_output4)
yfs_complete_larvae5 <- varid_match(yfs_complete_larvae5, temp_output5, salt_output5)
yfs_complete_larvae6 <- varid_match(yfs_complete_larvae6, temp_output6, salt_output6)
yfs_complete_larvae7 <- varid_match(yfs_complete_larvae7, temp_output7, salt_output7)

# combine to create whole datasets for each
yfs_complete_egg <- rbind(yfs_complete_egg1, yfs_complete_egg2, yfs_complete_egg3,
                      yfs_complete_egg4, yfs_complete_egg5, yfs_complete_egg6,
                      yfs_complete_egg7)
yfs_complete_larvae <- rbind(yfs_complete_larvae1, yfs_complete_larvae2, yfs_complete_larvae3,
                      yfs_complete_larvae4, yfs_complete_larvae5, yfs_complete_larvae6,
                      yfs_complete_larvae7)

