# Title: Import and cleaning of species data
# Date: complete 5/25/2021

### Libraries and functions ----
library(marmap)
library(raster)
library(ncdf4)
library(spacetime)
library(fields)
library(here)
library(tidyverse)
library(lubridate)
library(date)
source(here('code/functions', 'distance_function.R'))

### Bering 10K model output ----
# Using avg surface temperatures & salinity
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

# open cesm file
cesm_dfs <- readRDS('D:/data/cesm_dfs_trim.rds')

# check to see contents
print(bering_model_temp1)

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

# close netcdf's
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

# Load historical
temps_cesm_historical <- readRDS(here('data', 'temps_cesm_historical.rds'))
salts_cesm_historical <- readRDS(here('data', 'salts_cesm_historical.rds'))
temps_gfdl_historical <- readRDS(here('data', 'temps_gfdl_historical.rds'))
salts_gfdl_historical <- readRDS(here('data', 'salts_gfdl_historical.rds'))
temps_miroc_historical <- readRDS(here('data', 'temps_miroc_historical.rds'))
salts_miroc_historical <- readRDS(here('data', 'salts_miroc_historical.rds'))

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

### Egg and larval data import ----
# From Steve email (8/02/19): The files are for 60 cm bongo catches and include both net 1 and 2 
# in case net 1 was a fail. For the case where catches for both net 1 and 2 are included at the 
# same station, I would suggest using net 1 because that is the net typically used for quantitative catch.
yfs_egg_raw <- read_csv(here('data/species_data', 'YFSole_Egg_Catch.csv'))
yfs_larvae_raw <- read_csv(here('data/species_data', 'YFSole_Larvae_Catch.csv'))
akp_egg_raw <- read_csv(here('data/species_data', 'BS_PlaiceEggCatch.csv'))
akp_larvae_raw <- read_csv(here('data/species_data', 'BS_PlaiceLarvaeCatch.csv'))
fhs_egg_raw <- read_csv(here('data/species_data', 'BS_FlatheadEggCatch.csv'))
fhs_larvae_raw <- read_csv(here('data/species_data', 'BS_FlatheadLarvaeCatch.csv'))
pk_egg_raw <- read_csv(here('data/species_data', 'BS_PollockEggCatch.csv'))
pk_larvae_raw <- read_csv(here('data/species_data', 'Bs_PollockLarvaeCatch.csv'))
pcod_larvae_raw <- read_csv(here('data/species_data', 'Bs_PacificCodLarvaeCatch.csv'))
nrs_larvae_raw <- read_csv(here('data/species_data', 'Bs_NorthernRockSoleLarvaeCatch.csv'))

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
  # All net 3 samples are from a 153 um mesh size. They need to be discarded.
  data_two <- data_formatted[data_formatted$net < 3, ] # Let us now check the net 2 samples
  # There are 659 stations with net = 2 recorded
  data_tmp <- table(data_two$id)
  # Of these 84 have been sampled twice, which means that we only need to retain one of the two nets from these 84 stations
  data_net2 <- names(data_tmp)[data_tmp == 2] # ID of stations with 2 net records
  data_net1 <- data_formatted[data_formatted$net == 1, ] # These are data for records with net 1 only (n=5495)
  data_both <- data_two[data_two$id %in% data_net2 &
                          data_two$net == 1, ] # These are the data in which two nets were recorded, but I am only retaining net 1. We now need to remove data that were recorded with 2 nets, and then re-attach data with only net 1 (n=84)
  data_one <- data_two[!data_two$id %in% data_net2, ] # These are data where only 1 net was recorded, either 1 or 2 (n=5495+660-84*2=5987)
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
                           c('larvalcatchper1000m3', 'larvalcatchper10m2', 
                             'year', 'lat', 'lon', 'doy', 'date',
                             'volume_filtered')]
  names(data_subset) <- c('larvalcatchper1000m3', 'larvalcatchper10m2',
                          'year', 'lat', 'lon', 'doy', 'date',
                          'volume_filtered')
  return(data_subset)
}
varid_match_proj <- function(data, model_output1, model_output2, 
                             model_output3, model_output4, model_output5,
                             model_output6, model_output7, model_output8, list){
  data$roms_date <- NA
  data$roms_temperature <- NA
  data$roms_salinity <- NA
  data$roms_temp_hist_cesm <- NA
  data$roms_salt_hist_cesm <- NA
  data$roms_temp_hist_gfdl <- NA
  data$roms_salt_hist_gfdl <- NA
  data$roms_temp_hist_miroc <- NA
  data$roms_salt_hist_miroc <- NA
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
    data$roms_temp_hist_cesm[i] <- c(model_output3[[list]][[4]][, , idx_time])[idx_grid]
    data$roms_salt_hist_cesm[i] <- c(model_output4[[list]][[4]][, , idx_time])[idx_grid]
    data$roms_temp_hist_gfdl[i] <- c(model_output5[[list]][[4]][, , idx_time])[idx_grid]
    data$roms_salt_hist_gfdl[i] <- c(model_output6[[list]][[4]][, , idx_time])[idx_grid]
    data$roms_temp_hist_miroc[i] <- c(model_output7[[list]][[4]][, , idx_time])[idx_grid]
    data$roms_salt_hist_miroc[i] <- c(model_output8[[list]][[4]][, , idx_time])[idx_grid]
  }
  return(data)
} 

varid_match <- function(data, model_output1, model_output2, list){
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

### Yellowfin Sole ----
# change to lowercase
yfs_egg_formatted <- format_data(yfs_egg_raw)
yfs_larvae_formatted <- format_data(yfs_larvae_raw)

# Check attributes of egg and larval data
data_check(yfs_egg_formatted, yfs_larvae_formatted)

table(yfs_egg_formatted$haul_performance)
table(yfs_larvae_formatted$haul_performance)

# Clean up data to remove unnecessary hauls
yfs_egg_clean <- clean_data(yfs_egg_formatted)
yfs_larvae_clean <- clean_data(yfs_larvae_formatted)

table(yfs_egg_clean$primary_net)
table(yfs_larvae_clean$primary_net)

# Trim egg and larval data
# Year: 1988 forward, remove 1997 and 2011 for sporadic sampling
# Month: remove May & June for larvae
# Latitude: all
yfs_egg_trim <- trim_data(yfs_egg_clean)
yfs_egg_trim <- filter(yfs_egg_trim,
                       month > 7,
                       year != 1997, year != 2011)
yfs_larvae_trim <- trim_data(yfs_larvae_clean)
yfs_larvae_trim <- filter(yfs_larvae_trim,
                          month > 7,
                          year != 1997, year != 2011)

# Inspect new data
data_check(yfs_egg_trim, yfs_larvae_trim)

# Select data for constrained analyses, including stations that are <30km away from the closest positive catch
yfs_subset_egg <- final_data(yfs_egg_trim)
yfs_subset_larvae <- final_data(yfs_larvae_trim)
table(yfs_subset_egg$year)
table(yfs_subset_larvae$year)

# Split egg and larvae data into separate datasets by date
yfs_subset_egg1 <- filter(yfs_subset_egg, year < 1990)
yfs_complete_larvae1 <- filter(yfs_subset_larvae, year < 1990)
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

# Add Bering10K model temperatures and salinities
yfs_complete_egg1 <- varid_match(yfs_subset_egg1, temp_output1, salt_output1, 
                                 temps_cesm_historical, salts_cesm_historical,
                                 temps_gfdl_historical, salts_gfdl_historical,
                                 temps_miroc_historical, salts_miroc_historical, 1)
yfs_complete_egg2 <- varid_match(yfs_subset_egg2, temp_output2, salt_output2, 
                                 temps_cesm_historical, salts_cesm_historical,
                                 temps_gfdl_historical, salts_gfdl_historical,
                                 temps_miroc_historical, salts_miroc_historical, 2)
yfs_complete_egg3 <- varid_match(yfs_subset_egg3, temp_output3, salt_output3, 
                                 temps_cesm_historical, salts_cesm_historical,
                                 temps_gfdl_historical, salts_gfdl_historical,
                                 temps_miroc_historical, salts_miroc_historical, 3)
yfs_complete_egg4 <- varid_match(yfs_subset_egg4, temp_output4, salt_output4, 
                                 temps_cesm_historical, salts_cesm_historical,
                                 temps_gfdl_historical, salts_gfdl_historical,
                                 temps_miroc_historical, salts_miroc_historical, 4)
yfs_complete_egg5 <- varid_match(yfs_subset_egg5, temp_output5, salt_output5, 
                                 temps_cesm_historical, salts_cesm_historical,
                                 temps_gfdl_historical, salts_gfdl_historical,
                                 temps_miroc_historical, salts_miroc_historical, 5)
yfs_complete_egg6 <- varid_match(yfs_subset_egg6, temp_output6, salt_output6, 
                                 temps_cesm_historical, salts_cesm_historical,
                                 temps_gfdl_historical, salts_gfdl_historical,
                                 temps_miroc_historical, salts_miroc_historical, 6)
yfs_complete_egg7 <- varid_match(yfs_subset_egg7, temp_output7, salt_output7, 
                                 temps_cesm_historical, salts_cesm_historical,
                                 temps_gfdl_historical, salts_gfdl_historical,
                                 temps_miroc_historical, salts_miroc_historical, 7)

yfs_complete_larvae1 <- varid_match(yfs_complete_larvae1, temp_output1, salt_output1, 
                                    temps_cesm_historical, salts_cesm_historical,
                                    temps_gfdl_historical, salts_gfdl_historical,
                                    temps_miroc_historical, salts_miroc_historical, 1)
yfs_complete_larvae2 <- varid_match(yfs_complete_larvae2, temp_output2, salt_output2, 
                                    temps_cesm_historical, salts_cesm_historical,
                                    temps_gfdl_historical, salts_gfdl_historical,
                                    temps_miroc_historical, salts_miroc_historical, 2)
yfs_complete_larvae3 <- varid_match(yfs_complete_larvae3, temp_output3, salt_output3, 
                                    temps_cesm_historical, salts_cesm_historical,
                                    temps_gfdl_historical, salts_gfdl_historical,
                                    temps_miroc_historical, salts_miroc_historical, 3)
yfs_complete_larvae4 <- varid_match(yfs_complete_larvae4, temp_output4, salt_output4, 
                                    temps_cesm_historical, salts_cesm_historical,
                                    temps_gfdl_historical, salts_gfdl_historical,
                                    temps_miroc_historical, salts_miroc_historical, 4)
yfs_complete_larvae5 <- varid_match(yfs_complete_larvae5, temp_output5, salt_output5, 
                                    temps_cesm_historical, salts_cesm_historical,
                                    temps_gfdl_historical, salts_gfdl_historical,
                                    temps_miroc_historical, salts_miroc_historical, 5)
yfs_complete_larvae6 <- varid_match(yfs_complete_larvae6, temp_output6, salt_output6, 
                                    temps_cesm_historical, salts_cesm_historical,
                                    temps_gfdl_historical, salts_gfdl_historical,
                                    temps_miroc_historical, salts_miroc_historical, 6)
yfs_complete_larvae7 <- varid_match(yfs_complete_larvae7, temp_output7, salt_output7, 
                                    temps_cesm_historical, salts_cesm_historical,
                                    temps_gfdl_historical, salts_gfdl_historical,
                                    temps_miroc_historical, salts_miroc_historical, 7)

# combine to create whole datasets for each
yfs_complete_egg <- rbind(yfs_complete_egg1, yfs_complete_egg2, yfs_complete_egg3,
                          yfs_complete_egg4, yfs_complete_egg5, yfs_complete_egg6,
                          yfs_complete_egg7)
yfs_complete_larvae <- rbind(yfs_complete_larvae1, yfs_complete_larvae2, yfs_complete_larvae3,
                             yfs_complete_larvae4, yfs_complete_larvae5, yfs_complete_larvae6,
                             yfs_complete_larvae7)

# add counts
yfs_complete_egg$count <- round(yfs_complete_egg$larvalcatchper1000m3 * 
                                  yfs_complete_egg$volume_filtered / 1000, 0)
yfs_complete_larvae$count <- round(yfs_complete_larvae$larvalcatchper1000m3 * 
                                     yfs_complete_larvae$volume_filtered / 1000, 0)

## bias correction
# calculate means by week for hindcast and historical
yfs_complete_egg$week <- week(yfs_complete_egg$date)
yfs_egg_data_baseline <- yfs_complete_egg %>%
  group_by(week, lat, lon) %>%
  summarise(week_baseline_temp = mean(roms_temperature),
            week_baseline_salt = mean(roms_salinity),
            week_cesm_temp = mean(roms_temp_hist_cesm),
            week_cesm_salt = mean(roms_salt_hist_cesm),
            week_gfdl_temp = mean(roms_temp_hist_gfdl),
            week_gfdl_salt = mean(roms_salt_hist_gfdl),
            week_miroc_temp = mean(roms_temp_hist_miroc),
            week_miroc_salt = mean(roms_salt_hist_miroc))

# Match back to the main data frame


# import the forecast, summarize means by year, lat, lon, month
temps_cesm_ssp126 <- readRDS(here('data', 'temps_cesm_ssp126.rds'))
salts_cesm_ssp126 <- readRDS(here('data', 'salts_cesm_ssp126.rds'))

# create dataframes of hindcast, historical, and forecast
temp_vec <- as.vector(temps_cesm_ssp126[[1]][[4]])
temp_mat <- matrix(temp_vec, nrow = 182 * 258, ncol = 261)
lonlat <- as.matrix(expand.grid(temps_cesm_ssp126[[1]][[1]], temps_cesm_ssp126[[1]][[2]]))
temp_df <- data.frame(cbind(lonlat, temp_mat))

# subtract 


# save
saveRDS(yfs_complete_egg, file = here('data', 'yfs_egg.rds'))
saveRDS(yfs_complete_larvae, file = here('data', 'yfs_larvae.rds'))


### Alaska Plaice ----
# change to lowercase
akp_egg_formatted <- format_data(akp_egg_raw)
akp_larvae_formatted <- format_data(akp_larvae_raw)

# Check attributes of egg and larval data
data_check(akp_egg_formatted, akp_larvae_formatted)

table(akp_egg_formatted$haul_performance)
table(akp_larvae_formatted$haul_performance)

# Clean up data to remove unnecessary hauls
akp_egg_clean <- clean_data(akp_egg_formatted)
akp_larvae_clean <- clean_data(akp_larvae_formatted)

table(akp_egg_clean$primary_net)
table(akp_larvae_clean$primary_net)

# Trim egg and larval data
# Year: 1988 forward, 1997 and 2011 for sporadic sampling
# Month: 4 - 6 for eggs, 5 - 7 for larvae
# Latitude: all
akp_egg_trim <- trim_data(akp_egg_clean)
akp_egg_trim <- filter(akp_egg_trim,
                       month > 3, month < 7,
                       year != 1997, year != 2011)
akp_larvae_trim <- trim_data(akp_larvae_clean)
akp_larvae_trim <- filter(akp_larvae_trim,
                          month > 4, month < 8,
                          year != 1997, year != 2011)

# Inspect new data
data_check(akp_egg_trim, akp_larvae_trim)

# Select data for constrained analyses, including stations that are <30km away from the closest positive catch
akp_subset_egg <- final_data(akp_egg_trim)
akp_subset_larvae <- final_data(akp_larvae_trim)
table(akp_subset_egg$year)
table(akp_subset_larvae$year)

# Split egg and larvae data into separate datasets by date
akp_subset_egg1 <- filter(akp_subset_egg, year < 1990)
akp_complete_larvae1 <- filter(akp_subset_larvae, year < 1990)
akp_subset_egg2 <- filter(akp_subset_egg, year > 1989 & year < 1995)
akp_complete_larvae2 <- filter(akp_subset_larvae, year > 1989 & year < 1995)
akp_subset_egg3 <- filter(akp_subset_egg, year > 1994 & year < 2000)
akp_complete_larvae3 <- filter(akp_subset_larvae, year > 1994 & year < 2000)
akp_subset_egg4 <- filter(akp_subset_egg, year > 1999 & year < 2005)
akp_complete_larvae4 <- filter(akp_subset_larvae, year > 1999 & year < 2005)
akp_subset_egg5 <- filter(akp_subset_egg, year > 2004 & year < 2010)
akp_complete_larvae5 <- filter(akp_subset_larvae, year > 2004 & year < 2010)
akp_subset_egg6 <- filter(akp_subset_egg, year > 2009 & year < 2015)
akp_complete_larvae6 <- filter(akp_subset_larvae, year > 2009 & year < 2015)
akp_subset_egg7 <- filter(akp_subset_egg, year > 2014 & year < 2020)
akp_complete_larvae7 <- filter(akp_subset_larvae, year > 2014 & year < 2020)

# Add Bering10K model temperatures and salinities
akp_complete_egg1 <- varid_match(akp_subset_egg1, temp_output1, salt_output1)
akp_complete_egg2 <- varid_match(akp_subset_egg2, temp_output2, salt_output2)
akp_complete_egg3 <- varid_match(akp_subset_egg3, temp_output3, salt_output3)
akp_complete_egg4 <- varid_match(akp_subset_egg4, temp_output4, salt_output4)
akp_complete_egg5 <- varid_match(akp_subset_egg5, temp_output5, salt_output5)
akp_complete_egg6 <- varid_match(akp_subset_egg6, temp_output6, salt_output6)
akp_complete_egg7 <- varid_match(akp_subset_egg7, temp_output7, salt_output7)

akp_complete_larvae1 <- varid_match(akp_complete_larvae1, temp_output1, salt_output1)
akp_complete_larvae2 <- varid_match(akp_complete_larvae2, temp_output2, salt_output2)
akp_complete_larvae3 <- varid_match(akp_complete_larvae3, temp_output3, salt_output3)
akp_complete_larvae4 <- varid_match(akp_complete_larvae4, temp_output4, salt_output4)
akp_complete_larvae5 <- varid_match(akp_complete_larvae5, temp_output5, salt_output5)
akp_complete_larvae6 <- varid_match(akp_complete_larvae6, temp_output6, salt_output6)
akp_complete_larvae7 <- varid_match(akp_complete_larvae7, temp_output7, salt_output7)

# combine to create whole datasets for each
akp_complete_egg <- rbind(akp_complete_egg1, akp_complete_egg2, akp_complete_egg3,
                          akp_complete_egg4, akp_complete_egg5, akp_complete_egg6,
                          akp_complete_egg7)
akp_complete_larvae <- rbind(akp_complete_larvae1, akp_complete_larvae2, akp_complete_larvae3,
                             akp_complete_larvae4, akp_complete_larvae5, akp_complete_larvae6,
                             akp_complete_larvae7)

# add counts
akp_complete_egg$count <- round(akp_complete_egg$larvalcatchper1000m3 * 
                                  akp_complete_egg$volume_filtered / 1000, 0)
akp_complete_larvae$count <- round(akp_complete_larvae$larvalcatchper1000m3 * 
                                     akp_complete_larvae$volume_filtered / 1000, 0)

saveRDS(akp_complete_egg, file = here('data', 'akp_egg.rds'))
saveRDS(akp_complete_larvae, file = here('data', 'akp_larvae.rds'))


### Flathead Sole ----
# change to lowercase
fhs_egg_formatted <- format_data(fhs_egg_raw)
fhs_larvae_formatted <- format_data(fhs_larvae_raw)

# Check attributes of egg and larval data
data_check(fhs_egg_formatted, fhs_larvae_formatted)

table(fhs_egg_formatted$haul_performance)
table(fhs_larvae_formatted$haul_performance)

# Clean up data to remove unnecessary hauls
fhs_egg_clean <- clean_data(fhs_egg_formatted)
fhs_larvae_clean <- clean_data(fhs_larvae_formatted)

table(fhs_egg_clean$primary_net)
table(fhs_larvae_clean$primary_net)

# Trim egg and larval data
# Year: 1988 forward, remove 1997 and 2011 for sporadic samping
# Month: 5 - 9 for eggs, 4 - 7 for larvae
# Latitude: all
fhs_egg_trim <- trim_data(fhs_egg_clean)
fhs_egg_trim <- filter(fhs_egg_trim,
                       month > 4, month < 10,
                       year != 1997, year != 2011)
fhs_larvae_trim <- trim_data(fhs_larvae_clean)
fhs_larvae_trim <- filter(fhs_larvae_trim,
                          month > 4, month < 8,
                          year != 1997, year != 2011)

# Inspect new data
data_check(fhs_egg_trim, fhs_larvae_trim)

# Select data for constrained analyses, including stations that are <30km away from the closest positive catch
fhs_subset_egg <- final_data(fhs_egg_trim)
fhs_subset_larvae <- final_data(fhs_larvae_trim)
table(fhs_subset_egg$year)
table(fhs_subset_larvae$year)

# Split egg and larvae data into separate datasets by date
fhs_subset_egg1 <- filter(fhs_subset_egg, year < 1990)
fhs_complete_larvae1 <- filter(fhs_subset_larvae, year < 1990)
fhs_subset_egg2 <- filter(fhs_subset_egg, year > 1989 & year < 1995)
fhs_complete_larvae2 <- filter(fhs_subset_larvae, year > 1989 & year < 1995)
fhs_subset_egg3 <- filter(fhs_subset_egg, year > 1994 & year < 2000)
fhs_complete_larvae3 <- filter(fhs_subset_larvae, year > 1994 & year < 2000)
fhs_subset_egg4 <- filter(fhs_subset_egg, year > 1999 & year < 2005)
fhs_complete_larvae4 <- filter(fhs_subset_larvae, year > 1999 & year < 2005)
fhs_subset_egg5 <- filter(fhs_subset_egg, year > 2004 & year < 2010)
fhs_complete_larvae5 <- filter(fhs_subset_larvae, year > 2004 & year < 2010)
fhs_subset_egg6 <- filter(fhs_subset_egg, year > 2009 & year < 2015)
fhs_complete_larvae6 <- filter(fhs_subset_larvae, year > 2009 & year < 2015)
fhs_subset_egg7 <- filter(fhs_subset_egg, year > 2014 & year < 2020)
fhs_complete_larvae7 <- filter(fhs_subset_larvae, year > 2014 & year < 2020)

# Add Bering10K model temperatures and salinities
fhs_complete_egg1 <- varid_match(fhs_subset_egg1, temp_output1, salt_output1)
fhs_complete_egg2 <- varid_match(fhs_subset_egg2, temp_output2, salt_output2)
fhs_complete_egg3 <- varid_match(fhs_subset_egg3, temp_output3, salt_output3)
fhs_complete_egg4 <- varid_match(fhs_subset_egg4, temp_output4, salt_output4)
fhs_complete_egg5 <- varid_match(fhs_subset_egg5, temp_output5, salt_output5)
fhs_complete_egg6 <- varid_match(fhs_subset_egg6, temp_output6, salt_output6)
fhs_complete_egg7 <- varid_match(fhs_subset_egg7, temp_output7, salt_output7)

fhs_complete_larvae1 <- varid_match(fhs_complete_larvae1, temp_output1, salt_output1)
fhs_complete_larvae2 <- varid_match(fhs_complete_larvae2, temp_output2, salt_output2)
fhs_complete_larvae3 <- varid_match(fhs_complete_larvae3, temp_output3, salt_output3)
fhs_complete_larvae4 <- varid_match(fhs_complete_larvae4, temp_output4, salt_output4)
fhs_complete_larvae5 <- varid_match(fhs_complete_larvae5, temp_output5, salt_output5)
fhs_complete_larvae6 <- varid_match(fhs_complete_larvae6, temp_output6, salt_output6)
fhs_complete_larvae7 <- varid_match(fhs_complete_larvae7, temp_output7, salt_output7)

# combine to create whole datasets for each
fhs_complete_egg <- rbind(fhs_complete_egg1, fhs_complete_egg2, fhs_complete_egg3,
                          fhs_complete_egg4, fhs_complete_egg5, fhs_complete_egg6,
                          fhs_complete_egg7)
fhs_complete_larvae <- rbind(fhs_complete_larvae1, fhs_complete_larvae2, fhs_complete_larvae3,
                             fhs_complete_larvae4, fhs_complete_larvae5, fhs_complete_larvae6,
                             fhs_complete_larvae7)

# add counts
fhs_complete_egg$count <- round(fhs_complete_egg$larvalcatchper1000m3 * 
                                  fhs_complete_egg$volume_filtered / 1000, 0)
fhs_complete_larvae$count <- round(fhs_complete_larvae$larvalcatchper1000m3 * 
                                     fhs_complete_larvae$volume_filtered / 1000, 0)

saveRDS(fhs_complete_egg, file = here('data', 'fhs_egg.rds'))
saveRDS(fhs_complete_larvae, file = here('data', 'fhs_larvae.rds'))


### Walleye Pollock ----
# change to lowercase
pk_egg_formatted <- format_data(pk_egg_raw)
pk_larvae_formatted <- format_data(pk_larvae_raw)

# Check attributes of egg and larval data
data_check(pk_egg_formatted, pk_larvae_formatted)

table(pk_egg_formatted$haul_performance)
table(pk_larvae_formatted$haul_performance)

# Clean up data to remove unnecessary hauls
pk_egg_clean <- clean_data(pk_egg_formatted)
pk_larvae_clean <- clean_data(pk_larvae_formatted)

table(pk_egg_clean$primary_net)
table(pk_larvae_clean$primary_net)

# Trim egg and larval data
# Year: 1988 forward, remove 1997 and 2011 as well
# Month: 3 - 7
# Latitude: all
pk_egg_trim <- trim_data(pk_egg_clean)
pk_egg_trim <- filter(pk_egg_trim,
                      month > 3, month < 7, 
                      year != 1996, 
                      year != 1997, 
                      year != 1998, 
                      year != 2000,
                      year != 2013)
pk_larvae_trim <- trim_data(pk_larvae_clean)
pk_larvae_trim <- filter(pk_larvae_trim,
                         month > 3, month < 7,
                         year != 1996, 
                         year != 1997, 
                         year != 1998, 
                         year != 2000,
                         year != 2013)

# Inspect new data
data_check(pk_egg_trim, pk_larvae_trim)

# Select data for constrained analyses, including stations that are <30km away from the closest positive catch
pk_subset_egg <- final_data(pk_egg_trim)
pk_subset_larvae <- final_data(pk_larvae_trim)
table(pk_subset_egg$year)
table(pk_subset_larvae$year)

# Split egg and larvae data into separate datasets by date
pk_subset_egg1 <- filter(pk_subset_egg, year < 1990)
pk_complete_larvae1 <- filter(pk_subset_larvae, year < 1990)
pk_subset_egg2 <- filter(pk_subset_egg, year > 1989 & year < 1995)
pk_complete_larvae2 <- filter(pk_subset_larvae, year > 1989 & year < 1995)
pk_subset_egg3 <- filter(pk_subset_egg, year > 1994 & year < 2000)
pk_complete_larvae3 <- filter(pk_subset_larvae, year > 1994 & year < 2000)
pk_subset_egg4 <- filter(pk_subset_egg, year > 1999 & year < 2005)
pk_complete_larvae4 <- filter(pk_subset_larvae, year > 1999 & year < 2005)
pk_subset_egg5 <- filter(pk_subset_egg, year > 2004 & year < 2010)
pk_complete_larvae5 <- filter(pk_subset_larvae, year > 2004 & year < 2010)
pk_subset_egg6 <- filter(pk_subset_egg, year > 2009 & year < 2015)
pk_complete_larvae6 <- filter(pk_subset_larvae, year > 2009 & year < 2015)
pk_subset_egg7 <- filter(pk_subset_egg, year > 2014 & year < 2020)
pk_complete_larvae7 <- filter(pk_subset_larvae, year > 2014 & year < 2020)

# Add Bering10K model temperatures and salinities
pk_complete_egg1 <- varid_match(pk_subset_egg1, temp_output1, salt_output1)
pk_complete_egg2 <- varid_match(pk_subset_egg2, temp_output2, salt_output2)
pk_complete_egg3 <- varid_match(pk_subset_egg3, temp_output3, salt_output3)
pk_complete_egg4 <- varid_match(pk_subset_egg4, temp_output4, salt_output4)
pk_complete_egg5 <- varid_match(pk_subset_egg5, temp_output5, salt_output5)
pk_complete_egg6 <- varid_match(pk_subset_egg6, temp_output6, salt_output6)
pk_complete_egg7 <- varid_match(pk_subset_egg7, temp_output7, salt_output7)

pk_complete_larvae1 <- varid_match(pk_complete_larvae1, temp_output1, salt_output1)
pk_complete_larvae2 <- varid_match(pk_complete_larvae2, temp_output2, salt_output2)
pk_complete_larvae3 <- varid_match(pk_complete_larvae3, temp_output3, salt_output3)
pk_complete_larvae4 <- varid_match(pk_complete_larvae4, temp_output4, salt_output4)
pk_complete_larvae5 <- varid_match(pk_complete_larvae5, temp_output5, salt_output5)
pk_complete_larvae6 <- varid_match(pk_complete_larvae6, temp_output6, salt_output6)
pk_complete_larvae7 <- varid_match(pk_complete_larvae7, temp_output7, salt_output7)

# combine to create whole datasets for each
pk_complete_egg <- rbind(pk_complete_egg1, pk_complete_egg2, pk_complete_egg3,
                         pk_complete_egg4, pk_complete_egg5, pk_complete_egg6,
                         pk_complete_egg7)
pk_complete_larvae <- rbind(pk_complete_larvae1, pk_complete_larvae2, pk_complete_larvae3,
                            pk_complete_larvae4, pk_complete_larvae5, pk_complete_larvae6,
                            pk_complete_larvae7)

# add counts
pk_complete_egg$count <- round(pk_complete_egg$larvalcatchper1000m3 * 
                                 pk_complete_egg$volume_filtered / 1000, 0)
pk_complete_larvae$count <- round(pk_complete_larvae$larvalcatchper1000m3 * 
                                    pk_complete_larvae$volume_filtered / 1000, 0)

saveRDS(pk_complete_egg, file = here('data', 'pk_egg.rds'))
saveRDS(pk_complete_larvae, file = here('data', 'pk_larvae.rds'))


### Pacific Cod ----
# change to lowercase
pcod_larvae_formatted <- format_data(pcod_larvae_raw)

# Check attributes of egg and larval data
data_check(pk_egg_formatted, pcod_larvae_formatted)

table(pcod_larvae_formatted$haul_performance)

# Clean up data to remove unnecessary hauls
pcod_larvae_clean <- clean_data(pcod_larvae_formatted)

table(pcod_larvae_clean$primary_net)

# Trim larval data
# Year: 1988 forward, remove 1997 and 2011 as well
# Month: 4 - 6
# Latitude: all
pcod_larvae_trim <- trim_data(pcod_larvae_clean)
pcod_larvae_trim <- filter(pcod_larvae_trim,
                         month > 3, month < 7,
                         year != 1997, year != 2011)

# Inspect new data
data_check(pk_egg_trim, pcod_larvae_trim)

# Select data for constrained analyses, including stations that are <30km away from the closest positive catch
pcod_subset_larvae <- final_data(pcod_larvae_trim)
table(pcod_subset_larvae$year)

# Split egg and larvae data into separate datasets by date
pcod_complete_larvae1 <- filter(pcod_subset_larvae, year < 1990)
pcod_complete_larvae2 <- filter(pcod_subset_larvae, year > 1989 & year < 1995)
pcod_complete_larvae3 <- filter(pcod_subset_larvae, year > 1994 & year < 2000)
pcod_complete_larvae4 <- filter(pcod_subset_larvae, year > 1999 & year < 2005)
pcod_complete_larvae5 <- filter(pcod_subset_larvae, year > 2004 & year < 2010)
pcod_complete_larvae6 <- filter(pcod_subset_larvae, year > 2009 & year < 2015)
pcod_complete_larvae7 <- filter(pcod_subset_larvae, year > 2014 & year < 2020)

# Add Bering10K model temperatures and salinities
pcod_complete_larvae1 <- varid_match(pcod_complete_larvae1, temp_output1, salt_output1)
pcod_complete_larvae2 <- varid_match(pcod_complete_larvae2, temp_output2, salt_output2)
pcod_complete_larvae3 <- varid_match(pcod_complete_larvae3, temp_output3, salt_output3)
pcod_complete_larvae4 <- varid_match(pcod_complete_larvae4, temp_output4, salt_output4)
pcod_complete_larvae5 <- varid_match(pcod_complete_larvae5, temp_output5, salt_output5)
pcod_complete_larvae6 <- varid_match(pcod_complete_larvae6, temp_output6, salt_output6)
pcod_complete_larvae7 <- varid_match(pcod_complete_larvae7, temp_output7, salt_output7)

# combine to create whole datasets for each
pcod_complete_larvae <- rbind(pcod_complete_larvae1, pcod_complete_larvae2, pcod_complete_larvae3,
                            pcod_complete_larvae4, pcod_complete_larvae5, pcod_complete_larvae6,
                            pcod_complete_larvae7)

# add counts
pcod_complete_larvae$count <- round(pcod_complete_larvae$larvalcatchper1000m3 * 
                                    pcod_complete_larvae$volume_filtered / 1000, 0)

saveRDS(pcod_complete_larvae, file = here('data', 'pcod_larvae.rds'))

### Northern Rock Sole ----
# change to lowercase
nrs_larvae_formatted <- format_data(nrs_larvae_raw)

# Check attributes of egg and larval data
data_check(pk_egg_formatted, nrs_larvae_formatted)

table(nrs_larvae_formatted$haul_performance)

# Clean up data to remove unnecessary hauls
nrs_larvae_clean <- clean_data(nrs_larvae_formatted)

table(nrs_larvae_clean$primary_net)

# Trim larval data
# Year: 1988 forward, remove 1997 and 2011 as well
# Month: 4 - 6
# Latitude: all
nrs_larvae_trim <- trim_data(nrs_larvae_clean)
nrs_larvae_trim <- filter(nrs_larvae_trim,
                           month > 3, month < 7,
                           year != 1997, year != 2011)

# Inspect new data
data_check(pk_egg_trim, nrs_larvae_trim)

# Select data for constrained analyses, including stations that are <30km away from the closest positive catch
nrs_subset_larvae <- final_data(nrs_larvae_trim)
table(nrs_subset_larvae$year)

# Split egg and larvae data into separate datasets by date
nrs_complete_larvae1 <- filter(nrs_subset_larvae, year < 1990)
nrs_complete_larvae2 <- filter(nrs_subset_larvae, year > 1989 & year < 1995)
nrs_complete_larvae3 <- filter(nrs_subset_larvae, year > 1994 & year < 2000)
nrs_complete_larvae4 <- filter(nrs_subset_larvae, year > 1999 & year < 2005)
nrs_complete_larvae5 <- filter(nrs_subset_larvae, year > 2004 & year < 2010)
nrs_complete_larvae6 <- filter(nrs_subset_larvae, year > 2009 & year < 2015)
nrs_complete_larvae7 <- filter(nrs_subset_larvae, year > 2014 & year < 2020)

# Add Bering10K model temperatures and salinities
nrs_complete_larvae1 <- varid_match(nrs_complete_larvae1, temp_output1, salt_output1)
nrs_complete_larvae2 <- varid_match(nrs_complete_larvae2, temp_output2, salt_output2)
nrs_complete_larvae3 <- varid_match(nrs_complete_larvae3, temp_output3, salt_output3)
nrs_complete_larvae4 <- varid_match(nrs_complete_larvae4, temp_output4, salt_output4)
nrs_complete_larvae5 <- varid_match(nrs_complete_larvae5, temp_output5, salt_output5)
nrs_complete_larvae6 <- varid_match(nrs_complete_larvae6, temp_output6, salt_output6)
nrs_complete_larvae7 <- varid_match(nrs_complete_larvae7, temp_output7, salt_output7)

# combine to create whole datasets for each
nrs_complete_larvae <- rbind(nrs_complete_larvae1, nrs_complete_larvae2, nrs_complete_larvae3,
                              nrs_complete_larvae4, nrs_complete_larvae5, nrs_complete_larvae6,
                              nrs_complete_larvae7)

# add counts
nrs_complete_larvae$count <- round(nrs_complete_larvae$larvalcatchper1000m3 * 
                                      nrs_complete_larvae$volume_filtered / 1000, 0)

saveRDS(nrs_complete_larvae, file = here('data', 'nrs_larvae.rds'))
