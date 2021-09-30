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

# close netcdf's
nc_close(bering_model_temp1)
nc_close(bering_model_temp2)
nc_close(bering_model_temp3)
nc_close(bering_model_temp4)
nc_close(bering_model_temp5)
nc_close(bering_model_temp6)
nc_close(bering_model_temp7)

### Generate monthly values by year from surface temps -----
roms_temps <- data.frame(year = as.numeric(c('1988', '1991', '1992', '1993',
                                             '1994', '1995', '1996', '1998',
                                             '1999', '2000', '2001', '2002',
                                             '2003', '2004', '2005', '2006',
                                             '2007', '2008', '2009', '2010',
                                             '2012', '2013', '2014', '2015', 
                                             '2016')),
                         jan = as.numeric(NA),
                         feb = as.numeric(NA),
                         mar = as.numeric(NA),
                         apr = as.numeric(NA),
                         may = as.numeric(NA),
                         jun = as.numeric(NA),
                         jul = as.numeric(NA),
                         aug = as.numeric(NA),
                         sep = as.numeric(NA),
                         oct = as.numeric(NA),
                         nov = as.numeric(NA),
                         dec = as.numeric(NA))

year <- unique(substr(temp_output1[[3]], 1, 4))
month_year <- substr(temp_output1[[3]], 1, 7)
month <- "02"

roms_year <- unique(roms_temps$year)
i <- '1991'

for (i in 1:length(roms_year)) {
  idx_time <- (1:length(month_year))[month_year == paste(roms_year[i], month, sep = "-")][1]
  # Get temp
  temp_plot <- temp[, , idx_time]
  # Get start and end date of data fields
  plot_day1 <- as.character(time1[idx_time])
  # Make data frame for loess
  data_loess <- na.exclude(data.frame(z = c(temp_plot),
                                      x = c(lon1),
                                      y = c(lat)))

}

# Make box from -170 to -165 lon and 59 to 55 lat
polygon <- c(-170, -165, 59, 55)
temp1_data <- temp_output1[sapply(temp_output1, function(x) x[[1]] >= -170 | x[[1]] <= -150 )]