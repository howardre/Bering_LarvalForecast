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
library(lattice)
library(tidync)
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

# Bering 10K model output: avg bottom temp 2015-2019
bering_model <- nc_open(here('data/temperature_netcdf', 'B10K-K20_CORECFS_2015-2019_average_temp_surface5m.nc'))

lon <- ncvar_get(bering_model, varid = 'lon_rho')
lon1 <- ifelse(lon >= 180, lon -360, lon)
lat <- ncvar_get(bering_model, varid = 'lat_rho')
time <- ncvar_get(bering_model, varid = 'ocean_time')
time1 <- as.Date(time / (60 * 60 * 24), origin = "1900-01-01 00:00:00")
fillvalue_t <- ncatt_get(bering_model, 'temp', "_FillValue")
temp <- ncvar_get(bering_model, varid = 'temp')
temp[temp == fillvalue_t$value] <- NA

nc_close(bering_model)

### Generate mean values by year from surface temps -----
# Select out time period
time_index <- time1 >= "2016-02-01" & time1 <= "2016-04-30"
temp_array <- temp[, , time_index]

# Select out box on the shelf
temp_data <- as.data.frame(cbind(lon = as.vector(lon1), lat = as.vector(lat), temp = as.vector(temp_array)))
temp_filtered <- location %>% filter(lon >= -170 & lon <= -165, lat >= 55 & lat <= 59)

# Create empty data frame
roms_temps <- data.frame(year = as.numeric(c('1988', '1991', '1992', '1993',
                                             '1994', '1995', '1996', '1998',
                                             '1999', '2000', '2001', '2002',
                                             '2003', '2004', '2005', '2006',
                                             '2007', '2008', '2009', '2010',
                                             '2012', '2013', '2014', '2015', 
                                             '2016')),
                         'mean' = as.numeric(NA))


# Check to make sure logical
temp_filtered %>%
  slice(1:10000) %>%
ggplot(aes(x = lon, y = lat, color = temp)) +
  geom_point(size = 5, alpha = 0.5) +
  scale_color_viridis_b()

