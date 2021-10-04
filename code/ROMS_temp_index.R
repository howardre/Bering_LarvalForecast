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
library(ggplot2)

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

### Generate mean values by year from surface temps -----
mean_extract <- function(data, start_date, end_date){
  # Select out time period
  time_index <- data[[3]] >= start_date & data[[3]] <= end_date
  temp_array <- data[[4]][, , time_index]
  
  # Select out box on the shelf
  temp_data <- as.data.frame(cbind(lon = as.vector(data[[1]]), 
                                   lat = as.vector(data[[2]]), 
                                   temp = as.vector(temp_array)))
  temp_filtered <- temp_data %>% filter(lon >= -170 & lon <= -165, lat >= 56 & lat <= 58)
  mean <- mean(temp_filtered$temp, na.rm = T)
  return(mean)
  }

# Create empty data frame
roms_temps <- data.frame(year = as.numeric(c('1988', '1991', '1992', '1993',
                                             '1994', '1995', '1996', '1998',
                                             '1999', '2000', '2001', '2002',
                                             '2003', '2004', '2005', '2006',
                                             '2007', '2008', '2009', '2010',
                                             '2012', '2013', '2014', '2015', 
                                             '2016')),
                         'mean' = as.numeric(NA))

roms_temps[roms_temps$year == 1988, 'mean'] <- mean_extract(temp_output1, "1988-02-01", "1988-04-30") 
roms_temps[roms_temps$year == 1991, 'mean'] <- mean_extract(temp_output2, "1991-02-01", "1991-04-30") 
roms_temps[roms_temps$year == 1992, 'mean'] <- mean_extract(temp_output2, "1992-02-01", "1992-04-30") 
roms_temps[roms_temps$year == 1993, 'mean'] <- mean_extract(temp_output2, "1993-02-01", "1993-04-30") 
roms_temps[roms_temps$year == 1994, 'mean'] <- mean_extract(temp_output2, "1994-02-01", "1994-04-30") 
roms_temps[roms_temps$year == 1995, 'mean'] <- mean_extract(temp_output3, "1995-02-01", "1995-04-30") 
roms_temps[roms_temps$year == 1996, 'mean'] <- mean_extract(temp_output3, "1996-02-01", "1996-04-30") 
roms_temps[roms_temps$year == 1998, 'mean'] <- mean_extract(temp_output3, "1998-02-01", "1998-04-30") 
roms_temps[roms_temps$year == 1999, 'mean'] <- mean_extract(temp_output3, "1999-02-01", "1999-04-30")
roms_temps[roms_temps$year == 2000, 'mean'] <- mean_extract(temp_output4, "2000-02-01", "2000-04-30")
roms_temps[roms_temps$year == 2001, 'mean'] <- mean_extract(temp_output4, "2001-02-01", "2001-04-30") 
roms_temps[roms_temps$year == 2002, 'mean'] <- mean_extract(temp_output4, "2002-02-01", "2002-04-30") 
roms_temps[roms_temps$year == 2003, 'mean'] <- mean_extract(temp_output4, "2003-02-01", "2003-04-30") 
roms_temps[roms_temps$year == 2004, 'mean'] <- mean_extract(temp_output4, "2004-02-01", "2004-04-30") 
roms_temps[roms_temps$year == 2005, 'mean'] <- mean_extract(temp_output5, "2005-02-01", "2005-04-30") 
roms_temps[roms_temps$year == 2006, 'mean'] <- mean_extract(temp_output5, "2006-02-01", "2006-04-30") 
roms_temps[roms_temps$year == 2007, 'mean'] <- mean_extract(temp_output5, "2007-02-01", "2007-04-30")
roms_temps[roms_temps$year == 2008, 'mean'] <- mean_extract(temp_output5, "2008-02-01", "2008-04-30") 
roms_temps[roms_temps$year == 2009, 'mean'] <- mean_extract(temp_output5, "2009-02-01", "2009-04-30")
roms_temps[roms_temps$year == 2010, 'mean'] <- mean_extract(temp_output6, "2010-02-01", "2010-04-30")
roms_temps[roms_temps$year == 2011, 'mean'] <- mean_extract(temp_output6, "2011-02-01", "2011-04-30") 
roms_temps[roms_temps$year == 2012, 'mean'] <- mean_extract(temp_output6, "2012-02-01", "2012-04-30") 
roms_temps[roms_temps$year == 2013, 'mean'] <- mean_extract(temp_output6, "2013-02-01", "2013-04-30") 
roms_temps[roms_temps$year == 2014, 'mean'] <- mean_extract(temp_output6, "2014-02-01", "2014-04-30") 
roms_temps[roms_temps$year == 2015, 'mean'] <- mean_extract(temp_output7, "2015-02-01", "2015-04-30") 
roms_temps[roms_temps$year == 2016, 'mean'] <- mean_extract(temp_output7, "2016-02-01", "2016-04-30") 

saveRDS(roms_temps, file = here('data', 'roms_temps.rds'))

# Check to make sure logical - use code in function to get the temp_filtered df
# temp_filtered %>%
#   slice(1:10000) %>%
# ggplot(aes(x = lon, y = lat, color = temp)) +
#   geom_point(size = 5, alpha = 0.5) +
#   scale_color_viridis_b()

# Plot the mean temps over year
ggplot(data = roms_temps) +
  geom_point(aes(x = year, y = mean), color = "aquamarine4", size = 3) +
  labs(title = "Mean ROMS Surface Temperature",
       y = "Temperature (C)",
       x = "Year") +
  theme_bw() +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15),
        plot.title = element_text(size = 24),
        text = element_text(family = "serif"))

dev.copy(jpeg, here('results', 'mean_roms_temps.jpg'), 
         height = 10, width = 10, units = 'in', res = 200)
dev.off()
