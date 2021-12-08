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

get_temp_filepath <- function(years){
  filepath = paste('D:/B10K-K20P19_CMIP6_gfdl_ssp126/Level1/B10K-K20P19_CMIP6_gfdl_ssp126_', years, '_average_temp.nc', sep = '')
  return(filepath)
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

years_list <- list('2015-2019', '2020-2024', '2025-2029', '2030-2034',
                   '2035-2039', '2040-2044', '2045-2049', '2050-2054',
                   '2055-2059', '2060-2064', '2065-2069', '2070-2074',
                   '2075-2079', '2080-2084', '2085-2089', '2090-2094',
                   '2095-2099')

# Function to get temperature datasets
get_netcdf <- function(years_list, y){
  bering_model_temp = nc_open(get_temp_filepath(years_list[y]))
  temp_output <- nc_extract(bering_model_temp, temp, 'temp')
  nc_close(bering_model_temp)
  return(temp_output)
}

# Get mean temperatures for each year
mean_extract <- function(data, start_date, end_date){
  # Select out time period
  time_index <- data[[3]] >= start_date & data[[3]] <= end_date
  temp_array <- data[[4]][, , , time_index]
  
  # Select out box on the shelf
  temp_data <- as.data.frame(cbind(lon = as.vector(data[[1]]), 
                                   lat = as.vector(data[[2]]), 
                                   temp = as.vector(temp_array)))
  temp_filtered <- temp_data %>% filter(lon >= -170 & lon <= -165, lat >= 56 & lat <= 58)
  mean <- mean(temp_filtered$temp, na.rm = T)
  return(mean)
}

### Bering 10K model output ----
# Using avg surface temperatures
temp_output1 <- get_netcdf(years_list, 1)
temp_output2 <- get_netcdf(years_list, 2)
temp_output3 <- get_netcdf(years_list, 3)
temp_output4 <- get_netcdf(years_list, 4)
temp_output5 <- get_netcdf(years_list, 5)
temp_output6 <- get_netcdf(years_list, 6)
temp_output7 <- get_netcdf(years_list, 7)
temp_output8 <- get_netcdf(years_list, 8)
temp_output9 <- get_netcdf(years_list, 9)
temp_output10 <- get_netcdf(years_list, 10)
temp_output11 <- get_netcdf(years_list, 11)
temp_output12 <- get_netcdf(years_list, 12)
temp_output13 <- get_netcdf(years_list, 13)
temp_output14 <- get_netcdf(years_list, 14)
temp_output15 <- get_netcdf(years_list, 15)
temp_output16 <- get_netcdf(years_list, 16)
temp_output17 <- get_netcdf(years_list, 17)

### Generate mean values by year from surface temps -----
# Create empty data frame
forecast_temps <- data.frame(year = as.numeric(c('2015', '2016', '2017', '2018', '2019', '2020',
                                             '2021', '2022', '2023', '2024', '2025', '2026',
                                             '2027', '2028', '2029', '2030', '2031', '2032',
                                             '2033', '2034', '2035')),
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
