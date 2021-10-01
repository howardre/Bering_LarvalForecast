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

bering_model <- tidync(here('data/temperature_netcdf', 'B10K-K20_CORECFS_2015-2019_average_temp_surface5m.nc'))
time_ex <- bering_model %>% 
  activate("D2") %>% 
  hyper_array()

location_ex <- bering_model %>%
  activate("D1,D0") %>%
  hyper_array()

time1 <- as.Date(time_ex$ocean_time / (60 * 60 * 24), origin = "1900-01-01 00:00:00")
bering_data <- bering_model %>% activate("D1,D0") %>% 
  hyper_array()
transform <- attr(bering_data, "transforms")

latrange <- c(55, 59)

bering_slice <- bering_model %>% 
  activate("D1,D0") %>%
  hyper_filter(eta_rho = eta_rho > -170,
               xi_rho = xi_rho > latrange[1] & xi_rho <= latrange[2])
bering_slice_data <- bering_slice %>% hyper_array()
slice_transfrom <- attr(bering_slice_data, "transforms")
lon <- slice_transfrom$lon_rho %>% filter(selected)
lat <- slice_transfrom$xi_rho %>% filter(selected)

bering_model %>% hyper_filter()
bering_model %>% hyper_filter(xi_rho = xi_rho < 59 & xi_rho > 55)
bering_filtered <- bering_model %>% hyper_filter(xi_rho = dplyr::between(index, 55, 59),
                              eta_rho = index > -170)
bering_array <- bering_filtered %>% hyper_array()
str(bering_array)

tran <- attr(bering_array, "transforms")
lat <- tran$xi_rho %>% filter(selected)
lon <- tran$eta_rho %>% filter(selected)

image(lon$xi_rho)

hyper_vars(bering_filtered)




location <- cbind(lon = as.vector(lon1), lat = as.vector(lat), temp = as.vector(temp_plot))
loc_filter <- location %>% filter(lon > -170)



lon <- ncvar_get(bering_model, varid = 'lon_rho')
lon1 <- ifelse(lon >= 180, lon -360, lon)
lat <- ncvar_get(bering_model, varid = 'lat_rho')
time <- ncvar_get(bering_model, varid = 'ocean_time')
time1 <- as.Date(time / (60 * 60 * 24), origin = "1900-01-01 00:00:00")
fillvalue_t <- ncatt_get(bering_model, 'temp', "_FillValue")
temp <- ncvar_get(bering_model, varid = 'temp')
temp[temp == fillvalue_t$value] <- NA
dimnames(temp) <- list(lon = lon1, lat = lat, time = time1)

dim(temp)
range(temp, na.rm = T)
dim(lat)
range(lat, na.rm = T)
dim(lon)
range(lon, na.rm = T)
dim(time1)
range(time1, na.rm = T)
nc_close(bering_model)

### Generate monthly values by year from surface temps -----
roms_temps <- data.frame(year = as.numeric(c('1988', '1991', '1992', '1993',
                                             '1994', '1995', '1996', '1998',
                                             '1999', '2000', '2001', '2002',
                                             '2003', '2004', '2005', '2006',
                                             '2007', '2008', '2009', '2010',
                                             '2012', '2013', '2014', '2015', 
                                             '2016')),
                         '01' = as.numeric(NA),
                         '02' = as.numeric(NA),
                         '03' = as.numeric(NA),
                         '04' = as.numeric(NA),
                         '05' = as.numeric(NA),
                         '06' = as.numeric(NA),
                         '07' = as.numeric(NA),
                         '08' = as.numeric(NA),
                         '09' = as.numeric(NA),
                         '10' = as.numeric(NA),
                         '11' = as.numeric(NA),
                         '12' = as.numeric(NA))

year <- unique(substr(temp_output1[[3]], 1, 4)) # contains unique years in file
month_year <- substr(temp_output1[[3]], 1, 7) # contains month year combinations
month <- "03"

time_index <- time1 >= "2015-02-01" & time1 <= "2015-04-30"
lon_index <- lon1 <= -150
lat_index <- lat >= 55 & lat <= 59
temp_plot <- temp[, , time_index]

location <- cbind(lon = as.vector(lon1), lat = as.vector(lat), temp = as.vector(temp_plot))
loc_filter <- location %>% select(lon > -170)

dim(temp)
dim(lat_index)
dim(temp_plot)

lon_plot <- c(lon1[lon_index])
lat_plot <- c(lat[lat_index])

image.plot(lon_plot,
           sort(lat_plot),
           temp_plot[, order(lat_plot)])

# attempt to extract for month
date <- "1988-05-01"
sst <- temp_output1[[4]][, , which(temp_output1[[3]] == date)] %>%
  melt(varnames = c("lon", "lat")) %>%
  subset(!is.na(value))

data_frame <- data.frame(date = temp_output1[[3]],
                         dp = temp_output1[[4]][as.character()])


windows(width = 9, height = 10)
par(mfrow = c(3, 2))
for (i in 1:length(year)) {
  starttime <- (1:length(month_year))[month_year == paste(year[i], month, sep = "-")][1]
  
  # Get temp
  tempplot <- temp_output1[[4]][, , starttime]
  
  # Get date of image plot
  plotday1 <- as.character(temp_output1[[3]][starttime])
  
  
  
  # show Temp
  image.plot(lon,
             lat,
             tempplot,
             ylab = "Latitude (north)",
             col = viridis(100),
             xlab = "Longitude (east)",
             main = plotday1,
             xlim = c(-180,-155) + 360,
             ylim = c(54, 64),
             zlim = c(-2, 17.5))
  maps::map("world2",
            fill = T,
            col = "grey",
            add = T)
}


roms_year <- unique(roms_temps$year)

for (i in 1:length(roms_year)) {
  idx_time <- (1:length(month_year))[month_year == paste(year[i], month, sep = "-")][1]
  # Get temp
  temp_plot <- temp[, , idx_time]
  # Get start and end date of data fields
  plot_day1 <- as.character(time1[idx_time])
  # Make data frame for loess
  temp_data <- na.exclude(data.frame(z = c(temp_plot),
                                      x = c(lon1),
                                      y = c(lat)))
  
  roms_temps$loess_temperature[roms_temps$year == year[i]] <- predict(temp_loess,
                                                              newdata = pred_loess)
  roms_temps$loess_date[roms_temps$year == year[i]] <- plot_day1

}

# Make box from -170 to -165 lon and 59 to 55 lat
polygon <- c(-170, -165, 59, 55)
temp1_data <- temp_output1[sapply(temp_output1, function(x) x[[1]] >= -170 | x[[1]] <= -150 )]