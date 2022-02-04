# Title: Forecasting
# Date: complete 7/07/2021

### Libraries, functions, and data ----
library(raster)
library(ncdf4)
library(spacetime)
library(fields)
library(here)
library(tidyverse)
library(lubridate)
library(date)


# Function to extract data from .nc files
nc_extract <- function(file, variable, varid_name){
  lon <- ncvar_get(file, varid = 'lon_rho')
  lon1 <- ifelse(lon >= 180, lon -360, lon)
  lat <- ncvar_get(file, varid = 'lat_rho')
  time <- ncvar_get(file, varid = 'ocean_time')
  time1 <- as.Date(time / (60 * 60 * 24), origin = "1900-01-01 00:00:00")
  density <- ncvar_get(file, varid = 's_rho')
  
  month <- as.numeric(substr(time1, 6, 7))

  fillvalue_t <- ncatt_get(file, varid_name, "_FillValue")
  
  variable <- ncvar_get(file, varid = varid_name)
  variable[variable == fillvalue_t$value] <- NA
  
  variable <- variable[, , 30, ]
  
  return_list <- list("lon1" = lon1, "lat" = lat, "time1" = time1, "variable" = variable)
} 

# Function to extract netcdf
years_list <- list('2010-2014', '2015-2019', '2020-2024', '2025-2029', 
                   '2030-2034', '2035-2039', '2040-2044', '2045-2049', 
                   '2050-2054', '2055-2059', '2060-2064', '2065-2069', 
                   '2070-2074', '2075-2079', '2080-2084', '2085-2089', 
                   '2090-2094', '2095-2099')

years_list_historical <- list('1980-1984', '1985-1989', '1990-1994', '1995-1999',
                              '2000-2004', '2005-2009', '2010-2014')

get_temp_filepath <- function(model, years){
  filepath = paste('D:/B10K-K20P19_CMIP6_', model, '/Level1/B10K-K20P19_CMIP6_', model, '_', years, '_average_temp.nc', sep = '')
  return(filepath)
}

get_salt_filepath <- function(model, years){
  filepath = paste('D:/B10K-K20P19_CMIP6_', model, '/Level1/B10K-K20P19_CMIP6_', model, '_', years, '_average_salt.nc', sep = '')
  return(filepath)
}

### Create lists -----
# May need to remove from environment in order to get through them all
# Depends on amount of memory available
#### CESM SSP126 ----
temps_cesm_ssp126<- list()
for(j in 2:18){
  bering_model_temp = nc_open(get_temp_filepath('cesm_ssp126', years_list[j]))
  temp_output <- nc_extract(bering_model_temp, temp, 'temp')
  nc_close(bering_model_temp)
  gc()
  temps_cesm_ssp126[[paste("temp_output", j, sep = "")]] <- temp_output
}
saveRDS(temps_cesm_ssp126, file = here('data', 'temps_cesm_ssp126.rds'))
rm(temps_cesm_ssp126)

salts_cesm_ssp126<- list()
for(j in 2:18){
  bering_model_salt = nc_open(get_salt_filepath('cesm_ssp126', years_list[j]))
  salt_output <- nc_extract(bering_model_salt, salt, 'salt')
  nc_close(bering_model_salt)
  gc()
  salts_cesm_ssp126[[paste("salt_output", j, sep = "")]] <- salt_output
}
saveRDS(salts_cesm_ssp126, file = here('data', 'salts_cesm_ssp126.rds'))
rm(salts_cesm_ssp126)

#### CESM SSP585 ----
temps_cesm_ssp585<- list()
for(j in 2:18){
  bering_model_temp = nc_open(get_temp_filepath('cesm_ssp585', years_list[j]))
  temp_output <- nc_extract(bering_model_temp, temp, 'temp')
  nc_close(bering_model_temp)
  gc()
  temps_cesm_ssp585[[paste("temp_output", j, sep = "")]] <- temp_output
}
saveRDS(temps_cesm_ssp585, file = here('data', 'temps_cesm_ssp585.rds'))
rm(temps_cesm_ssp585)

salts_cesm_ssp585<- list()
for(j in 2:18){
  bering_model_salt = nc_open(get_salt_filepath('cesm_ssp585', years_list[j]))
  salt_output <- nc_extract(bering_model_salt, salt, 'salt')
  nc_close(bering_model_salt)
  gc()
  salts_cesm_ssp585[[paste("salt_output", j, sep = "")]] <- salt_output
}
saveRDS(salts_cesm_ssp585, file = here('data', 'salts_cesm_ssp585.rds'))
rm(salts_cesm_ssp585)

#### CESM Historical ----
temps_cesm_historical<- list()
for(j in 1:7){
  bering_model_temp = nc_open(get_temp_filepath('cesm_historical', years_list_historical[j]))
  temp_output <- nc_extract(bering_model_temp, temp, 'temp')
  nc_close(bering_model_temp)
  gc()
  temps_cesm_historical[[paste("temp_output", j, sep = "")]] <- temp_output
}
saveRDS(temps_cesm_historical, file = here('data', 'temps_cesm_historical.rds'))
rm(temps_cesm_historical)

salts_cesm_historical<- list()
for(j in 1:7){
  bering_model_salt = nc_open(get_salt_filepath('cesm_historical', years_list_historical[j]))
  salt_output <- nc_extract(bering_model_salt, salt, 'salt')
  nc_close(bering_model_salt)
  gc()
  salts_cesm_historical[[paste("salt_output", j, sep = "")]] <- salt_output
}
saveRDS(salts_cesm_historical, file = here('data', 'salts_cesm_historical.rds'))
rm(salts_cesm_historical)

#### GFDL SSP126 ----
temps_gfdl_ssp126<- list()
for(j in 2:18){
  bering_model_temp = nc_open(get_temp_filepath('gfdl_ssp126', years_list[j]))
  temp_output <- nc_extract(bering_model_temp, temp, 'temp')
  nc_close(bering_model_temp)
  gc()
  temps_gfdl_ssp126[[paste("temp_output", j, sep = "")]] <- temp_output
}
saveRDS(temps_gfdl_ssp126, file = here('data', 'temps_gfdl_ssp126.rds'))
rm(temps_gfdl_ssp126)

salts_gfdl_ssp126<- list()
for(j in 2:18){
  bering_model_salt = nc_open(get_salt_filepath('gfdl_ssp126', years_list[j]))
  salt_output <- nc_extract(bering_model_salt, salt, 'salt')
  nc_close(bering_model_salt)
  gc()
  salts_gfdl_ssp126[[paste("salt_output", j, sep = "")]] <- salt_output
}
saveRDS(salts_gfdl_ssp126, file = here('data', 'salts_gfdl_ssp126.rds'))
rm(salts_gfdl_ssp126)

#### GFDL SSP585 ----
temps_gfdl_ssp585<- list()
for(j in 2:18){
  bering_model_temp = nc_open(get_temp_filepath('gfdl_ssp585', years_list[j]))
  temp_output <- nc_extract(bering_model_temp, temp, 'temp')
  nc_close(bering_model_temp)
  gc()
  temps_gfdl_ssp585[[paste("temp_output", j, sep = "")]] <- temp_output
}
saveRDS(temps_gfdl_ssp585, file = here('data', 'temps_gfdl_ssp585.rds'))
rm(temps_gfdl_ssp585)

salts_gfdl_ssp585<- list()
for(j in 2:18){
  bering_model_salt = nc_open(get_salt_filepath('gfdl_ssp585', years_list[j]))
  salt_output <- nc_extract(bering_model_salt, salt, 'salt')
  nc_close(bering_model_salt)
  gc()
  salts_gfdl_ssp585[[paste("salt_output", j, sep = "")]] <- salt_output
}
saveRDS(salts_gfdl_ssp585, file = here('data', 'salts_gfdl_ssp585.rds'))
rm(salts_gfdl_ssp585)

#### GFDL Historical ----
temps_gfdl_historical<- list()
for(j in 1:7){
  bering_model_temp = nc_open(get_temp_filepath('gfdl_historical', years_list_historical[j]))
  temp_output <- nc_extract(bering_model_temp, temp, 'temp')
  nc_close(bering_model_temp)
  gc()
  temps_gfdl_historical[[paste("temp_output", j, sep = "")]] <- temp_output
}
saveRDS(temps_gfdl_historical, file = here('data', 'temps_gfdl_historical.rds'))
rm(temps_gfdl_historical)

salts_gfdl_historical<- list()
for(j in 1:7){
  bering_model_salt = nc_open(get_salt_filepath('gfdl_historical', years_list_historical[j]))
  salt_output <- nc_extract(bering_model_salt, salt, 'salt')
  nc_close(bering_model_salt)
  gc()
  salts_gfdl_historical[[paste("salt_output", j, sep = "")]] <- salt_output
}
saveRDS(salts_gfdl_historical, file = here('data', 'salts_gfdl_historical.rds'))
rm(salts_gfdl_historical)

#### MIROC SSP126 ----
temps_miroc_ssp126<- list()
for(j in 2:18){
  bering_model_temp = nc_open(get_temp_filepath('miroc_ssp126', years_list[j]))
  temp_output <- nc_extract(bering_model_temp, temp, 'temp')
  nc_close(bering_model_temp)
  gc()
  temps_miroc_ssp126[[paste("temp_output", j, sep = "")]] <- temp_output
}
saveRDS(temps_miroc_ssp126, file = here('data', 'temps_miroc_ssp126.rds'))
rm(temps_miroc_ssp126)

salts_miroc_ssp126<- list()
for(j in 2:18){
  bering_model_salt = nc_open(get_salt_filepath('miroc_ssp126', years_list[j]))
  salt_output <- nc_extract(bering_model_salt, salt, 'salt')
  nc_close(bering_model_salt)
  gc()
  salts_miroc_ssp126[[paste("salt_output", j, sep = "")]] <- salt_output
}
saveRDS(salts_miroc_ssp126, file = here('data', 'salts_miroc_ssp126.rds'))
rm(salts_miroc_ssp126)

#### MIROC SSP585 ----
temps_miroc_ssp585<- list()
for(j in 2:18){
  bering_model_temp = nc_open(get_temp_filepath('miroc_ssp585', years_list[j]))
  temp_output <- nc_extract(bering_model_temp, temp, 'temp')
  nc_close(bering_model_temp)
  gc()
  temps_miroc_ssp585[[paste("temp_output", j, sep = "")]] <- temp_output
}
saveRDS(temps_miroc_ssp585, file = here('data', 'temps_miroc_ssp585.rds'))
rm(temps_miroc_ssp585)

salts_miroc_ssp585<- list()
for(j in 2:18){
  bering_model_salt = nc_open(get_salt_filepath('miroc_ssp585', years_list[j]))
  salt_output <- nc_extract(bering_model_salt, salt, 'salt')
  nc_close(bering_model_salt)
  gc()
  salts_miroc_ssp585[[paste("salt_output", j, sep = "")]] <- salt_output
}
saveRDS(salts_miroc_ssp585, file = here('data', 'salts_miroc_ssp585.rds'))
rm(salts_miroc_ssp585)

#### MIROC Historical ----
temps_miroc_historical<- list()
for(j in 1:7){
  bering_model_temp = nc_open(get_temp_filepath('miroc_historical', years_list_historical[j]))
  temp_output <- nc_extract(bering_model_temp, temp, 'temp')
  nc_close(bering_model_temp)
  gc()
  temps_miroc_historical[[paste("temp_output", j, sep = "")]] <- temp_output
}
saveRDS(temps_miroc_historical, file = here('data', 'temps_miroc_historical.rds'))
rm(temps_miroc_historical)

salts_miroc_historical<- list()
for(j in 1:7){
  bering_model_salt = nc_open(get_salt_filepath('miroc_historical', years_list_historical[j]))
  salt_output <- nc_extract(bering_model_salt, salt, 'salt')
  nc_close(bering_model_salt)
  gc()
  salts_miroc_historical[[paste("salt_output", j, sep = "")]] <- salt_output
}
saveRDS(salts_miroc_historical, file = here('data', 'salts_miroc_historical.rds'))
rm(salts_miroc_historical)












### TidyNC Method------------------------------------------------------------------------------------------------------------------------------
#### Load netcdf with tidync
library(here)
library(ncdf4)
library(tidync)
require(tidyverse)
library(data.table)
library(sf)

# download from server
url_base <- "https://data.pmel.noaa.gov/aclim/thredds/"
opendap_area  <- "dodsC/extended_grid/Bering10K_extended_grid.nc"

nc <- nc_open(paste(url_base, opendap_area, sep = ""))

# create objects for known lats and longs and xi and eta axes
lats <- ncvar_get(nc, "lat_rho")
lons <- ncvar_get(nc, "lon_rho")

nc_close(nc)

#### Hindcast ####
# Temperature
# read in files
hindcast_temp_file_list <- list.files(path = "F:/data/temperature_netcdf")

prestring <- "F:/data/temperature_netcdf/"

hindcast_temp_dat_list <- list()

for(i in hindcast_temp_file_list){
  hindcast_temp_dat_list[[i]] <- paste0(prestring, i)
  hindcast_temp_dat_list
}

hindcast_temp_df_list <- list()
for(i in hindcast_temp_dat_list){
  hindcast_temp_df_list[[i]] <- tidync(i) %>%
    hyper_tibble(select_var = "temp")
  hindcast_temp_df_list
}

hindcast_temp_dfs <- bind_rows(hindcast_temp_df_list)

# add in lat/longs matched to xi/eta 
hindcast_temp_dfs$lon <- lons[cbind(hindcast_temp_dfs$xi_rho, hindcast_temp_dfs$eta_rho)]
hindcast_temp_dfs$lat <- lats[cbind(hindcast_temp_dfs$xi_rho, hindcast_temp_dfs$eta_rho)]

# create object for time axis
hindcast_temp_dfs$DateTime <- as.POSIXct(hindcast_temp_dfs$ocean_time, origin = "1900-01-01", tz = "GMT")

saveRDS(hindcast_temp_dfs, file = "F:/data/hindcast_temp_dfs.rds")

# remove objects
rm(hindcast_temp_dat_list, 
   hindcast_temp_df_list,
   hindcast_temp_dfs)


# Salinity
# read in files
hindcast_salt_file_list <- list.files(path = "F:/data/salinity_netcdf")

prestring <- "F:/data/salinity_netcdf/"

hindcast_salt_dat_list <- list()

for(i in hindcast_salt_file_list){
  hindcast_salt_dat_list[[i]] <- paste0(prestring, i)
  hindcast_salt_dat_list
}

hindcast_salt_df_list <- list()
for(i in hindcast_salt_dat_list){
  hindcast_salt_df_list[[i]] <- tidync(i) %>%
    hyper_tibble(select_var = "salt")
  hindcast_salt_df_list
}

hindcast_salt_dfs <- bind_rows(hindcast_salt_df_list)

# add in lat/longs matched to xi/eta 
hindcast_salt_dfs$lon <- lons[cbind(hindcast_salt_dfs$xi_rho, hindcast_salt_dfs$eta_rho)]
hindcast_salt_dfs$lat <- lats[cbind(hindcast_salt_dfs$xi_rho, hindcast_salt_dfs$eta_rho)]

# create object for time axis
hindcast_salt_dfs$DateTime <- as.POSIXct(hindcast_salt_dfs$ocean_time, origin = "1900-01-01", tz = "GMT")

saveRDS(hindcast_salt_dfs, file = "F:/data/hindcast_salt_dfs.rds")

# remove objects
rm(hindcast_salt_dat_list, 
   hindcast_salt_df_list,
   hindcast_salt_dfs)

# Functions
get_values <- function(path, projection, type){
  file_list <- list.files(path = path) # read in files
  
  baseline_file_list <- file_list[str_detect(file_list, projection)]
  
  prestring <- path
  
  dat_list <- list()
  
  for(i in baseline_file_list){
    dat_list[[i]] <- paste0(prestring, i)
    dat_list
  }
  
  df_list <- list()
  for(i in dat_list){
    df_list[[i]] <- tidync(i) %>%
      hyper_filter(s_rho = index > 29) %>% # extract only the top layer
      hyper_tibble(select_var = type) # get just temperature
    df_list
  }
  
  dfs <- bind_rows(df_list)
  
  # add in lat/longs matched to xi/eta 
  dfs$lon <- lons[cbind(dfs$xi_rho, dfs$eta_rho)]
  dfs$lat <- lats[cbind(dfs$xi_rho, dfs$eta_rho)]
  
  # create object for time axis
  dfs$DateTime <- as.POSIXct(dfs$ocean_time, origin = "1900-01-01", tz = "GMT")
  return(dfs)
}

join_values <- function(ssp126_dfs, ssp585_dfs, hist_dfs) {
  ssp126_dfs <- select(ssp126_dfs,-c(xi_rho, eta_rho, s_rho, ocean_time))
  ssp585_dfs <- select(ssp585_dfs,-c(xi_rho, eta_rho, s_rho, ocean_time))
  hist_dfs <- select(hist_dfs,-c(xi_rho, eta_rho, s_rho, ocean_time))
  
  ssp126_dfs$projection <- "ssp126"
  ssp585_dfs$projection <- "ssp585"
  hist_dfs$projection <- "historical"
  
  # likely need to increase memory to do this
  dfs <- bind_rows(hist_dfs, ssp126_dfs, ssp585_dfs)
  
  # separate date column into components
  dfs$date <- as.Date(dfs$DateTime) # date in Date format
  dfs$month <- month(dfs$date) # month of year
  dfs$week <- week(dfs$date) # week of year
  dfs$year <- year(dfs$date)
  
  # remove all months aside from Feb - Aug
  months <- c(2:9)
  
  dfs_trim <- dfs %>%
    filter(., month %in% months)
  
  dfs_trim <- select(dfs_trim,-DateTime)
  
  # trim df to those lat/lons in hindcast df
  dfs_trim <- na.omit(dfs_trim)
  return(dfs_trim)
}


#### CESM simulations ####
# Temperature
cesm_temp_hist_dfs <- get_values("F:/CMIP6_relevant/cesm/temp/", "historical", "temp")
saveRDS(cesm_temp_hist_dfs, file = 'F:/data/cesm_temp_hist_dfs.rds')

cesm_temp_ssp126_dfs <- get_values("F:/CMIP6_relevant/cesm/temp/", "ssp126", "temp")
saveRDS(cesm_temp_hist_dfs, file = 'F:/data/cesm_temp_ssp126_dfs.rds')

cesm_temp_ssp585_dfs <- get_values("F:/CMIP6_relevant/cesm/temp/", "ssp585", "temp")
saveRDS(cesm_temp_hist_dfs, file = 'F:/data/cesm_temp_ssp585_dfs.rds')

cesm_temp_dfs_trim <- join_values(cesm_temp_ssp126_dfs, cesm_temp_ssp585_dfs, cesm_temp_hist_dfs)
saveRDS(cesm_temp_dfs_trim, file = 'F:/data/cesm_temp_dfs_trim.rds')

rm(cesm_temp_hist_dfs, cesm_temp_ssp126_dfs, cesm_temp_ssp585_dfs, cesm_temp_dfs_trim)

# Salinity
cesm_salt_hist_dfs <- get_values("F:/CMIP6_relevant/cesm/salt/", "historical", "salt")
saveRDS(cesm_salt_hist_dfs, file = 'F:/data/cesm_salt_hist_dfs.rds')

cesm_salt_ssp126_dfs <- get_values("F:/CMIP6_relevant/cesm/salt/", "ssp126", "salt")
saveRDS(cesm_salt_hist_dfs, file = 'F:/data/cesm_salt_ssp126_dfs.rds')

cesm_salt_ssp585_dfs <- get_values("F:/CMIP6_relevant/cesm/salt/", "ssp585", "salt")
saveRDS(cesm_salt_hist_dfs, file = 'F:/data/cesm_salt_ssp585_dfs.rds')

cesm_salt_dfs_trim <- join_values(cesm_salt_ssp126_dfs, cesm_salt_ssp585_dfs, cesm_salt_hist_dfs)
saveRDS(cesm_salt_dfs_trim, file = 'F:/data/cesm_salt_dfs_trim.rds')

rm(cesm_salt_hist_dfs, cesm_salt_ssp126_dfs, cesm_salt_ssp585_dfs, cesm_salt_dfs_trim)


#### GFDL simulations ####
# Temperature
gfdl_temp_hist_dfs <- get_values("F:/CMIP6_relevant/gfdl/temp/", "historical", "temp")
saveRDS(gfdl_temp_hist_dfs, file = 'F:/data/gfdl_temp_hist_dfs.rds')

gfdl_temp_ssp126_dfs <- get_values("F:/CMIP6_relevant/gfdl/temp/", "ssp126", "temp")
saveRDS(gfdl_temp_hist_dfs, file = 'F:/data/gfdl_temp_ssp126_dfs.rds')

gfdl_temp_ssp585_dfs <- get_values("F:/CMIP6_relevant/gfdl/temp/", "ssp585", "temp")
saveRDS(gfdl_temp_hist_dfs, file = 'F:/data/gfdl_temp_ssp585_dfs.rds')

gfdl_temp_dfs_trim <- join_values(gfdl_temp_ssp126_dfs, gfdl_temp_ssp585_dfs, gfdl_temp_hist_dfs)
saveRDS(gfdl_temp_dfs_trim, file = 'F:/data/gfdl_temp_dfs_trim.rds')

rm(gfdl_temp_hist_dfs, gfdl_temp_ssp126_dfs, gfdl_temp_ssp585_dfs, gfdl_temp_dfs_trim)

# Salinity
gfdl_salt_hist_dfs <- get_values("F:/CMIP6_relevant/gfdl/salt/", "historical", "salt")
saveRDS(gfdl_salt_hist_dfs, file = 'F:/data/gfdl_salt_hist_dfs.rds')

gfdl_salt_ssp126_dfs <- get_values("F:/CMIP6_relevant/gfdl/salt/", "ssp126", "salt")
saveRDS(gfdl_salt_hist_dfs, file = 'F:/data/gfdl_salt_ssp126_dfs.rds')

gfdl_salt_ssp585_dfs <- get_values("F:/CMIP6_relevant/gfdl/salt/", "ssp585", "salt")
saveRDS(gfdl_salt_hist_dfs, file = 'F:/data/gfdl_salt_ssp585_dfs.rds')

gfdl_salt_dfs_trim <- join_values(gfdl_salt_ssp126_dfs, gfdl_salt_ssp585_dfs, gfdl_salt_hist_dfs)
saveRDS(gfdl_salt_dfs_trim, file = 'F:/data/gfdl_salt_dfs_trim.rds')

rm(gfdl_salt_hist_dfs, gfdl_salt_ssp126_dfs, gfdl_salt_ssp585_dfs, gfdl_salt_dfs_trim)


#### MIROC simulations ####
# Temperature
miroc_temp_hist_dfs <- get_values("F:/CMIP6_relevant/miroc/temp/", "historical", "temp")
saveRDS(miroc_temp_hist_dfs, file = 'F:/data/miroc_temp_hist_dfs.rds')

miroc_temp_ssp126_dfs <- get_values("F:/CMIP6_relevant/miroc/temp/", "ssp126", "temp")
saveRDS(miroc_temp_hist_dfs, file = 'F:/data/miroc_temp_ssp126_dfs.rds')

miroc_temp_ssp585_dfs <- get_values("F:/CMIP6_relevant/miroc/temp/", "ssp585", "temp")
saveRDS(miroc_temp_hist_dfs, file = 'F:/data/miroc_temp_ssp585_dfs.rds')

miroc_temp_dfs_trim <- join_values(miroc_temp_ssp126_dfs, miroc_temp_ssp585_dfs, miroc_temp_hist_dfs)
saveRDS(miroc_temp_dfs_trim, file = 'F:/data/miroc_temp_dfs_trim.rds')

rm(miroc_temp_hist_dfs, miroc_temp_ssp126_dfs, miroc_temp_ssp585_dfs, miroc_temp_dfs_trim)

# Salinity
miroc_salt_hist_dfs <- get_values("F:/CMIP6_relevant/miroc/salt/", "historical", "salt")
saveRDS(miroc_salt_hist_dfs, file = 'F:/data/miroc_salt_hist_dfs.rds')

miroc_salt_ssp126_dfs <- get_values("F:/CMIP6_relevant/miroc/salt/", "ssp126", "salt")
saveRDS(miroc_salt_hist_dfs, file = 'F:/data/miroc_salt_ssp126_dfs.rds')

miroc_salt_ssp585_dfs <- get_values("F:/CMIP6_relevant/miroc/salt/", "ssp585", "salt")
saveRDS(miroc_salt_hist_dfs, file = 'F:/data/miroc_salt_ssp585_dfs.rds')

miroc_salt_dfs_trim <- join_values(miroc_salt_ssp126_dfs, miroc_salt_ssp585_dfs, miroc_salt_hist_dfs)
saveRDS(miroc_salt_dfs_trim, file = 'F:/data/miroc_salt_dfs_trim.rds')

rm(miroc_salt_hist_dfs, miroc_salt_ssp126_dfs, miroc_salt_ssp585_dfs, miroc_salt_dfs_trim)

