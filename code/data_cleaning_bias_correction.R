# Title: Bias correction of ROMS forecast
# Date: 1/31/2022

### Libraries and functions ----
library(spacetime)
library(tidyverse)
library(lubridate)
library(date)

### Functions ----
hindcast_mean_temp <- function(hindcast_dfs){
  
  baseline_years <- 1980:2014 # define baseline/ref years (here based on Cheng et al 2021)
  
  hindcast_dfs$year <- year(hindcast_dfs$DateTime)
  hindcast_dfs$month <- month(hindcast_dfs$DateTime)
  
  ROMS_baseline_dat <- hindcast_dfs %>% # select ref yrs from df
    filter(., year %in% baseline_years)
  
  # estimate a monthly-avg temp for each grid cell for each month
  # for the reference period
  # (so an avg temp for each month at each grid cell averaged across 1980 - 2014)
  ROMS_baseline_dat_mo <- ROMS_baseline_dat %>%
    mutate(lon = lon,
           lat = lat) %>%
    group_by(month, lat, lon) %>%
    summarize(mo_baseline_value = mean(temp))
  return(ROMS_baseline_dat_mo)
}
hindcast_mean_salt <- function(hindcast_dfs){
  
  baseline_years <- 1980:2014 # define baseline/ref years (here based on Cheng et al 2021)
  
  hindcast_dfs$year <- year(hindcast_dfs$DateTime)
  hindcast_dfs$month <- month(hindcast_dfs$DateTime)
  
  ROMS_baseline_dat <- hindcast_dfs %>% # select ref yrs from df
    filter(., year %in% baseline_years)
  
  # estimate a monthly-avg temp for each grid cell for each month
  # for the reference period
  # (so an avg temp for each month at each grid cell averaged across 1980 - 2014)
  ROMS_baseline_dat_mo <- ROMS_baseline_dat %>%
    mutate(lon = lon,
           lat = lat) %>%
    group_by(month, lat, lon) %>%
    summarize(mo_baseline_value = mean(salt))
  return(ROMS_baseline_dat_mo)
}
baseline_mean_temp <- function(baseline_dfs){
  baseline_years <- 1980:2014
  baseline_dat <- baseline_dfs %>% # select ref yrs from df
    filter(., year %in% baseline_years)
  
  # estimate a monthly-avg temp for each grid cell for each month
  # for the ref period
  # (so an avg temp for each month at each grid cell averaged across 1980 - 2014)
  baseline_dat_mo <- baseline_dat %>%
    group_by(month, lat, lon) %>%
    summarize(mean_proj_baseline = mean(temp))
  return(baseline_dat_mo)
}
baseline_mean_salt <- function(baseline_dfs){
  baseline_years <- 1980:2014
  baseline_dat <- baseline_dfs %>% # select ref yrs from df
    filter(., year %in% baseline_years)
  
  # estimate a monthly-avg temp for each grid cell for each month
  # for the ref period
  # (so an avg temp for each month at each grid cell averaged across 1980 - 2014)
  baseline_dat_mo <- baseline_dat %>%
    group_by(month, lat, lon) %>%
    summarize(mean_proj_baseline = mean(salt))
  return(baseline_dat_mo)
}
forecast_deltas_temp <- function(forecast_dfs, projection_years,  
                                 baseline_means, hindcast_means){
  proj_dat <- forecast_dfs %>%
    filter(., year %in% projection_years)
  
  proj_dat <- proj_dat %>%
    group_by(projection, year, month, lat, lon) %>%
    summarise(mo_avg_proj = mean(temp))
  
  #4 calculate deltas (difference btw raw projected temp and mean proj temp across
  # ref period)
  
  # combine the monthly means for historical period and projected df into one df
  delta_dat <- merge(proj_dat, baseline_means,
                     by = c("lat", "lon", "month"))
  
  delta_dat <- delta_dat %>%
    mutate(delta = (mo_avg_proj - mean_proj_baseline))
  
  #5 add deltas mean of the hindcast during the reference years (step 1)
  bcs <- merge(hindcast_means, delta_dat,
               by = c("lat", "lon", "month"))
  
  bcs <- bcs %>%
    mutate(bc = delta + mo_baseline_value)
  return(bcs)
}
forecast_deltas_salt <- function(forecast_dfs, projection_years,  
                                 baseline_means, hindcast_means){
  proj_dat <- forecast_dfs %>%
    filter(., year %in% projection_years)
  
  proj_dat <- proj_dat %>%
    group_by(projection, year, month, lat, lon) %>%
    summarise(mo_avg_proj = mean(salt))
  
  #4 calculate deltas (difference btw raw projected temp and mean proj temp across
  # ref period)
  
  # combine the monthly means for historical period and projected df into one df
  delta_dat <- merge(proj_dat, baseline_means,
                     by = c("lat", "lon", "month"))
  
  delta_dat <- delta_dat %>%
    mutate(delta = (mo_avg_proj - mean_proj_baseline))
  
  #5 add deltas mean of the hindcast during the reference years (step 1)
  bcs <- merge(hindcast_means, delta_dat,
               by = c("lat", "lon", "month"))
  
  bcs <- bcs %>%
    mutate(bc = delta + mo_baseline_value)
  return(bcs)
}

# To run the forecast function, memory limit likely needs to be increased
# The steps for bias correction are outlined in the functions (steps 1-5)

memory.limit(300000)

### CESM ----
# Temperature
hindcast_temp_dfs <- readRDS('F:/data/hindcast_temp_dfs.rds')
cesm_temp_dfs <- readRDS('F:/data/cesm_temp_dfs_trim.rds')

cesm_hindcast_temps <- hindcast_mean_temp(hindcast_temp_dfs)
cesm_baseline_temps <- baseline_mean_temp(cesm_temp_dfs)
cesm_forecast_temp1 <- forecast_deltas_temp(cesm_temp_dfs, 2015:2039,
                                            cesm_baseline_temps, 
                                            cesm_hindcast_temps)
saveRDS(cesm_forecast_temp1, 'F:/data/cesm_forecast_temp1.rds')
rm(cesm_forecast_temp1)
cesm_forecast_temp2 <- forecast_deltas_temp(cesm_temp_dfs, 2040:2069,
                                            cesm_baseline_temps, 
                                            cesm_hindcast_temps)
saveRDS(cesm_forecast_temp2, 'F:/data/cesm_forecast_temp2.rds')
rm(cesm_forecast_temp2)
cesm_forecast_temp3 <- forecast_deltas_temp(cesm_temp_dfs, 2070:2099,
                                            cesm_baseline_temps, 
                                            cesm_hindcast_temps)
saveRDS(cesm_forecast_temp3, 'F:/data/cesm_forecast_temp3.rds')
rm(cesm_forecast_temp3, hindcast_temp_dfs, cesm_temp_dfs)
gc()

# Salinity
hindcast_salt_dfs <- readRDS('F:/data/hindcast_salt_dfs.rds')
cesm_salt_dfs <- readRDS('F:/data/cesm_salt_dfs_trim.rds')

cesm_hindcast_salts <- hindcast_mean_salt(hindcast_salt_dfs)
cesm_baseline_salts <- baseline_mean_salt(cesm_salt_dfs)
cesm_forecast_salt1 <- forecast_deltas_salt(cesm_salt_dfs, 2015:2039,
                                            cesm_baseline_salts, 
                                            cesm_hindcast_salts)
saveRDS(cesm_forecast_salt1, 'F:/data/cesm_forecast_salt1.rds')
rm(cesm_forecast_salt1)
cesm_forecast_salt2 <- forecast_deltas_salt(cesm_salt_dfs, 2040:2069,
                                            cesm_baseline_salts, 
                                            cesm_hindcast_salts)
saveRDS(cesm_forecast_salt2, 'F:/data/cesm_forecast_salt2.rds')
rm(cesm_forecast_salt2)
cesm_forecast_salt3 <- forecast_deltas_salt(cesm_salt_dfs, 2070:2099,
                                            cesm_baseline_salts, 
                                            cesm_hindcast_salts)
saveRDS(cesm_forecast_salt3, 'F:/data/cesm_forecast_salt3.rds')
rm(cesm_forecast_salt3, hindcast_salt_dfs, cesm_salt_dfs)
gc()


### GFDL ----
# Temperature
hindcast_temp_dfs <- readRDS('F:/data/hindcast_temp_dfs.rds')
gfdl_temp_dfs <- readRDS('F:/data/gfdl_temp_dfs_trim.rds')

gfdl_hindcast_temps <- hindcast_mean_temp(hindcast_temp_dfs)
gfdl_baseline_temps <- baseline_mean_temp(gfdl_temp_dfs)
gfdl_forecast_temp1 <- forecast_deltas_temp(gfdl_temp_dfs, 2015:2039,
                                            gfdl_baseline_temps, 
                                            gfdl_hindcast_temps)
saveRDS(gfdl_forecast_temp1, 'F:/data/gfdl_forecast_temp1.rds')
rm(gfdl_forecast_temp1)
gfdl_forecast_temp2 <- forecast_deltas_temp(gfdl_temp_dfs, 2040:2069,
                                            gfdl_baseline_temps, 
                                            gfdl_hindcast_temps)
saveRDS(gfdl_forecast_temp2, 'F:/data/gfdl_forecast_temp2.rds')
rm(gfdl_forecast_temp2)
gfdl_forecast_temp3 <- forecast_deltas_temp(gfdl_temp_dfs, 2070:2099,
                                            gfdl_baseline_temps, 
                                            gfdl_hindcast_temps)
saveRDS(gfdl_forecast_temp3, 'F:/data/gfdl_forecast_temp3.rds')
rm(gfdl_forecast_temp3, hindcast_temp_dfs, gfdl_temp_dfs)
gc()

# Salinity
hindcast_salt_dfs <- readRDS('F:/data/hindcast_salt_dfs.rds')
gfdl_salt_dfs <- readRDS('F:/data/gfdl_salt_dfs_trim.rds')

gfdl_hindcast_salts <- hindcast_mean_salt(hindcast_salt_dfs)
gfdl_baseline_salts <- baseline_mean_salt(gfdl_salt_dfs)
gfdl_forecast_salt1 <- forecast_deltas_salt(gfdl_salt_dfs, 2015:2039,
                                            gfdl_baseline_salts, 
                                            gfdl_hindcast_salts)
saveRDS(gfdl_forecast_salt1, 'F:/data/gfdl_forecast_salt1.rds')
rm(gfdl_forecast_salt1)
gfdl_forecast_salt2 <- forecast_deltas_salt(gfdl_salt_dfs, 2040:2069,
                                            gfdl_baseline_salts, 
                                            gfdl_hindcast_salts)
saveRDS(gfdl_forecast_salt2, 'F:/data/gfdl_forecast_salt2.rds')
rm(gfdl_forecast_salt2)
gfdl_forecast_salt3 <- forecast_deltas_salt(gfdl_salt_dfs, 2070:2099,
                                            gfdl_baseline_salts, 
                                            gfdl_hindcast_salts)
saveRDS(gfdl_forecast_salt3, 'F:/data/gfdl_forecast_salt3.rds')
rm(gfdl_forecast_salt3, hindcast_salt_dfs, gfdl_salt_dfs)
gc()


### MIROC ----
# Temperature
hindcast_temp_dfs <- readRDS('F:/data/hindcast_temp_dfs.rds')
miroc_temp_dfs <- readRDS('F:/data/miroc_temp_dfs_trim.rds')

miroc_hindcast_temps <- hindcast_mean_temp(hindcast_temp_dfs)
miroc_baseline_temps <- baseline_mean_temp(miroc_temp_dfs)
miroc_forecast_temp1 <- forecast_deltas_temp(miroc_temp_dfs, 2015:2039,
                                            miroc_baseline_temps, 
                                            miroc_hindcast_temps)
saveRDS(miroc_forecast_temp1, 'F:/data/miroc_forecast_temp1.rds')
rm(miroc_forecast_temp1)
miroc_forecast_temp2 <- forecast_deltas_temp(miroc_temp_dfs, 2040:2069,
                                            miroc_baseline_temps, 
                                            miroc_hindcast_temps)
saveRDS(miroc_forecast_temp2, 'F:/data/miroc_forecast_temp2.rds')
rm(miroc_forecast_temp2)
miroc_forecast_temp3 <- forecast_deltas_temp(miroc_temp_dfs, 2070:2099,
                                            miroc_baseline_temps, 
                                            miroc_hindcast_temps)
saveRDS(miroc_forecast_temp3, 'F:/data/miroc_forecast_temp3.rds')
rm(miroc_forecast_temp3, hindcast_temp_dfs, miroc_temp_dfs)
gc()

# Salinity
hindcast_salt_dfs <- readRDS('F:/data/hindcast_salt_dfs.rds')
miroc_salt_dfs <- readRDS('F:/data/miroc_salt_dfs_trim.rds')

miroc_hindcast_salts <- hindcast_mean_salt(hindcast_salt_dfs)
miroc_baseline_salts <- baseline_mean_salt(miroc_salt_dfs)
miroc_forecast_salt1 <- forecast_deltas_salt(miroc_salt_dfs, 2015:2039,
                                            miroc_baseline_salts, 
                                            miroc_hindcast_salts)
saveRDS(miroc_forecast_salt1, 'F:/data/miroc_forecast_salt1.rds')
rm(miroc_forecast_salt1)
miroc_forecast_salt2 <- forecast_deltas_salt(miroc_salt_dfs, 2040:2069,
                                            miroc_baseline_salts, 
                                            miroc_hindcast_salts)
saveRDS(miroc_forecast_salt2, 'F:/data/miroc_forecast_salt2.rds')
rm(miroc_forecast_salt2)
miroc_forecast_salt3 <- forecast_deltas_salt(miroc_salt_dfs, 2070:2099,
                                            miroc_baseline_salts, 
                                            miroc_hindcast_salts)
saveRDS(miroc_forecast_salt3, 'F:/data/miroc_forecast_salt3.rds')
rm(miroc_forecast_salt3, hindcast_salt_dfs, miroc_salt_dfs)
gc()