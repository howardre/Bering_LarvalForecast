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

### Create vectors with variables and save as list
# May need to remove from environment in order to get through them all
# Depends on amount of memory available
# CESM SSP126
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

# CESM SSP585
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

# CESM Historical
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

# GFDL SSP126
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

# GFDL SSP585
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

# GFDL Historical
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

# MIROC SSP126
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

# MIROC SSP585
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

# MIROC Historical
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












####--------------------------------------------TidyNC Method----------------------------------------------------------------------------------
### Load netcdf with tidync
# read in netcdf files of projected temp from Bering 10k ROMS

library(here)
library(ncdf4)
library(tidync)
require(tidyverse)
library(data.table)
library(sf)

# set up lat/lons from area grid file

# download from server
url_base <- "https://data.pmel.noaa.gov/aclim/thredds/"
opendap_area  <- "dodsC/extended_grid/Bering10K_extended_grid.nc"

nc <- nc_open(paste(url_base, opendap_area, sep = ""))

# create objects for known lats and longs and xi and eta axes
lats <- ncvar_get(nc, "lat_rho")
lons <- ncvar_get(nc, "lon_rho")

nc_close(nc)

# CESM simulations ####

# read in files
cesm_file_list <- list.files(path = "D:/CMIP6_relevant/cesm/temp/")

# historical baseline period ####
cesm_historical_baseline_file_list <- cesm_file_list[str_detect(cesm_file_list, "historical")]

prestring <- "D:/CMIP6_relevant/cesm/temp/"

cesm_hist_dat_list <- list()

for(i in cesm_historical_baseline_file_list){
  cesm_hist_dat_list[[i]] <- paste0(prestring, i)
  cesm_hist_dat_list
}

cesm_hist_df_list <- list()
for(i in cesm_hist_dat_list){
  cesm_hist_df_list[[i]] <- tidync(i) %>%
    hyper_filter(s_rho = index > 29) %>%
    hyper_tibble(select_var = "temp")
  cesm_hist_df_list
}

cesm_hist_dfs <- bind_rows(cesm_hist_df_list)

# add in lat/longs matched to xi/eta 
cesm_hist_dfs$Lon <- lons[cbind(cesm_hist_dfs$xi_rho, cesm_hist_dfs$eta_rho)]
cesm_hist_dfs$Lat <- lats[cbind(cesm_hist_dfs$xi_rho, cesm_hist_dfs$eta_rho)]

# create object for time axis
cesm_hist_dfs$DateTime <- as.POSIXct(cesm_hist_dfs$ocean_time, origin = "1900-01-01", tz = "GMT")


# ssp 126 projection ####

cesm_ssp126_file_list <- cesm_file_list[str_detect(cesm_file_list, "ssp126")]

prestring <- "D:/CMIP6_relevant/cesm/temp/"

cesm_ssp126_dat_list <- list()

for(i in cesm_ssp126_file_list){
  cesm_ssp126_dat_list[[i]] <- paste0(prestring, i)
  cesm_ssp126_dat_list
}

cesm_ssp126_df_list <- list()
for(i in cesm_ssp126_dat_list){
  cesm_ssp126_df_list[[i]] <- tidync(i) %>% 
    hyper_filter(s_rho = index > 29) %>%
    hyper_tibble(select_var = "temp")
  cesm_ssp126_df_list
}

cesm_ssp126_dfs <- bind_rows(cesm_ssp126_df_list)

# add in lat/longs matched to xi/eta 
cesm_ssp126_dfs$Lon <- lons[cbind(cesm_ssp126_dfs$xi_rho, cesm_ssp126_dfs$eta_rho)]
cesm_ssp126_dfs$Lat <- lats[cbind(cesm_ssp126_dfs$xi_rho, cesm_ssp126_dfs$eta_rho)]

# create object for time axis
cesm_ssp126_dfs$DateTime <- as.POSIXct(cesm_ssp126_dfs$ocean_time, origin = "1900-01-01", tz = "GMT")

saveRDS(cesm_ssp126_dfs, file = here('data', 'cesm_ssp126_dfs.rds'))

# ssp585 projection ####

cesm_ssp585_file_list <- cesm_file_list[str_detect(cesm_file_list, "ssp585")]

prestring <- "D:/CMIP6_relevant/cesm/temp/"

cesm_ssp585_dat_list <- list()
for(i in cesm_ssp585_file_list){
  cesm_ssp585_dat_list[[i]] <- paste0(prestring, i)
  cesm_ssp585_dat_list
}

cesm_ssp585_df_list <- list()
for(i in cesm_ssp585_dat_list){
  cesm_ssp585_df_list[[i]] <- tidync(i) %>% 
    hyper_filter(s_rho = index > 29) %>%
    hyper_tibble(select_var = "temp")
  cesm_ssp585_df_list
}

cesm_ssp585_dfs <- bind_rows(cesm_ssp585_df_list)

# add in lat/longs matched to xi/eta 
cesm_ssp585_dfs$Lon <- lons[cbind(cesm_ssp585_dfs$xi_rho, cesm_ssp585_dfs$eta_rho)]
cesm_ssp585_dfs$Lat <- lats[cbind(cesm_ssp585_dfs$xi_rho, cesm_ssp585_dfs$eta_rho)]

# create object for time axis
cesm_ssp585_dfs$DateTime <- as.POSIXct(cesm_ssp585_dfs$ocean_time, origin = "1900-01-01", tz = "GMT")

cesm_ssp126_dfs <- select(cesm_ssp126_dfs, -c(xi_rho, eta_rho, s_rho, ocean_time, station))
cesm_ssp585_dfs <- select(cesm_ssp585_dfs, -c(xi_rho, eta_rho, s_rho, ocean_time, station))
cesm_hist_dfs <- select(cesm_hist_dfs, -c(xi_rho, eta_rho, s_rho, ocean_time, station))

# join together
cesm_ssp126_dfs$projection <- "ssp126"
cesm_ssp585_dfs$projection <- "ssp585"
cesm_hist_dfs$projection <- "historical"


cesm_dfs <- bind_rows(cesm_hist_dfs, cesm_ssp126_dfs, cesm_ssp585_dfs) 

# separate date column into components
cesm_dfs$date <- as.Date(cesm_dfs$DateTime) # date in Date format
cesm_dfs$month <- month(cesm_dfs$date) # month of year
cesm_dfs$week <- week(cesm_dfs$date) # week of year
cesm_dfs$year <- year(cesm_dfs$date)

# remove all months aside from Feb - Aug
months <- c(2:9)

cesm_dfs_trim <- cesm_dfs %>%
  filter(., month %in% months)

cesm_dfs_trim <- select(cesm_dfs_trim, -DateTime)

# trim df to those lat/lons in hindcast df 

# summarize by lat/lon and convert to sf object
cesm_dfs_trim_sum <- cesm_dfs_trim %>%
  group_by(Lat, Lon) %>%
  summarize(mean_temp = mean(temp)) %>%
  mutate(latitude = Lat,
         long_not_360 = case_when(
           Lon >= 180 ~ Lon - 360,
           Lon < 180 ~ Lon)) %>%
  st_as_sf(coords = c("long_not_360", "latitude"), crs = 4326)

cesm_dfs_trim_sum <- cesm_dfs_trim_sum %>%
  rename(latitude = Lat,
         longitude = Lon)

# make a summary object of the hindcast data for intersecting the lat/lons
ROMS_hindcast_dat_sum <- ROMS_hindcast_dat %>%
  group_by(latitude, longitude) %>%
  summarise(mean_temp = mean(temp)) %>%
  mutate(long_not_360 = case_when(
    longitude >= 180 ~ longitude - 360,
    longitude < 180 ~ longitude)) %>%
  st_as_sf(coords = c("long_not_360", "latitude"), crs = 4326)

dat_ints <- st_intersection(cesm_dfs_trim_sum, ROMS_hindcast_dat_sum)

cesm_dat_trim <- cesm_dfs_trim %>% 
  filter(., Lon %in% dat_ints$longitude) %>%
  filter(., Lat %in% dat_ints$latitude)

fwrite(cesm_dat_trim, "./data/cesm_dat_trim.csv")


# plots ####

# yearly plot

# summarize by year
cesm_dat_trim_sum <- cesm_dat_trim %>%
  filter(., projection != "historical") %>%
  group_by(projection, Lat, Lon, year) %>%
  summarise(mean_temp = mean(temp))

# convert to sf object
cesm_dat_trim_sum_sf <- cesm_dat_trim_sum %>% 
  mutate(latitude = Lat,
         long_not_360 = case_when(
           Lon >= 180 ~ Lon - 360,
           Lon < 180 ~ Lon)) %>%
  st_as_sf(coords = c("long_not_360", "latitude"), crs = 4326)


cesm_yr_plot_func <- function(x){
  
  new_dat <- cesm_dfs_trim_sum_sf %>% filter(., year == x)
  
  plot <- 
    ggplot() +
    geom_sf(data = new_dat, aes(color = mean_temp))  +
    geom_sf(data = world_map_data, fill = "grey", lwd = 0) +
    facet_wrap(~ projection) +
    coord_sf(crs = 3338) +
    scale_color_viridis_c() +
    scale_x_continuous(
      breaks = c(-175, -170, -165, -160),
      labels = c("-175˚", "-170˚", "-165˚", "-160˚"),
      name = "Longitude",
      limits = c(-1400000, -150000)
    ) +
    scale_y_continuous(
      breaks = c(55, 60),
      limits = c(470000, 1900000),
      name = "Latitude",
    ) +
    labs(colour = "bottom temp C") +
    theme_bw() +
    theme(
      axis.text = element_text(size = 12),	
      axis.title = element_text(size = 14),
      legend.title.align=0.5)
  
  plot
  
}

years <- unique(cesm_dfs_trim_sum_sf$year)

yr_plot_list <- lapply(years, cesm_yr_plot_func)

yr_name_func_year <- function(x){
  year_name <- paste0(x, "_cesm_yr_btemp.png")
}

names_year <- sapply(years, yr_name_func_year)

yr_name_func_file <- function(x){
  year_file_name <- paste0("/Users/jenniferbigman/My Drive/NOAA AFSC Postdoc/Pcod Bering Sea Habitat Suitability/Pcod-Bering-Sea/output/plots/yearly plots/projected temps/cesm plots/", x)
}

yearly_names <- sapply(names_year, yr_name_func_file)

plot_list <- mapply(ggsave_func, x = yr_plot_list, y = yearly_names)


# plot: summary by decade

cesm_dat_trim_sum_decade_sf <- cesm_dat_trim_sum_sf %>%
  mutate(decade = case_when(
    between(year, 2010, 2019) ~ "2010s",
    between(year, 2020, 2029) ~ "2020s",
    between(year, 2030, 2039) ~ "2030s",
    between(year, 2040, 2049) ~ "2040s",
    between(year, 2050, 2059) ~ "2050s",
    between(year, 2060, 2069) ~ "2060s",
    between(year, 2070, 2079) ~ "2070s",
    between(year, 2080, 2089) ~ "2080s",
    between(year, 2090, 2099) ~ "2090s"))

plot_cesm_decade <- 
  ggplot() +
  geom_sf(data = cesm_dat_trim_sum_decade_sf, 
          aes(color = mean_temp))  +
  geom_sf(data = world_map_data, fill = "grey", lwd = 0) +
  coord_sf(crs = 3338) +
  facet_grid(projection ~ decade) +
  scale_color_viridis_c() +
  scale_x_continuous(
    breaks = c(-175, -170, -165, -160),
    labels = c("-175˚", "-170˚", "-165˚", "-160˚"),
    name = "Longitude",
    limits = c(-1400000, -150000)
  ) +
  scale_y_continuous(
    breaks = c(55, 60),
    limits = c(470000, 1900000),
    name = "Latitude",
  ) +
  labs(colour = "bottom temp C") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    strip.background = element_blank(),
    axis.text = element_text(size = 12),	
    axis.title = element_text(size = 14),
    legend.title = element_text(hjust = 0.5, size = 12),
    plot.title = element_text(size = 18, face = "bold")
  )

ggsave("./output/plots/plot_cesm_decade.png",
       plot_cesm_decade,
       width = 15, height = 10, units = "in")



