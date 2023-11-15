# Title: Bias correction of ROMS forecast
# Date: 1/31/2022

### Libraries and functions ----
library(spacetime)
library(tidyverse)
library(lubridate)
library(date)
library(ggplot2)

### Functions ----
# This calculates the mean for the hindcast during the baseline years
# Jenny's paper shows that the set of years selected really doens't matter
# Part 1 of the calculation for temperature and salinity
hindcast_mean_temp <- function(hindcast_dfs){
  
  baseline_years <- 1985:2014 # define baseline/ref years (here based on Cheng et al 2021)
  
  hindcast_dfs$year <- year(hindcast_dfs$DateTime)
  hindcast_dfs$month <- month(hindcast_dfs$DateTime)
  
  ROMS_baseline_dat <- hindcast_dfs %>% # select ref yrs from df
    filter(., year %in% baseline_years)
  
  # estimate a monthly-avg temp for each grid cell for each month for the reference period
  # (so an avg temp for each month at each grid cell averaged across 1985 - 2014)
  ROMS_baseline_dat_mo <- ROMS_baseline_dat %>%
    mutate(lon = lon,
           lat = lat) %>%
    group_by(month, lat, lon) %>%
    summarize(mo_baseline_value = mean(temp))
  return(ROMS_baseline_dat_mo)
}

hindcast_mean_salt <- function(hindcast_dfs){
  
  baseline_years <- 1985:2014 # define baseline/ref years (here based on Cheng et al 2021)
  
  hindcast_dfs$year <- year(hindcast_dfs$DateTime)
  hindcast_dfs$month <- month(hindcast_dfs$DateTime)
  
  ROMS_baseline_dat <- hindcast_dfs %>% # select ref yrs from df
    filter(., year %in% baseline_years)
  
  # estimate a monthly-avg temp for each grid cell for each month
  # for the reference period
  # (so an avg temp for each month at each grid cell averaged across 1985 - 2014)
  ROMS_baseline_dat_mo <- ROMS_baseline_dat %>%
    mutate(lon = lon,
           lat = lat) %>%
    group_by(month, lat, lon) %>%
    summarize(mo_baseline_value = mean(salt))
  return(ROMS_baseline_dat_mo)
}

# These two functions are for the projections (CESM, GFDL, MIROC)
# Part 2 of the equation
baseline_mean_temp <- function(baseline_dfs){
  baseline_years <- 1985:2014
  baseline_dat <- baseline_dfs %>% # select ref yrs from df
    filter(., year %in% baseline_years)
  
  # estimate a monthly-avg temp for each grid cell for each month
  # for the ref period
  # (so an avg temp for each month at each grid cell averaged across 1985 - 2014)
  baseline_dat_mo <- baseline_dat %>%
    group_by(month, lat, lon) %>%
    summarize(mean_proj_baseline = mean(temp))
  return(baseline_dat_mo)
}

baseline_mean_salt <- function(baseline_dfs){
  baseline_years <- 1985:2014
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

# One function for parts 3-5
# Part 3 of the equation
# This calculates the average salinity and temperature per grid cell
forecast_deltas_temp <- function(forecast_dfs, projection_years,  
                                 baseline_means, hindcast_means){
  proj_dat <- forecast_dfs %>%
    filter(., year %in% projection_years)
  
  proj_dat <- proj_dat %>%
    group_by(projection, year, month, lat, lon) %>%
    summarise(mo_avg_proj = mean(temp))
  
  # Part 4 calculate deltas 
  # This is the difference between the raw projected temp/salt and mean projected temp/salt across ref period
  
  # combine the monthly means for historical period and projected df into one df
  delta_dat <- merge(proj_dat, baseline_means,
                     by = c("lat", "lon", "month"))
  
  delta_dat <- delta_dat %>%
    mutate(delta = (mo_avg_proj - mean_proj_baseline))
  
  # Part 5 is to add deltas mean of the hindcast during the reference years (step 1)
  bcs <- merge(hindcast_means, delta_dat,
               by = c("lat", "lon", "month"))
  
  bcs <- bcs %>%
    mutate(bc = delta + mo_baseline_value)
  return(bcs)
}


forecast_deltas_salt <- function(forecast_dfs, projection_years,  
                                 baseline_means, hindcast_means){
  # Part 3 of the equation
  # This calculates the average salinity and temperature per grid cell
  proj_dat <- forecast_dfs %>%
    filter(., year %in% projection_years)
  
  proj_dat <- proj_dat %>%
    group_by(projection, year, month, lat, lon) %>%
    summarise(mo_avg_proj = mean(salt))
  
  # Part 4 calculate deltas 
  # This is the difference between the raw projected temp/salt and mean projected temp/salt across ref period
  
  # combine the monthly means for historical period and projected df into one df
  delta_dat <- merge(proj_dat, baseline_means,
                     by = c("lat", "lon", "month"))
  
  delta_dat <- delta_dat %>%
    mutate(delta = (mo_avg_proj - mean_proj_baseline))
  
  # Part 5 is to add deltas mean of the hindcast during the reference years (step 1)
  bcs <- merge(hindcast_means, delta_dat,
               by = c("lat", "lon", "month"))
  
  bcs <- bcs %>%
    mutate(bc = delta + mo_baseline_value)
  return(bcs)
}

# The steps for bias correction are outlined in the functions (steps 1-5)

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

### Figures ----
summarise_years <- function(df){
  temp_sum_y <- df %>% 
    filter(lon <= 170 & lon >= 165,
           lat >= 56 & lat <= 58,
           month == 4) %>%
    group_by(projection, lat, lon, year) %>%
    summarise(mean_temp = mean(bc), .groups = 'keep')
  
  sum_year <- temp_sum_y %>%
    group_by(projection, year) %>%
    mutate(value = mean(mean_temp),
           min = min(mean_temp),
           max = max(mean_temp))
  return(sum_year)
}

summarise_years_hist <- function(df){
  temp_sum_y <- df %>%
    filter(lon <= 170 & lon >= 165,
           lat >= 56 & lat <= 58,
           month == 4) %>%
    group_by(projection, lat, lon, year) %>%
    summarise(mean_temp = mean(temp), .groups = 'keep')
  
  sum_year <- temp_sum_y %>%
    group_by(projection, year) %>%
    mutate(value = mean(mean_temp),
           min = min(mean_temp),
           max = max(mean_temp))
}

summarise_years_hist_salt <- function(df){
  salt_sum_y <- df %>%
    filter(lon <= 170 & lon >= 165,
           lat >= 56 & lat <= 58,
           month == 4) %>%
    group_by(projection, lat, lon, year) %>%
    summarise(mean_salt = mean(salt), .groups = 'keep')
  
  sum_year <- salt_sum_y %>%
    group_by(projection, year) %>%
    mutate(value = mean(mean_salt),
           min = min(mean_salt),
           max = max(mean_salt))
}

# CESM
cesm_temp_dfs1 <- readRDS('F:/data/cesm_forecast_temp1.rds')
cesm_temp_dfs2 <- readRDS('F:/data/cesm_forecast_temp2.rds')
cesm_temp_dfs3 <- readRDS('F:/data/cesm_forecast_temp3.rds')
cesm_temp_hist <- filter(readRDS('F:/data/cesm_temp_dfs_trim.rds'),
                         projection == 'historical')

cesm_temp_sum1 <- summarise_years(cesm_temp_dfs1)
cesm_temp_sum2 <- summarise_years(cesm_temp_dfs2)
cesm_temp_sum3 <- summarise_years(cesm_temp_dfs3)
cesm_temp_hist_sum <- summarise_years_hist(cesm_temp_hist)

cesm_dfs <- bind_rows(cesm_temp_sum1, 
                      cesm_temp_sum2,
                      cesm_temp_sum3,
                      cesm_temp_hist_sum)

ggplot(cesm_dfs, aes(x = year,
                     y = value,
                     group = projection,
                     fill = projection)) +
  geom_line(aes(group = projection,
                color = projection),
            show.legend = FALSE) +
  geom_ribbon(aes(group = projection,
                  fill = projection,
                  ymin = min,
                  ymax = max),
              alpha = 0.2) +
  labs(title = "CESM Projections",
       y = "Temperature (\u00B0C)",
       x = "Year",
       fill = "Projection") +
  scale_fill_manual(values = c("goldenrod3", "darkslateblue", "coral2"),
                    labels = c('Historical', "SSP1-2.6", "SSP5-8.5")) +
  scale_color_manual(values = c("goldenrod3", "darkslateblue", "coral2")) +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 25),
        plot.title = element_text(size = 40, hjust = 0.5),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 28),
        text = element_text(family = "serif"))
dev.copy(jpeg, 'F:/results/cesm_roms_temps.jpg', 
         height = 10, width = 20, units = 'in', res = 200)
dev.off()

# GFDL
gfdl_temp_dfs1 <- readRDS('F:/data/gfdl_forecast_temp1.rds')
gfdl_temp_dfs2 <- readRDS('F:/data/gfdl_forecast_temp2.rds')
gfdl_temp_dfs3 <- readRDS('F:/data/gfdl_forecast_temp3.rds')
gfdl_temp_hist <- filter(readRDS('F:/data/gfdl_temp_dfs_trim.rds'),
                         projection == 'historical')

gfdl_temp_sum1 <- summarise_years(gfdl_temp_dfs1)
gfdl_temp_sum2 <- summarise_years(gfdl_temp_dfs2)
gfdl_temp_sum3 <- summarise_years(gfdl_temp_dfs3)
gfdl_temp_hist_sum <- summarise_years_hist(gfdl_temp_hist)

gfdl_dfs <- bind_rows(gfdl_temp_sum1, 
                      gfdl_temp_sum2,
                      gfdl_temp_sum3,
                      gfdl_temp_hist_sum)

ggplot(gfdl_dfs, aes(x = year,
                     y = value,
                     group = projection,
                     fill = projection)) +
  geom_line(aes(group = projection,
                color = projection),
            show.legend = FALSE) +
  geom_ribbon(aes(group = projection,
                  fill = projection,
                  ymin = min,
                  ymax = max),
              alpha = 0.2) +
  labs(title = "GFDL Projections",
       y = "Temperature (\u00B0C)",
       x = "Year",
       fill = "Projection") +
  scale_fill_manual(values = c("goldenrod3", "darkslateblue", "coral2"),
                    labels = c('Historical', "SSP1-2.6", "SSP5-8.5")) +
  scale_color_manual(values = c("goldenrod3", "darkslateblue", "coral2")) +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 25),
        plot.title = element_text(size = 40, hjust = 0.5),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 28),
        text = element_text(family = "serif"))
dev.copy(jpeg, 'F:/results/gfdl_roms_temps.jpg', 
         height = 10, width = 20, units = 'in', res = 200)
dev.off()

# MIROC
miroc_temp_dfs1 <- readRDS('F:/data/miroc_forecast_temp1.rds')
miroc_temp_dfs2 <- readRDS('F:/data/miroc_forecast_temp2.rds')
miroc_temp_dfs3 <- readRDS('F:/data/miroc_forecast_temp3.rds')
miroc_temp_hist <- filter(readRDS('F:/data/miroc_temp_dfs_trim.rds'),
                          projection == 'historical')

miroc_temp_sum1 <- summarise_years(miroc_temp_dfs1)
miroc_temp_sum2 <- summarise_years(miroc_temp_dfs2)
miroc_temp_sum3 <- summarise_years(miroc_temp_dfs3)
miroc_temp_hist_sum <- summarise_years_hist(miroc_temp_hist)

miroc_dfs <- bind_rows(miroc_temp_sum1, 
                       miroc_temp_sum2,
                       miroc_temp_sum3,
                       miroc_temp_hist_sum)

ggplot(miroc_dfs, aes(x = year,
                      y = value,
                      group = projection,
                      fill = projection)) +
  geom_line(aes(group = projection,
                color = projection),
            show.legend = FALSE) +
  geom_ribbon(aes(group = projection,
                  fill = projection,
                  ymin = min,
                  ymax = max),
              alpha = 0.2) +
  labs(title = "MIROC Projections",
       y = "Temperature (\u00B0C)",
       x = "Year",
       fill = "Projection") +
  scale_fill_manual(values = c("goldenrod3", "darkslateblue", "coral2"),
                    labels = c('Historical', "SSP1-2.6", "SSP5-8.5")) +
  scale_color_manual(values = c("goldenrod3", "darkslateblue", "coral2")) +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 25),
        plot.title = element_text(size = 40, hjust = 0.5),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 28),
        text = element_text(family = "serif"))
dev.copy(jpeg, 'F:/results/miroc_roms_temps.jpg', 
         height = 10, width = 20, units = 'in', res = 200)
dev.off()

## Salinity
# CESM
cesm_salt_dfs1 <- readRDS('F:/data/cesm_forecast_salt1.rds')
cesm_salt_dfs2 <- readRDS('F:/data/cesm_forecast_salt2.rds')
cesm_salt_dfs3 <- readRDS('F:/data/cesm_forecast_salt3.rds')
cesm_salt_hist <- filter(readRDS('F:/data/cesm_salt_dfs_trim.rds'),
                         projection == 'historical')

cesm_salt_sum1 <- summarise_years(cesm_salt_dfs1)
cesm_salt_sum2 <- summarise_years(cesm_salt_dfs2)
cesm_salt_sum3 <- summarise_years(cesm_salt_dfs3)
cesm_salt_hist_sum <- summarise_years_hist_salt(cesm_salt_hist)

cesm_dfs <- bind_rows(cesm_salt_sum1, 
                      cesm_salt_sum2,
                      cesm_salt_sum3,
                      cesm_salt_hist_sum)

ggplot(cesm_dfs, aes(x = year,
                     y = value,
                     group = projection,
                     fill = projection)) +
  geom_line(aes(group = projection,
                color = projection)) +
  geom_ribbon(aes(group = projection,
                  fill = projection,
                  ymin = min,
                  ymax = max),
              alpha = 0.2) +
  labs(title = "CESM Projections",
       y = "Salinity",
       x = "Year") +
  theme_bw() +
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        plot.title = element_text(size = 28),
        legend.text = element_text(size = 19),
        legend.title = element_text(size = 20),
        text = element_text(family = "serif"))

dev.copy(jpeg, 'F:/results/cesm_roms_salts.jpg', 
         height = 10, width = 20, units = 'in', res = 200)
dev.off()

# GFDL
gfdl_salt_dfs1 <- readRDS('F:/data/gfdl_forecast_salt1.rds')
gfdl_salt_dfs2 <- readRDS('F:/data/gfdl_forecast_salt2.rds')
gfdl_salt_dfs3 <- readRDS('F:/data/gfdl_forecast_salt3.rds')
gfdl_salt_hist <- filter(readRDS('F:/data/gfdl_salt_dfs_trim.rds'),
                         projection == 'historical')

gfdl_salt_sum1 <- summarise_years(gfdl_salt_dfs1)
gfdl_salt_sum2 <- summarise_years(gfdl_salt_dfs2)
gfdl_salt_sum3 <- summarise_years(gfdl_salt_dfs3)
gfdl_salt_hist_sum <- summarise_years_hist_salt(gfdl_salt_hist)

gfdl_dfs <- bind_rows(gfdl_salt_sum1, 
                      gfdl_salt_sum2,
                      gfdl_salt_sum3,
                      gfdl_salt_hist_sum)

ggplot(gfdl_dfs, aes(x = year,
                     y = value,
                     group = projection,
                     fill = projection)) +
  geom_line(aes(group = projection,
                color = projection)) +
  geom_ribbon(aes(group = projection,
                  fill = projection,
                  ymin = min,
                  ymax = max),
              alpha = 0.2) +
  labs(title = "GFDL Projections",
       y = "Salinity",
       x = "Year") +
  theme_bw() +
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        plot.title = element_text(size = 28),
        legend.text = element_text(size = 19),
        legend.title = element_text(size = 20),
        text = element_text(family = "serif"))

dev.copy(jpeg, 'F:/results/gfdl_roms_salts.jpg', 
         height = 10, width = 20, units = 'in', res = 200)
dev.off()

# MIROC
miroc_salt_dfs1 <- readRDS('F:/data/miroc_forecast_salt1.rds')
miroc_salt_dfs2 <- readRDS('F:/data/miroc_forecast_salt2.rds')
miroc_salt_dfs3 <- readRDS('F:/data/miroc_forecast_salt3.rds')
miroc_salt_hist <- filter(readRDS('F:/data/miroc_salt_dfs_trim.rds'),
                          projection == 'historical')

miroc_salt_sum1 <- summarise_years(miroc_salt_dfs1)
miroc_salt_sum2 <- summarise_years(miroc_salt_dfs2)
miroc_salt_sum3 <- summarise_years(miroc_salt_dfs3)
miroc_salt_hist_sum <- summarise_years_hist_salt(miroc_salt_hist)

miroc_dfs <- bind_rows(miroc_salt_sum1, 
                       miroc_salt_sum2,
                       miroc_salt_sum3,
                       miroc_salt_hist_sum)

ggplot(miroc_dfs, aes(x = year,
                      y = value,
                      group = projection,
                      fill = projection)) +
  geom_line(aes(group = projection,
                color = projection)) +
  geom_ribbon(aes(group = projection,
                  fill = projection,
                  ymin = min,
                  ymax = max),
              alpha = 0.2) +
  labs(title = "MIROC Projections",
       y = "Salinity",
       x = "Year") +
  theme_bw() +
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        plot.title = element_text(size = 28),
        legend.text = element_text(size = 19),
        legend.title = element_text(size = 20),
        text = element_text(family = "serif"))

dev.copy(jpeg, 'F:/results/miroc_roms_salts.jpg', 
         height = 10, width = 20, units = 'in', res = 200)
dev.off()