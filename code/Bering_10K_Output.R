# Title: Examine primary attributes of Bering 10K model.
# Authors: Lorenzo Ciannelli, Rebecca Howard
# 3/18/21

# See email from Kelly Kearney from 12/3/20
# See also B10K-K20_CORECFS_readme.pdf for explanation of netcdf files and roms 10k, and Bering 10K ROMS coordinates.pdf for working with coordinates

# Load libraries ----
library(maps)
library(mapdata)
library(marmap)
library(raster)
library(ncdf4)
library(fields)
library(pracma)
library(viridis)
library(itsadug)
library(date)
library(here)
library(tidyverse)
library(lubridate)

# Load data and necessary functions ----
source(here('code/functions', 'distance_function.R'))

# Depth data obtained from NGDC grid extract tool for ETOPO1
str_name <- (here('data', 'bering_bathy.tiff'))
bering_bathy <- as.bathy(raster(str_name) * -1)
bathy_lat <- as.numeric(colnames(bering_bathy))
bathy_lon <- as.numeric(rownames(bering_bathy))
bathy_ylim = range(bathy_lat)
bathy_xlim = range(bathy_lon)
bering_bathy[bering_bathy <= -1] <- NA

# Bering 10K model output: avg bottom temp 2015-2019
bering_model <-
  nc_open(here('data', 'B10K-K20_CORECFS_2015-2019_average_temp_bottom5m.nc'))
print(bering_model)

lon <- ncvar_get(bering_model, varid = 'lon_rho')
lon1 <- ifelse(lon >= 180, lon -360, lon)
lat <- ncvar_get(bering_model, varid = 'lat_rho')
time <- ncvar_get(bering_model, varid = 'ocean_time')
time1 <- as.Date(time / (60 * 60 * 24), origin = "1900-01-01 00:00:00")
fillvalue_t <- ncatt_get(bering_model, 'temp', "_FillValue")
temp <- ncvar_get(bering_model, varid = 'temp')
temp[temp == fillvalue_t$value] <- NA

dim(temp)
range(temp, na.rm = T)
dim(lat)
range(lat, na.rm = T)
dim(lon)
range(lon, na.rm = T)
dim(time1)
range(time1, na.rm = T)
nc_close(bering_model)

# Rotation by delta angle. Not needed
# delta <- 300
# newlon = lon * cos(delta * pi / 180) - lat * sin(delta * pi / 180)
# 
# newlat = lat * cos(delta * pi / 180) + lon * sin(delta * pi / 180)
# 
# range(newlon)
# range(newlat)

# Create figure of temperatures for 2015 - 2019 ----
# Get index for "month"
year <- unique(substr(time1, 1, 4))
month_year <- substr(time1, 1, 7)
month <- "06"

# Plot bottom temp in 'month'
windows(width = 9, height = 10)
par(mfrow = c(3, 2))
for (i in 1:length(year)) {
  starttime <-
    (1:length(month_year))[month_year == paste(year[i], month, sep = "-")][1]
  
  # Get temp
  tempplot <- temp[, , starttime]
  
  # Get date of image plot
  plotday1 <- as.character(time1[starttime])
  
  # show Temp
  image.plot(
    lon,
    lat,
    tempplot,
    ylab = "Latitude (north)",
    col = viridis(100),
    xlab = "Longitude (east)",
    main = plotday1,
    xlim = c(-180, -155) + 360,
    ylim = c(54, 64),
    zlim = c(-2, 17.5)
  )
  map("world2",
      fill = T,
      col = "grey",
      add = T)
}
dev.copy(
  jpeg,
  'results/Bottom_temp_2015_2019.jpg',
  height = 10,
  width = 9,
  res = 200,
  units = 'in'
)
dev.off()

# Add in pollock data ----
# Obtain temp values in correspondence of sampled stations from groundfish survey 
pk_adults_catch <- read_csv(here('data', 'EBS_POLL.csv'))
pk_adults_catch <- separate(pk_adults_catch, DAY_MONTH, c("day", "month"))
names(pk_adults_catch)[-2:-3] <- tolower(names(pk_adults_catch)[-2:-3])
pk_adults_catch <- pk_adults_catch %>% 
  mutate(date = make_date(year, month, day),
         cpue_noha = number_fish / (distance_fished * 1000 * net_width) * 10000)

head(pk_adults_catch)
pk_dat <- pk_adults_catch[pk_adults_catch$year > 2014, ]
pk_dat$roms_date <- NA
pk_dat$roms_temp <- NA

# Two versions: 1. Closest grid/time point; 2. loess interpolations

#1. Closest grid/time point
for (i in 1:nrow(pk_dat)) {
  idx_time <- order(abs(time1 - pk_dat$date[i]))[1]
  pk_dat$roms_date[i] <- time1[idx_time]
  idx_grid <- order(distance_function(
    pk_dat$start_latitude[i],
    pk_dat$start_longitude[i],
    c(lat),
    c(lon1)
  ))[1]
  pk_dat$roms_temp[i] <- c(temp[, , idx_time])[idx_grid]
}

#2. Loess interpolation by year, in the month of July
#Get index for "month"
year <- unique(pk_dat$year)
month_year <- substr(time1, 1, 7)
month <- "07"
pk_dat$loess_temp <- NA
pk_dat$loess_date <- NA
for (i in 1:length(year)) {
  idx.time <-
    (1:length(month_year))[month_year == paste(year[i], month, sep = "-")][1]
  # Get temp
  temp_plot <- temp[, , idx_time]
  # Get start and end date of data fields
  plot_day1 <- as.character(time.1[idx_time])
  # Make data frame for loess
  data_loess <-
    na.exclude(data.frame(
      z = c(temp_plot),
      x = c(lon1),
      y = c(lat)
    ))
  # Fit loess and check predictions
  temp.loess <- loess(z ~ x * y,
                      data = data.loess,
                      span = 0.002,
                      degree = 2)
  summary(lm(predict(temp.loess) ~ data.loess$z))
  # Make prediction data frame
  pred.loess <-
    data.frame(x = pk_dat$START_LONGITUDE[pk_dat$YEAR == year[i]],
               y = pk_dat$START_LATITUDE[pk_dat$YEAR == year[i]])
  # Add predictions to fish data
  pk_dat$loess.bt[pk_dat$YEAR == year[i]] <- predict(temp.loess,
                                                     newdata = pred.loess)
  pk_dat$loess.date[pk_dat$YEAR == year[i]] <- plot.day1
}
head(pk_dat)

#Check predictions
quartz(height = 8, width = 8)
par(mfrow = c(2, 2))
plot(
  pk_dat$roms.bt,
  pk_dat$GEAR_TEMPERATURE,
  xlab = "Closest ROMS Temp",
  ylab = "Gear temp",
  main = "Observation vs ROMS_close"
)
abline(0, 1, col = "red", lwd = 2)
plot(
  pk_dat$loess.bt,
  pk_dat$GEAR_TEMPERATURE,
  xlab = "LOESS Temp",
  ylab = "Gear temp",
  main = "Observation vs ROMS_loess"
)
abline(0, 1, col = "red", lwd = 2)
plot(
  pk_dat$loess.bt,
  pk_dat$roms.bt,
  xlab = "LOESS Temp",
  ylab = "Closest ROMS",
  main = "ROMS_loess vs ROMS_close"
)
abline(0, 1, col = "red", lwd = 2)
dev.copy(
  jpeg,
  'Model_Obs_Bottom_Temp.jpg',
  height = 8,
  width = 8,
  res = 200,
  units = 'in'
)
dev.off()

#Fit GAMs using observed and predicted bottom temp
#1. Observed BT
gam_obs <- gam(
  log(CPUE_NOHA) ~ factor(YEAR) +
    s(START_LONGITUDE, START_LATITUDE) +
    s(GEAR_TEMPERATURE, k = 4),
  data = pk_dat[pk_dat$CPUE_NOHA > 0, ]
)
summary(gam_obs)
#R-sq.(adj) =  0.383   Deviance explained = 39.3%
#GCV = 1.8872  Scale est. = 1.854     n = 1864

#2. ROMS close
gam_roms <- gam(
  log(CPUE_NOHA) ~ factor(YEAR) +
    s(START_LONGITUDE, START_LATITUDE) +
    s(roms.bt, k = 4),
  data = pk_dat[pk_dat$CPUE_NOHA > 0, ]
)
summary(gam_roms)
#R-sq.(adj) =  0.344   Deviance explained = 35.6%
#GCV = 2.0006  Scale est. = 1.9656    n = 1863

#3. ROMS loess
gam_loess <- gam(
  log(CPUE_NOHA) ~ factor(YEAR) +
    s(START_LONGITUDE, START_LATITUDE) +
    s(loess.bt, k = 4),
  data = pk_dat[pk_dat$CPUE_NOHA > 0, ]
)
summary(gam_loess)
#R-sq.(adj) =  0.345   Deviance explained = 35.6%
#GCV = 2.0018  Scale est. = 1.968     n = 1864

#Plot all three temp effects
quartz(height = 8, width = 8)
par(mfrow = c(2, 2))
plot(
  gam_obs,
  select = 2,
  scale = 0,
  shade = T,
  res = T,
  main = "Observed temperature"
)
plot(
  gam_roms,
  select = 2,
  scale = 0,
  shade = T,
  res = T,
  main = "ROMS_close temperature"
)
plot(
  gam_loess,
  select = 2,
  scale = 0,
  shade = T,
  res = T,
  main = "ROMS_loess temperature"
)
dev.copy(
  jpeg,
  'GAM_Bottom_Temp.jpg',
  height = 8,
  width = 8,
  res = 200,
  units = 'in'
)
dev.off()
