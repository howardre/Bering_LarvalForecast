# Exploratory analysis of hindcast data
# 5/25/2021

### Load libraries & functions ----
library(marmap)
library(maps)
library(mapdata)
library(fields)
library(here)
library(mgcv)
library(tidyr)
library(dplyr)
library(ggplot2)
library(raster)
library(rgdal)
library(viridis)
source(here('code/functions', 'vis_gam_COLORS.R'))
source(here('code/functions', 'distance_function.R'))

# Function to make maps
map_phenology <- function(data, grid, grid2){
  nlat = 80
  nlon = 120
  latd = seq(min(grid$lat), max(grid$lat), length.out = nlat)
  lond = seq(min(grid$lon), max(grid$lon), length.out = nlon)
  my_color = colorRampPalette(c("#1C0D51", "#4C408E", "#7E77B0",
                                "#AFABCB", "#DAD9E5", "#F9F9F9",
                                "#FFDAB7", "#FFB377","#E18811",
                                "#AC6000", "#743700"))
  color_levels = 100
  max_absolute_value = max(abs(c(min(grid$pred, na.rm = T), max(grid$pred, na.rm = T)))) 
  color_sequence = seq(-max_absolute_value, max_absolute_value, 
                       length.out = color_levels + 1)
  n_in_class = hist(grid$pred, breaks = color_sequence, plot = F)$counts > 0
  col_to_include = min(which(n_in_class == T)):max(which(n_in_class == T))
  breaks_to_include = min(which(n_in_class == T)):(max(which(n_in_class == T)) + 1)
  image(lond,
        latd,
        t(matrix(grid$pred,
                 nrow = length(latd),
                 ncol = length(lond),
                 byrow = T)),
        xlim = c(-176.5, -156.5),
        ylim = c(52, 62),
        axes = FALSE,
        xlab = "",
        ylab = "")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "mintcream")
  par(new = TRUE)
  image(lond,
        latd,
        t(matrix(grid$pred,
                 nrow = length(latd),
                 ncol = length(lond),
                 byrow = T)),
        col = my_color(n = color_levels)[col_to_include], 
        breaks = color_sequence[breaks_to_include],
        ylab = "Latitude",
        xlab = "Longitude",
        xlim = c(-176.5, -156.5),
        ylim = c(52, 62),
        main = "Distribution",
        cex.main = 1.2,
        cex.lab = 1.1,
        cex.axis = 1.1)
  symbols(data$lon[data$larvalcatchper10m2 > 0],
          data$lat[data$larvalcatchper10m2 > 0],
          circles = log(data$larvalcatchper10m2 + 1)[data$larvalcatchper10m2 > 0],
          inches = 0.1,
          bg = alpha('grey', 0.1),
          fg = alpha('black', 0.05),
          add = T)
  points(data$lon[data$larvalcatchper10m2 == 0], 
         data$lat[data$larvalcatchper10m2 == 0], 
         pch = '')
  maps::map("worldHires",
            fill = T,
            col = "wheat4",
            add = T)
  image.plot(legend.only = T,
             col = my_color(n = color_levels)[col_to_include],
             legend.shrink = 0.2,
             smallplot = c(.765, .79, .25, .37),
             legend.cex = 0.4,
             axis.args = list(cex.axis = 0.6),
             legend.width = 0.5,
             legend.mar = 6,
             zlim = c(min(grid$pred, na.rm = T), max(grid$pred, na.rm = T)),
             legend.args = list("Predicted \n Change",
                                side = 2, cex = 0.7))
  plot(grid2$doy,
       grid2$pred,
       main = 'Phenology',
       type = 'l',
       ylab = 'Egg density ln(n/10m2)',
       xlab = 'Day of the year',
       cex.lab = 1.1,
       cex.axis = 1.1,
       cex.main = 1.2,
       xlim = c(min(grid2$doy, na.rm = T), max(grid2$doy, na.rm = T)),
       ylim = range(c(grid2$pred_up, grid2$pred_lw)),
       col = 'blue',
       lwd = 2)
  polygon(c(grid2$doy, rev(grid2$doy)),
          c(grid2$pred_lw, rev(grid2$pred_up)),
          col = alpha('gray', 0.6),
          lty = 0)
}

### Load fish data ----
yfs_egg <- readRDS(here('data', 'yfs_egg.rds'))
yfs_larvae <- readRDS(here('data', 'yfs_larvae.rds'))
akp_egg <- readRDS(here('data', 'akp_egg.rds'))
akp_larvae <- readRDS(here('data', 'akp_larvae.rds'))
fhs_egg <- readRDS(here('data', 'fhs_egg.rds'))
fhs_larvae <- readRDS(here('data', 'fhs_larvae.rds'))
pk_egg <- readRDS(here('data', 'pk_egg.rds'))
pk_larvae <- readRDS(here('data', 'pk_larvae.rds'))

### Load bathymetry data ----
str_name <- (here('data', 'bering_bathy.tiff'))
bering_bathy <- as.bathy(raster(str_name) * -1)
bathy_lat <- as.numeric(colnames(bering_bathy))
bathy_lon <- as.numeric(rownames(bering_bathy))
bathy_ylim = range(bathy_lat)
bathy_xlim = range(bathy_lon)
bering_bathy[bering_bathy <= -1] <- NA

## NCEP SST data ----
# https://www.esrl.noaa.gov/psd/data/timeseries/ 
# Latitude Range used: 56.2 to  56.2, Longitude Range used: 195.0 to 196.9
ncep_temp <- read.csv(here('data', 'NCEP_Temp.csv'))

### Data exploration ----
# Yellowfin
ggplot(yfs_egg) +
  geom_point(aes(larvalcatchper10m2, year)) # seems to be an outlier catch (~340) <- remove
ggplot(yfs_larvae) +
  geom_point(aes(larvalcatchper10m2, year)) # outlier above 30000 <- remove

# AKP
ggplot(akp_egg) +
  geom_point(aes(larvalcatchper10m2, year)) # possibly several outliers
# two really high above 7000
# 6 above 4000 <- likely remove these
ggplot(akp_larvae) +
  geom_point(aes(larvalcatchper10m2, year)) # remove above 500

# FHS
ggplot(fhs_egg) +
  geom_point(aes(larvalcatchper10m2, year)) # outliers above ~1250 <- remove these
ggplot(fhs_larvae) +
  geom_point(aes(larvalcatchper10m2, year)) # outlier above 9000 <- remove

# pollock
ggplot(pk_egg) +
  geom_point(aes(larvalcatchper10m2, year))
# 3 extreme above 3.5 * 10^5 <- likely remove these
# 7 above 1.25 * 10^5
ggplot(pk_larvae) +
  geom_point(aes(larvalcatchper10m2, year)) # remove above 1*10^5

# year
ggplot(yfs_egg) +
  geom_bar(aes(year))  
ggplot(pk_egg) +
  geom_bar(aes(year))
ggplot(akp_egg) +
  geom_bar(aes(year))
ggplot(fhs_egg) +
  geom_bar(aes(year))
# higher number of catches in recent years, post 2001

# day of year
ggplot(yfs_egg) +
  geom_bar(aes(doy)) 
ggplot(pk_egg) +
  geom_bar(aes(doy))
ggplot(fhs_egg) +
  geom_bar(aes(doy))
ggplot(akp_egg) +
  geom_bar(aes(doy))
# all have few catches prior to day 99, may not be too much of an issue

# temperature
ggplot(pk_larvae) +
  geom_point(aes(roms_temperature, year))
ggplot(yfs_larvae) +
  geom_point(aes(roms_temperature, year))
ggplot(fhs_larvae) +
  geom_point(aes(roms_temperature, year))
ggplot(akp_larvae) +
  geom_point(aes(roms_temperature, year))
ggplot(pk_egg) +
  geom_point(aes(roms_temperature, year))
ggplot(yfs_egg) +
  geom_point(aes(roms_temperature, year))
ggplot(fhs_egg) +
  geom_point(aes(roms_temperature, year))
ggplot(akp_egg) +
  geom_point(aes(roms_temperature, year))
# nothing of note

# salinity
ggplot(pk_larvae) +
  geom_point(aes(roms_salinity, year)) # outlier below 29?
ggplot(yfs_larvae) +
  geom_point(aes(roms_salinity, year)) # weird pattern
# clustered at high salinity, possible outliers below 25
ggplot(akp_larvae) +
  geom_point(aes(roms_salinity, year)) # outlier below 29
ggplot(fhs_larvae) +
  geom_point(aes(roms_salinity, year)) 
ggplot(pk_egg) +
  geom_point(aes(roms_salinity, year)) # outliers below 29?
ggplot(yfs_egg) +
  geom_point(aes(roms_salinity, year))
ggplot(fhs_egg) +
  geom_point(aes(roms_salinity, year)) # outlier below 29
ggplot(akp_egg) +
  geom_point(aes(roms_salinity, year))

# Remove outliers?
# yfs_egg <- filter(yfs_egg, larvalcatchper10m2 < 330)
# yfs_larvae <- filter(yfs_larvae, larvalcatchper10m2 < 30000,
#                      roms_salinity > 25)
# akp_egg <- filter(akp_egg, larvalcatchper10m2 < 4000)
# akp_larvae <- filter(akp_larvae, larvalcatchper10m2 < 500,
#                      roms_salinity > 29)
# fhs_egg <- filter(fhs_egg, larvalcatchper10m2 < 1400,
#                   roms_salinity > 29)
# fhs_larvae <- filter(fhs_larvae, larvalcatchper10m2 < 9000)
# pk_egg <- filter(pk_egg, larvalcatchper10m2 < 350000,
#                  roms_salinity > 29)
# pk_larvae <- filter(pk_larvae, larvalcatchper10m2 < 100000,
#                     roms_salinity > 29)

# only remove catches if there are few for a year
# use something like 50-60 for year
# keep early months?

# map of data
# Create bathymetry dataset
BS_bathy <- getNOAA.bathy(lon1 = -180, 
                          lon2 = -156, 
                          lat1 = 66, 
                          lat2 = 52, 
                          resolution = 1)
blues <- c("lightsteelblue4", 
           "lightsteelblue3", 
           "lightsteelblue2", 
           "lightsteelblue1")
greys <- c(grey(0.6), 
           grey(0.93), 
           grey(0.99))

windows()
plot.bathy(
  BS_bathy,
  image = T,
  axes = T,
  lwd = 0.03,
  land = T,
  n = 0,
  bpal = list(c(0, max(BS_bathy), greys), c(min(BS_bathy), 0, blues)),
  ylim = c(52, 66),
  xlim = c(-180,-156),
  ylab = "Latitude °N",
  xlab = "Longitude °W",
  cex.lab = 1,
  cex.main = 1,
  cex.axis = 1)
points(pk_egg$lon[pk_egg$larvalcatchper10m2 == 0],
       pk_egg$lat[pk_egg$larvalcatchper10m2 == 0],
       pch = 4,
       col = 'darkgray',
       cex = 1)
points(pk_egg$lon[pk_egg$larvalcatchper10m2 > 0],
       pk_egg$lat[pk_egg$larvalcatchper10m2 > 0],
       pch = 16,
       col = 'black',
       cex = 0.7)
# doesn't seem to have anything weird going on
# was fixed with the distance function during cleaning, don't need to do this for all species


### Yellowfin sole ----
#### Eggs ----
hist(yfs_egg$larvalcatchper10m2)

# Negative Binomial
yfs_egg_basenb <- gam(count ~ factor(year) + 
                       s(lon, lat) + 
                       s(doy) +
                       s(roms_salinity) +
                       s(roms_temperature),
                     data = yfs_egg,
                     family = nb(),
                     offset = log(volume_filtered))
summary(yfs_egg_basenb)
# R2: 0.142
# Deviance: 84.3%

yfs_egg_basenb$family$getTheta(TRUE) # 0.1022598

windows()
par(mfrow = c(2, 2))
plot(yfs_egg_basenb, select = 2, main = "DOY")
plot(yfs_egg_basenb, select = 3, main = "Salinity")
plot(yfs_egg_basenb, select = 4, main = "Temperature")
dev.copy(jpeg, here('results/yellowfin_hindcast', 'yellowfin_egg_basenb.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(yfs_egg_basenb)

# quasiPoisson
yfs_egg_basep <- gam(count ~ factor(year) + 
                      s(lon, lat) + 
                      s(doy) +
                      s(roms_salinity) +
                      s(roms_temperature),
                    data = yfs_egg,
                    family = quasipoisson(link = "log"),
                    offset = log(volume_filtered))
summary(yfs_egg_basep)
# R2: 0.996
# Deviance: 98.5%

windows()
par(mfrow = c(2, 2))
plot(yfs_egg_basep, select = 2, main = "DOY")
plot(yfs_egg_basep, select = 3, main = "Salinity")
plot(yfs_egg_basep, select = 4, main = "Temperature")
dev.copy(jpeg, here('results/yellowfin_hindcast', 'yellowfin_egg_basep.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(yfs_egg_basep)

# Tweedie
yfs_egg_baset <- gam(larvalcatchper10m2 ~ factor(year) + 
                      s(lon, lat) + 
                      s(doy) +
                      s(roms_salinity) +
                      s(roms_temperature),
                    data = yfs_egg,
                    family = tw(link = 'log'),
                    method = 'REML',
                    offset = log(volume_filtered))
summary(yfs_egg_baset)
# R2: 0.474
# Deviance: 81.5%

windows()
par(mfrow = c(2, 2))
plot(yfs_egg_baset, select = 2, main = "DOY")
plot(yfs_egg_baset, select = 3, main = "Salinity")
plot(yfs_egg_baset, select = 4, main = "Temperature")
dev.copy(jpeg, here('results/yellowfin_hindcast', 'yellowfin_egg_baset.jpg'), 
         height = 5, width = 5, units = 'in', res = 200)
dev.off()

par(mfrow = c(2, 2))
gam.check(yfs_egg_baset)

# Zero-inflated Poisson (1 stage and 2 stage)
yfs_egg_zip <- gam(count ~ s(year) +  # cannot use year as a factor
                    s(lon, lat) + 
                    s(doy) +
                    s(roms_salinity) +
                    s(roms_temperature),
                  data = yfs_egg,
                  family = ziP(),
                  offset = log(volume_filtered))
summary(yfs_egg_zip)
# Deviance: 89.9%

windows()
par(mfrow = c(2, 2))
plot(yfs_egg_zip, select = 3, main = "DOY")
plot(yfs_egg_zip, select = 4, main = "Salinity")
plot(yfs_egg_zip, select = 5, main = "Temperature")
dev.copy(jpeg, here('results/yellowfin_hindcast', 'yellowfin_egg_zip.jpg'), 
         height = 5, width = 5, units = 'in', res = 200)
dev.off()

par(mfrow = c(2, 2))
gam.check(yfs_egg_zip)

yfs_egg_ziplss <- gam(list(count ~ s(year) +
                            s(lon, lat) +
                            s(doy) +
                            s(roms_salinity) +
                            s(roms_temperature),
                          ~ s(year) +
                            s(lon, lat) +
                            s(doy)),
                     offset = log(volume_filtered),
                     data = yfs_egg,
                     family = ziplss()) #ziplss
summary(yfs_egg_ziplss)
# Deviance: 69.3%

windows()
par(mfrow = c(2, 2))
plot(yfs_egg_ziplss, select = 3, main = "DOY")
plot(yfs_egg_ziplss, select = 4, main = "Salinity")
plot(yfs_egg_ziplss, select = 5, main = "Temperature")
dev.copy(jpeg, here('results/yellowfin_hindcast', 'yellowfin_egg_ziplss.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(yfs_egg_ziplss)

# Two-part binomial and gaussian
# binomial
yfs_egg$presence <- 1 * (yfs_egg$count > 0)
yfs_egg_gam1 <- gam(presence ~ factor(year) +
                     s(doy) +
                     s(lon, lat),
                   data = yfs_egg,
                   family = "binomial",
                   offset = log(volume_filtered))
summary(yfs_egg_gam1)
# R2: 0.382
# Deviance: 61.1%

windows()
plot(yfs_egg_gam1, select = 1, main = "DOY")
dev.copy(jpeg, here('results/yellowfin_hindcast', 'yellowfin_egg_gam1.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(yfs_egg_gam1)

# gaussian
yfs_egg_gam2 <- gam(log(larvalcatchper10m2 + 1) ~ factor(year) +
                      s(doy) +
                      s(lon, lat, k = 4),
                    offset = log(volume_filtered),
                    data = yfs_egg[yfs_egg$larvalcatchper10m2 > 0, ])
summary(yfs_egg_gam2)
# R2: 0.0605
# Deviance: 64.2% 

windows()
par(mfrow = c(2, 2))
plot(yfs_egg_gam2, select = 1, main = "DOY")
plot(yfs_egg_gam2, select = 3, main = "Salinity")
plot(yfs_egg_gam2, select = 4, main = "Temperature")
dev.copy(jpeg, here('results/yellowfin_hindcast', 'yellowfin_egg_gam2.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(yfs_egg_gam2)


# Plot best model
# Plot results of average geography and phenology along with decrease of MSE
# Prediction grid
nlat = 80
nlon = 120
latd = seq(min(yfs_egg$lat), max(yfs_egg$lat), length.out = nlat)
lond = seq(min(yfs_egg$lon), max(yfs_egg$lon), length.out = nlon)
grid_extent_eggyfs <- expand.grid(lond, latd)
names(grid_extent_eggyfs) <- c('lon', 'lat')

# Calculate distance of each grid point to closest 'positive observation'
grid_extent_eggyfs$dist <- NA
for (k in 1:nrow(grid_extent_eggyfs)) {
  dist <- distance_function(grid_extent_eggyfs$lat[k],
                            grid_extent_eggyfs$lon[k],
                            yfs_egg$lat,
                            yfs_egg$lon)
  grid_extent_eggyfs$dist[k] <- min(dist)
}

# Assign a within sample year and doy to the grid data
yfs_egg_model1 <- yfs_egg_gam1
yfs_egg_model2 <- yfs_egg_gam2
grid_extent_eggyfs$year <- 2014
grid_extent_eggyfs$doy <- median(yfs_egg$doy)
grid_extent_eggyfs$pred1 <- predict(yfs_egg_model1, newdata = grid_extent_eggyfs)
grid_extent_eggyfs$pred2 <- predict(yfs_egg_model2, newdata = grid_extent_eggyfs)
grid_extent_eggyfs$pred <- grid_extent_eggyfs$pred1 * grid_extent_eggyfs$pred2
grid_extent_eggyfs$pred[grid_extent_eggyfs$dist > 30000] <- NA

# Phenology
grid_extent_eggyfs2 <- data.frame('lon' = rep(-170, 100),
                                 'lat' = rep(57, 100),
                                 'doy' = seq(min(yfs_egg$doy), 
                                             max(yfs_egg$doy), 
                                             length = 100),
                                 'year' = rep(2014, 100))
grid_extent_eggyfs2$pred1 <- predict(yfs_egg_model1, newdata = grid_extent_eggyfs2)
grid_extent_eggyfs2$pred2 <- predict(yfs_egg_model2, newdata = grid_extent_eggyfs2)
grid_extent_eggyfs2$pred <- grid_extent_eggyfs2$pred1 * grid_extent_eggyfs2$pred2
grid_extent_eggyfs2$se1 <- predict(yfs_egg_model1, newdata = grid_extent_eggyfs2, se = T)[[2]]
grid_extent_eggyfs2$se2 <- predict(yfs_egg_model2, newdata = grid_extent_eggyfs2, se = T)[[2]]
grid_extent_eggyfs2$se <- grid_extent_eggyfs2$se1 * grid_extent_eggyfs2$se2
grid_extent_eggyfs2$pred_up <- grid_extent_eggyfs2$pred + 1.96 * grid_extent_eggyfs2$se
grid_extent_eggyfs2$pred_lw <- grid_extent_eggyfs2$pred - 1.96 * grid_extent_eggyfs2$se

# Plot
windows(width = 8, height = 3.5)
par(mfrow = c(1, 2), 
    mai = c(0.8, 0.9, 0.5, 0.5))
map_phenology(yfs_egg, grid_extent_eggyfs, grid_extent_eggyfs2)
dev.copy(jpeg, here('results/yellowfin_hindcast', 'yellowfin_base_egg.jpg'), 
         height = 3.5, width = 8, res = 200, units = 'in')
dev.off()

#### Larvae ----
# Negative Binomial
yfs_larvae_basenb <- gam(count ~ factor(year) + 
                          s(lon, lat) + 
                          s(doy) +
                          s(roms_salinity) +
                          s(roms_temperature),
                        data = yfs_larvae,
                        family = nb(),
                        offset = log(volume_filtered))
summary(yfs_larvae_basenb)
# R2: -1.54e+03
# Deviance: 79.9%

yfs_larvae_basenb$family$getTheta(TRUE) # 0.3223544

windows()
par(mfrow = c(2, 2))
plot(yfs_larvae_basenb, select = 2, main = "DOY")
plot(yfs_larvae_basenb, select = 3, main = "Salinity")
plot(yfs_larvae_basenb, select = 4, main = "Temperature")
dev.copy(jpeg, here('results/yellowfin_hindcast', 'yellowfin_larvae_basenb.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(yfs_larvae_basenb)

# quasiPoisson
yfs_larvae_basep <- gam(count ~ factor(year) + 
                         s(lon, lat) + 
                         s(doy) +
                         s(roms_salinity) +
                         s(roms_temperature),
                       data = yfs_larvae,
                       family = quasipoisson(link = "log"),
                       offset = log(volume_filtered))
summary(yfs_larvae_basep)
# R2: -0.00972
# Deviance: 69.9%

windows()
par(mfrow = c(2, 2))
plot(yfs_larvae_basep, select = 2, main = "DOY")
plot(yfs_larvae_basep, select = 3, main = "Salinity")
plot(yfs_larvae_basep, select = 4, main = "Temperature")
dev.copy(jpeg, here('results/yellowfin_hindcast', 'yellowfin_larvae_basep.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(yfs_larvae_basep)

# Tweedie
yfs_larvae_baset <- gam(larvalcatchper10m2 ~ factor(year) + 
                         s(lon, lat) + 
                         s(doy) +
                         s(roms_salinity) +
                         s(roms_temperature),
                       data = yfs_larvae,
                       family = tw(link = 'log'),
                       method = 'REML')
summary(yfs_larvae_baset)
# R2: 0.0891
# Deviance: 74.1%

windows()
par(mfrow = c(2, 2))
plot(yfs_larvae_baset, select = 2, main = "DOY")
plot(yfs_larvae_baset, select = 3, main = "Salinity")
plot(yfs_larvae_baset, select = 4, main = "Temperature")
dev.copy(jpeg, here('results/yellowfin_hindcast', 'yellowfin_larvae_baset.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(yfs_larvae_baset)

# Zero-inflated Poisson (1 stage and 2 stage)
yfs_larvae_zip <- gam(count ~ s(year) +  # cannot use year as a factor
                       s(lon, lat) + 
                       s(doy) +
                       s(roms_salinity) +
                       s(roms_temperature),
                     data = yfs_larvae,
                     family = ziP(),
                     offset = log(volume_filtered))
summary(yfs_larvae_zip)
# Deviance: 26.6% ?

windows()
par(mfrow = c(2, 2))
plot(yfs_larvae_zip, select = 3, main = "DOY")
plot(yfs_larvae_zip, select = 4, main = "Salinity")
plot(yfs_larvae_zip, select = 5, main = "Temperature")
dev.copy(jpeg, here('results/yellowfin_hindcast', 'yellowfin_larvae_zip.jpg'), 
         height = 5, width = 5, units = 'in', res = 200)
dev.off()

par(mfrow = c(2, 2))
gam.check(yfs_larvae_zip)

yfs_larvae_ziplss <- gam(list(count ~ s(year) +
                               s(lon, lat) +
                               s(doy) +
                               s(roms_salinity) +
                               s(roms_temperature),
                             ~ s(year) +
                               s(lon, lat) +
                               s(doy)),
                        offset = log(volume_filtered),
                        data = yfs_larvae,
                        family = ziplss()) #ziplss
summary(yfs_larvae_ziplss)
# Deviance: 39.8%

windows()
par(mfrow = c(2, 2))
plot(yfs_larvae_ziplss, select = 3, main = "DOY")
plot(yfs_larvae_ziplss, select = 4, main = "Salinity")
plot(yfs_larvae_ziplss, select = 5, main = "Temperature")
dev.copy(jpeg, here('results/yellowfin_hindcast', 'yellowfin_larvae_ziplss.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(yfs_larvae_ziplss)

# Two-part binomial and gaussian
# binomial
yfs_larvae$presence <- 1 * (yfs_larvae$count > 0)
yfs_larvae_gam1 <- gam(presence ~ factor(year) +
                        s(doy) +
                        s(lon, lat),
                      data = yfs_larvae,
                      family = "binomial")
summary(yfs_larvae_gam1)
# R2: 0.557
# Deviance: 55.3%

windows()
plot(yfs_larvae_gam1, select = 1, main = "DOY")
dev.copy(jpeg, here('results/yellowfin_hindcast', 'yellowfin_larvae_gam1.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(yfs_larvae_gam1)

# gaussian
yfs_larvae_gam2 <- gam(log(larvalcatchper10m2 + 1) ~ factor(year) +
                         s(doy) +
                         s(lon, lat),
                      data = yfs_larvae[yfs_larvae$larvalcatchper10m2 > 0, ])
summary(yfs_larvae_gam2)
# R2: 0.339
# Deviance: 38.8% 

windows()
par(mfrow = c(2, 2))
plot(yfs_larvae_gam2, select = 1, main = "DOY")
plot(yfs_larvae_gam2, select = 3, main = "Salinity")
plot(yfs_larvae_gam2, select = 4, main = "Temperature")
dev.copy(jpeg, here('results/yellowfin_hindcast', 'yellowfin_larvae_gam2.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(yfs_larvae_gam2)

# Plot best model
# Plot results of average geography and phenology along with decrease of MSE
# Prediction grid
nlat = 80
nlon = 120
latd = seq(min(yfs_larvae$lat), max(yfs_larvae$lat), length.out = nlat)
lond = seq(min(yfs_larvae$lon), max(yfs_larvae$lon), length.out = nlon)
grid_extent_larvaeyfs <- expand.grid(lond, latd)
names(grid_extent_larvaeyfs) <- c('lon', 'lat')

# Calculate distance of each grid point to closest 'positive observation'
grid_extent_larvaeyfs$dist <- NA
for (k in 1:nrow(grid_extent_larvaeyfs)) {
  dist <- distance_function(grid_extent_larvaeyfs$lat[k],
                            grid_extent_larvaeyfs$lon[k],
                            yfs_larvae$lat,
                            yfs_larvae$lon)
  grid_extent_larvaeyfs$dist[k] <- min(dist)
}

# Assign a within sample year and doy to the grid data
yfs_larvae_model1 <- yfs_larvae_gam1
yfs_larvae_model2 <- yfs_larvae_gam2
grid_extent_larvaeyfs$year <- 2014
grid_extent_larvaeyfs$doy <- median(yfs_larvae$doy)
grid_extent_larvaeyfs$pred1 <- predict(yfs_larvae_model1, newdata = grid_extent_larvaeyfs)
grid_extent_larvaeyfs$pred2 <- predict(yfs_larvae_model2, newdata = grid_extent_larvaeyfs)
grid_extent_larvaeyfs$pred <- grid_extent_larvaeyfs$pred1 * grid_extent_larvaeyfs$pred2
grid_extent_larvaeyfs$pred[grid_extent_larvaeyfs$dist > 30000] <- NA

# Plot phenology
grid_extent_larvaeyfs2 <- data.frame('lon' = rep(-170, 100),
                                    'lat' = rep(57, 100),
                                    'doy' = seq(min(yfs_larvae$doy), max(yfs_larvae$doy), length = 100),
                                    'year' = rep(2014, 100))
grid_extent_larvaeyfs2$pred1 <- predict(yfs_larvae_model1, newdata = grid_extent_larvaeyfs2)
grid_extent_larvaeyfs2$pred2 <- predict(yfs_larvae_model2, newdata = grid_extent_larvaeyfs2)
grid_extent_larvaeyfs2$pred <- grid_extent_larvaeyfs2$pred1 * grid_extent_larvaeyfs2$pred2
grid_extent_larvaeyfs2$se1 <- predict(yfs_larvae_model1, newdata = grid_extent_larvaeyfs2, se = T)[[2]]
grid_extent_larvaeyfs2$se2 <- predict(yfs_larvae_model2, newdata = grid_extent_larvaeyfs2, se = T)[[2]]
grid_extent_larvaeyfs2$se <- grid_extent_larvaeyfs2$se1 * grid_extent_larvaeyfs2$se2 # need to change this, can't just be multiplied
grid_extent_larvaeyfs2$pred_up <- grid_extent_larvaeyfs2$pred + 1.96 * grid_extent_larvaeyfs2$se
grid_extent_larvaeyfs2$pred_lw <- grid_extent_larvaeyfs2$pred - 1.96 * grid_extent_larvaeyfs2$se

# Plot
windows(width = 8, height = 3.5)
par(mfrow = c(1, 2), 
    mai = c(0.8, 0.9, 0.5, 0.5))
map_phenology(yfs_larvae, grid_extent_larvaeyfs, grid_extent_larvaeyfs2)
dev.copy(jpeg, here('results/yellowfin_hindcast', 'yellowfin_base_larvae.jpg'), height = 3.5, width = 8, res = 200, units = 'in')
dev.off()

### Walleye Pollock ----
#### Eggs ----
hist(pk_egg$larvalcatchper10m2)

# Negative Binomial
pk_egg_basenb <- gam(count ~ factor(year) + 
                       s(lon, lat) + 
                       s(doy) +
                       s(roms_salinity) +
                       s(roms_temperature),
                     data = pk_egg,
                     family = nb(),
                     offset = log(volume_filtered))
summary(pk_egg_basenb)
# R2: -1.18
# Deviance: 68.7%

pk_egg_basenb$family$getTheta(TRUE) # 0.2818324

windows()
par(mfrow = c(2, 2))
plot(pk_egg_basenb, select = 2, main = "DOY")
plot(pk_egg_basenb, select = 3, main = "Salinity")
plot(pk_egg_basenb, select = 4, main = "Temperature")
dev.copy(jpeg, here('results/pollock_hindcast', 'pollock_egg_basenb.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(pk_egg_basenb)

# quasiPoisson
pk_egg_basep <- gam(count ~ factor(year) + 
                      s(lon, lat) + 
                      s(doy) +
                      s(roms_salinity) +
                      s(roms_temperature),
                    data = pk_egg,
                    family = quasipoisson(link = "log"),
                    offset = log(volume_filtered))
summary(pk_egg_basep)
# R2: 0.315
# Deviance: 68.5%

windows()
par(mfrow = c(2, 2))
plot(pk_egg_basep, select = 2, main = "DOY")
plot(pk_egg_basep, select = 3, main = "Salinity")
plot(pk_egg_basep, select = 4, main = "Temperature")
dev.copy(jpeg, here('results/pollock_hindcast', 'pollock_egg_basep.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(pk_egg_basep)

# Tweedie
pk_egg_baset <- gam(larvalcatchper10m2 ~ factor(year) + 
                      s(lon, lat) + 
                      s(doy) +
                      s(roms_salinity) +
                      s(roms_temperature),
                    data = pk_egg,
                    family = tw(link = 'log'),
                    method = 'REML',
                    offset = log(volume_filtered))
summary(pk_egg_baset)
# R2: -0.544
# Deviance: 66.5%

windows()
par(mfrow = c(2, 2))
plot(pk_egg_baset, select = 2, main = "DOY")
plot(pk_egg_baset, select = 3, main = "Salinity")
plot(pk_egg_baset, select = 4, main = "Temperature")
dev.copy(jpeg, here('results/pollock_hindcast', 'pollock_egg_baset.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(pk_egg_baset)

plot(pk_egg_baset, select = 2, main = "Tweedie")

# Zero-inflated Poisson (1 stage and 2 stage)
pk_egg_zip <- gam(count ~ s(year) +  # cannot use year as a factor
                    s(lon, lat) + 
                    s(doy) +
                    s(roms_salinity) +
                    s(roms_temperature),
                  data = pk_egg,
                  family = ziP(),
                  offset = log(volume_filtered))
summary(pk_egg_zip)
# Deviance: 100%

windows()
par(mfrow = c(2, 2))
plot(pk_egg_zip, select = 3, main = "DOY")
plot(pk_egg_zip, select = 4, main = "Salinity")
plot(pk_egg_zip, select = 5, main = "Temperature")
dev.copy(jpeg, here('results/pollock_hindcast', 'pollock_egg_zip.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(pk_egg_zip)

plot(pk_egg_zip, select = 3, main = "Zero-inflated Poisson (1-stage)")

pk_egg_ziplss <- gam(list(count ~ s(year) +
                            s(lon, lat) +
                            s(doy) +
                            s(roms_salinity) +
                            s(roms_temperature),
                          ~ s(year) +
                            s(lon, lat) +
                            s(doy)),
                         offset = log(volume_filtered),
                         data = pk_egg,
                         family = ziplss()) #ziplss
summary(pk_egg_ziplss)
# Deviance: 64.6%

windows()
par(mfrow = c(2, 2))
plot(pk_egg_ziplss, select = 2, main = "DOY")
plot(pk_egg_ziplss, select = 3, main = "Salinity")
plot(pk_egg_ziplss, select = 4, main = "Temperature")
dev.copy(jpeg, here('results/pollock_hindcast', 'pollock_egg_ziplss.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(pk_egg_ziplss)

# Two-part binomial and gaussian
# binomial
pk_egg$presence <- 1 * (pk_egg$count > 0)
pk_egg_gam1 <- gam(presence ~ factor(year) +
                     s(doy) +
                     s(lon, lat),
                   data = pk_egg,
                   family = "binomial")
summary(pk_egg_gam1)
# R2: 0.487
# Deviance: 42.1

par(mfrow = c(2, 2))
gam.check(pk_egg_gam1)

# gaussian
pk_egg_gam2 <- gam(log(larvalcatchper10m2 + 1) ~ factor(year) +
                     s(doy) +
                     s(lon, lat),
                   data = pk_egg[pk_egg$larvalcatchper10m2 > 0, ])
summary(pk_egg_gam2)
# R2: 0.411
# Deviance: 42.5% 

par(mfrow = c(2, 2))
gam.check(pk_egg_gam2)


# Plot best model
# Plot results of average geography and phenology along with decrease of MSE
# Prediction grid
nlat = 80
nlon = 120
latd = seq(min(pk_egg$lat), max(pk_egg$lat), length.out = nlat)
lond = seq(min(pk_egg$lon), max(pk_egg$lon), length.out = nlon)
grid_extent_eggpk <- expand.grid(lond, latd)
names(grid_extent_eggpk) <- c('lon', 'lat')

# Calculate distance of each grid point to closest 'positive observation'
grid_extent_eggpk$dist <- NA
for (k in 1:nrow(grid_extent_eggpk)) {
  dist <- distance_function(grid_extent_eggpk$lat[k],
                            grid_extent_eggpk$lon[k],
                            pk_egg$lat,
                            pk_egg$lon)
  grid_extent_eggpk$dist[k] <- min(dist)
}

# Assign a within sample year and doy to the grid data
pk_egg_model1 <- pk_egg_gam1
pk_egg_model2 <- pk_egg_gam2
grid_extent_eggpk$year <- 2014
grid_extent_eggpk$doy <- median(pk_egg$doy)
grid_extent_eggpk$pred1 <- predict(pk_egg_model1, newdata = grid_extent_eggpk)
grid_extent_eggpk$pred2 <- predict(pk_egg_model2, newdata = grid_extent_eggpk)
grid_extent_eggpk$pred <- grid_extent_eggpk$pred1 * grid_extent_eggpk$pred2
grid_extent_eggpk$pred[grid_extent_eggpk$dist > 30000] <- NA

# Phenology
grid_extent_eggpk2 <- data.frame('lon' = rep(-170, 100),
                           'lat' = rep(57, 100),
                           'doy' = seq(min(pk_egg$doy), 
                                       max(pk_egg$doy), 
                                       length = 100),
                           'year' = rep(2014, 100))
grid_extent_eggpk2$pred1 <- predict(pk_egg_model1, newdata = grid_extent_eggpk2)
grid_extent_eggpk2$pred2 <- predict(pk_egg_model2, newdata = grid_extent_eggpk2)
grid_extent_eggpk2$pred <- grid_extent_eggpk2$pred1 * grid_extent_eggpk2$pred2
grid_extent_eggpk2$se1 <- predict(pk_egg_model1, newdata = grid_extent_eggpk2, se = T)[[2]]
grid_extent_eggpk2$se2 <- predict(pk_egg_model2, newdata = grid_extent_eggpk2, se = T)[[2]]
grid_extent_eggpk2$se <- grid_extent_eggpk2$se1 * grid_extent_eggpk2$se2
grid_extent_eggpk2$pred_up <- grid_extent_eggpk2$pred + 1.96 * grid_extent_eggpk2$se
grid_extent_eggpk2$pred_lw <- grid_extent_eggpk2$pred - 1.96 * grid_extent_eggpk2$se

# Plot
windows(width = 8, height = 3.5)
par(mfrow = c(1, 2), 
    mai = c(0.8, 0.9, 0.5, 0.5))
map_phenology(pk_egg, grid_extent_eggpk, grid_extent_eggpk2)
dev.copy(jpeg, here('results/pollock_hindcast', 'pollock_base_egg.jpg'), 
         height = 3.5, width = 8, res = 200, units = 'in')
dev.off()

## Variable coefficient GAM using best models
pk_egg$sst_may <- ncep_temp$may[match(pk_egg$year, ncep_temp$year)]
range(pk_egg$sst_may)

# This will again use the presence GAM from the two stage models
# create variable phenology model by adding doy*sst
month_formula <- formula(log(larvalcatchper10m2 + 1) ~ factor(year) +
                           s(doy) +
                           s(lon, lat) +
                           s(doy, by = sst_may))
pk_egg_month <- gam(month_formula, data = pk_egg[pk_egg$larvalcatchper10m2 > 0, ])
summary(pk_egg_month)

par(mfrow = c(2,2))
plot(pk_egg_month)

par(mfrow = c(2, 2))
gam.check(pk_egg_month)

# Use the original presence model from before
var_ratio_pk_month <- (summary(pk_egg_gam2)$scale - summary(pk_egg_month)$scale) / 
  summary(pk_egg_gam2)$scale
var_ratio_pk_month

# Create variable geography using space*sst
space_formula <- formula(log(larvalcatchper10m2 + 1) ~ factor(year) +
                           s(doy) +
                           s(lon, lat) +
                           s(lon, lat, by = sst_may))
pk_egg_space <- gam(space_formula, data = pk_egg[pk_egg$larvalcatchper10m2 > 0, ])
summary(pk_egg_space)

par(mfrow = c(2,2))
plot(pk_egg_space)

par(mfrow = c(2, 2))
gam.check(pk_egg_space)

# Use the original presence model from before
var_ratio_pk_space <- (summary(pk_egg_gam2)$scale - summary(pk_egg_space)$scale) / 
  summary(pk_egg_gam2)$scale
var_ratio_pk_space

#### Larvae ----
# Negative Binomial
pk_larvae_basenb <- gam(count ~ factor(year) + 
                       s(lon, lat) + 
                       s(doy) +
                       s(roms_salinity) +
                       s(roms_temperature),
                     data = pk_larvae,
                     family = nb(),
                     offset = log(volume_filtered))
summary(pk_larvae_basenb)
# R2: 0.176
# Deviance: 73.6%

pk_larvae_basenb$family$getTheta(TRUE) # 0.3223544

windows()
par(mfrow = c(2, 2))
plot(pk_larvae_basenb, select = 2, main = "DOY")
plot(pk_larvae_basenb, select = 3, main = "Salinity")
plot(pk_larvae_basenb, select = 4, main = "Temperature")
dev.copy(jpeg, here('results/pollock_hindcast', 'pollock_larvae_basenb.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(pk_larvae_basenb)

# quasiPoisson
pk_larvae_basep <- gam(count ~ factor(year) + 
                      s(lon, lat) + 
                      s(doy) +
                      s(roms_salinity) +
                      s(roms_temperature),
                    data = pk_larvae,
                    family = quasipoisson(link = "log"),
                    offset = log(volume_filtered))
summary(pk_larvae_basep)
# R2: 0.502
# Deviance: 76.3%

windows()
par(mfrow = c(2, 2))
plot(pk_larvae_basep, select = 2, main = "DOY")
plot(pk_larvae_basep, select = 3, main = "Salinity")
plot(pk_larvae_basep, select = 4, main = "Temperature")
dev.copy(jpeg, here('results/pollock_hindcast', 'pollock_larvae_basep.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(pk_larvae_basep)

# Tweedie
pk_larvae_baset <- gam(larvalcatchper10m2 ~ factor(year) + 
                         s(lon, lat) + 
                         s(doy) +
                         s(roms_salinity) +
                         s(roms_temperature),
                       data = pk_larvae,
                       family = tw(link = 'log'),
                       method = 'REML')
summary(pk_larvae_baset)
# R2: 0.163
# Deviance: 77.1%

windows()
par(mfrow = c(2, 2))
plot(pk_larvae_baset, select = 2, main = "DOY")
plot(pk_larvae_baset, select = 3, main = "Salinity")
plot(pk_larvae_baset, select = 4, main = "Temperature")
dev.copy(jpeg, here('results/pollock_hindcast', 'pollock_larvae_baset.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(pk_larvae_baset)

# Zero-inflated Poisson (1 stage and 2 stage)
pk_larvae_zip <- gam(count ~ s(year) +  # cannot use year as a factor
                    s(lon, lat) + 
                    s(doy) +
                    s(roms_salinity) +
                    s(roms_temperature),
                  data = pk_larvae,
                  family = ziP(),
                  offset = log(volume_filtered))
summary(pk_larvae_zip)
# Deviance: 100% ?

par(mfrow = c(2, 2))
gam.check(pk_larvae_zip)

plot(pk_larvae_zip, select = 3, main = "Zero-inflated Poisson (1-stage)")

pk_larvae_ziplss <- gam(list(count ~ s(year) +
                            s(lon, lat) +
                            s(doy) +
                            s(roms_salinity) +
                            s(roms_temperature),
                          ~ s(year) +
                            s(lon, lat) +
                            s(doy)),
                     offset = log(volume_filtered),
                     data = pk_larvae,
                     family = ziplss()) #ziplss
summary(pk_larvae_ziplss)
# Deviance: 69.7%

par(mfrow = c(2, 2))
gam.check(pk_larvae_ziplss)

# Two-part binomial and gaussian
# binomial
pk_larvae$presence <- 1 * (pk_larvae$count > 0)
pk_larvae_gam1 <- gam(presence ~ factor(year) +
                     s(doy) +
                     s(lon, lat),
                   data = pk_larvae,
                   family = "binomial")
summary(pk_larvae_gam1)
# R2: 0.587
# Deviance: 51.5%

windows()
plot(pk_larvae_gam1, select = 1, main = "DOY")
dev.copy(jpeg, here('results/pollock_hindcast', 'pollock_larvae_gam1.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(pk_larvae_gam1)

# gaussian
pk_larvae_gam2 <- gam(log(larvalcatchper10m2 + 1) ~ factor(year) +
                        s(doy) +
                        s(lon, lat),
                      data = pk_larvae[pk_larvae$larvalcatchper10m2 > 0, ])
summary(pk_larvae_gam2)
# R2: 0.464
# Deviance: 48.2% 

windows()
par(mfrow = c(2, 2))
plot(pk_larvae_gam2, select = 1, main = "DOY")
plot(pk_larvae_gam2, select = 3, main = "Salinity")
plot(pk_larvae_gam2, select = 4, main = "Temperature")
dev.copy(jpeg, here('results/pollock_hindcast', 'pollock_larvae_gam2.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(pk_larvae_gam2)

# Plot best model
# Plot results of average geography and phenology along with decrease of MSE
# Prediction grid
latd = seq(min(pk_larvae$lat), max(pk_larvae$lat), length.out = nlat)
lond = seq(min(pk_larvae$lon), max(pk_larvae$lon), length.out = nlon)
grid_extent_larvaepk <- expand.grid(lond, latd)
names(grid_extent_larvaepk) <- c('lon', 'lat')

# Calculate distance of each grid point to closest 'positive observation'
grid_extent_larvaepk$dist <- NA
for (k in 1:nrow(grid_extent_larvaepk)) {
  dist <- distance_function(grid_extent_larvaepk$lat[k],
                            grid_extent_larvaepk$lon[k],
                            pk_larvae$lat,
                            pk_larvae$lon)
  grid_extent_larvaepk$dist[k] <- min(dist)
}

# Assign a within sample year and doy to the grid data
pk_larvae_model1 <- pk_larvae_gam1
pk_larvae_model2 <- pk_larvae_gam2
grid_extent_larvaepk$year <- 2014
grid_extent_larvaepk$doy <- median(pk_larvae$doy)
grid_extent_larvaepk$pred1 <- predict(pk_larvae_model1, newdata = grid_extent_larvaepk)
grid_extent_larvaepk$pred2 <- predict(pk_larvae_model2, newdata = grid_extent_larvaepk)
grid_extent_larvaepk$pred <- grid_extent_larvaepk$pred1 * grid_extent_larvaepk$pred2
grid_extent_larvaepk$pred[grid_extent_larvaepk$dist > 30000] <- NA

# Plot phenology
grid_extent_larvaepk2 <- data.frame('lon' = rep(-170, 100),
                           'lat' = rep(57, 100),
                           'doy' = seq(min(pk_larvae$doy), max(pk_larvae$doy), length = 100),
                           'year' = rep(2014, 100))
grid_extent_larvaepk2$pred1 <- predict(pk_larvae_model1, newdata = grid_extent_larvaepk2)
grid_extent_larvaepk2$pred2 <- predict(pk_larvae_model2, newdata = grid_extent_larvaepk2)
grid_extent_larvaepk2$pred <- grid_extent_larvaepk2$pred1 * grid_extent_larvaepk2$pred2
grid_extent_larvaepk2$se1 <- predict(pk_larvae_model1, newdata = grid_extent_larvaepk2, se = T)[[2]]
grid_extent_larvaepk2$se2 <- predict(pk_larvae_model2, newdata = grid_extent_larvaepk2, se = T)[[2]]
grid_extent_larvaepk2$se <- grid_extent_larvaepk2$se1 * grid_extent_larvaepk2$se2
grid_extent_larvaepk2$pred_up <- grid_extent_larvaepk2$pred + 1.96 * grid_extent_larvaepk2$se
grid_extent_larvaepk2$pred_lw <- grid_extent_larvaepk2$pred - 1.96 * grid_extent_larvaepk2$se

# Plot
windows(width = 8, height = 3.5)
par(mfrow = c(1, 2), 
    mai = c(0.8, 0.9, 0.5, 0.5))
map_phenology(pk_larvae, grid_extent_larvaepk, grid_extent_larvaepk2)
dev.copy(jpeg, here('results/pollock_hindcast', 'pollock_base_larvae.jpg'), height = 3.5, width = 8, res = 200, units = 'in')
dev.off()

### Flathead sole ----
#### Eggs ----
hist(fhs_egg$larvalcatchper10m2)

# Negative Binomial
fhs_egg_basenb <- gam(count ~ factor(year) + 
                        s(lon, lat) + 
                        s(doy) +
                        s(roms_salinity) +
                        s(roms_temperature),
                      data = fhs_egg,
                      family = nb(),
                      offset = log(volume_filtered))
summary(fhs_egg_basenb)
# R2: 0.15
# Deviance: 81.1%

fhs_egg_basenb$family$getTheta(TRUE) # 0.3481331

windows()
par(mfrow = c(2, 2))
plot(fhs_egg_basenb, select = 2, main = "DOY")
plot(fhs_egg_basenb, select = 3, main = "Salinity")
plot(fhs_egg_basenb, select = 4, main = "Temperature")
dev.copy(jpeg, here('results/flathead_hindcast', 'flathead_egg_basenb.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(fhs_egg_basenb)

# quasiPoisson
fhs_egg_basep <- gam(count ~ factor(year) + 
                       s(lon, lat) + 
                       s(doy) +
                       s(roms_salinity) +
                       s(roms_temperature),
                     data = fhs_egg,
                     family = quasipoisson(link = "log"),
                     offset = log(volume_filtered))
summary(fhs_egg_basep)
# R2: 0.483
# Deviance: 77.5%

windows()
par(mfrow = c(2, 2))
plot(fhs_egg_basep, select = 2, main = "DOY")
plot(fhs_egg_basep, select = 3, main = "Salinity")
plot(fhs_egg_basep, select = 4, main = "Temperature")
dev.copy(jpeg, here('results/flathead_hindcast', 'flathead_egg_basep.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(fhs_egg_basep)

# Tweedie
fhs_egg_baset <- gam(larvalcatchper10m2 ~ factor(year) + 
                       s(lon, lat) + 
                       s(doy) +
                       s(roms_salinity) +
                       s(roms_temperature),
                     data = fhs_egg,
                     family = tw(link = 'log'),
                     method = 'REML',
                     offset = log(volume_filtered))
summary(fhs_egg_baset)
# R2: 0.285
# Deviance: 77%

windows()
par(mfrow = c(2, 2))
plot(fhs_egg_baset, select = 2, main = "DOY")
plot(fhs_egg_baset, select = 3, main = "Salinity")
plot(fhs_egg_baset, select = 4, main = "Temperature")
dev.copy(jpeg, here('results/flathead_hindcast', 'flathead_egg_baset.jpg'), 
         height = 5, width = 5, units = 'in', res = 200)
dev.off()

par(mfrow = c(2, 2))
gam.check(fhs_egg_baset)

# Zero-inflated Poisson (1 stage and 2 stage)
fhs_egg_zip <- gam(count ~ s(year) +  # cannot use year as a factor
                     s(lon, lat) + 
                     s(doy) +
                     s(roms_salinity) +
                     s(roms_temperature),
                   data = fhs_egg,
                   family = ziP(),
                   offset = log(volume_filtered))
summary(fhs_egg_zip)
# Deviance: 65%

windows()
par(mfrow = c(2, 2))
plot(fhs_egg_zip, select = 3, main = "DOY")
plot(fhs_egg_zip, select = 4, main = "Salinity")
plot(fhs_egg_zip, select = 5, main = "Temperature")
dev.copy(jpeg, here('results/flathead_hindcast', 'flathead_egg_zip.jpg'), 
         height = 5, width = 5, units = 'in', res = 200)
dev.off()

par(mfrow = c(2, 2))
gam.check(fhs_egg_zip)

fhs_egg_ziplss <- gam(list(count ~ s(year) +
                             s(lon, lat) +
                             s(doy) +
                             s(roms_salinity) +
                             s(roms_temperature),
                           ~ s(year) +
                             s(lon, lat) +
                             s(doy)),
                      offset = log(volume_filtered),
                      data = fhs_egg,
                      family = ziplss()) #ziplss
summary(fhs_egg_ziplss)
# Deviance: 58.2%

windows()
par(mfrow = c(2, 2))
plot(fhs_egg_ziplss, select = 3, main = "DOY")
plot(fhs_egg_ziplss, select = 4, main = "Salinity")
plot(fhs_egg_ziplss, select = 5, main = "Temperature")
dev.copy(jpeg, here('results/flathead_hindcast', 'flathead_egg_ziplss.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(fhs_egg_ziplss)

# Two-part binomial and gaussian
# binomial
fhs_egg$presence <- 1 * (fhs_egg$count > 0)
fhs_egg_gam1 <- gam(presence ~ factor(year) +
                      s(doy) +
                      s(lon, lat, k = 4),
                    data = fhs_egg,
                    family = "binomial",
                    offset = log(volume_filtered))
summary(fhs_egg_gam1)
# R2: 0.535
# Deviance: 53.9%

windows()
plot(fhs_egg_gam1, select = 1, main = "DOY")
dev.copy(jpeg, here('results/flathead_hindcast', 'flathead_egg_gam1.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(fhs_egg_gam1)

# gaussian
fhs_egg_gam2 <- gam(log(larvalcatchper10m2 + 1) ~ factor(year) +
                      s(doy) +
                      s(lon, lat, k = 4),
                    offset = log(volume_filtered),
                    data = fhs_egg[fhs_egg$larvalcatchper10m2 > 0, ])
summary(fhs_egg_gam2)
# R2: 0.116
# Deviance: 19.7% 

windows()
par(mfrow = c(2, 2))
plot(fhs_egg_gam2, select = 1, main = "DOY")
plot(fhs_egg_gam2, select = 3, main = "Salinity")
plot(fhs_egg_gam2, select = 4, main = "Temperature")
dev.copy(jpeg, here('results/flathead_hindcast', 'flathead_egg_gam2.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(fhs_egg_gam2)


# Plot best model
# Plot results of average geography and phenology along with decrease of MSE
# Prediction grid
nlat = 80
nlon = 120
latd = seq(min(fhs_egg$lat), max(fhs_egg$lat), length.out = nlat)
lond = seq(min(fhs_egg$lon), max(fhs_egg$lon), length.out = nlon)
grid_extent_eggfhs <- expand.grid(lond, latd)
names(grid_extent_eggfhs) <- c('lon', 'lat')

# Calculate distance of each grid point to closest 'positive observation'
grid_extent_eggfhs$dist <- NA
for (k in 1:nrow(grid_extent_eggfhs)) {
  dist <- distance_function(grid_extent_eggfhs$lat[k],
                            grid_extent_eggfhs$lon[k],
                            fhs_egg$lat,
                            fhs_egg$lon)
  grid_extent_eggfhs$dist[k] <- min(dist)
}

# Assign a within sample year and doy to the grid data
fhs_egg_model1 <- fhs_egg_gam1
fhs_egg_model2 <- fhs_egg_gam2
grid_extent_eggfhs$year <- 2014
grid_extent_eggfhs$doy <- median(fhs_egg$doy)
grid_extent_eggfhs$pred1 <- predict(fhs_egg_model1, newdata = grid_extent_eggfhs)
grid_extent_eggfhs$pred2 <- predict(fhs_egg_model2, newdata = grid_extent_eggfhs)
grid_extent_eggfhs$pred <- grid_extent_eggfhs$pred1 * grid_extent_eggfhs$pred2
grid_extent_eggfhs$pred[grid_extent_eggfhs$dist > 30000] <- NA

# Plot phenology
grid_extent_eggfhs2 <- data.frame('lon' = rep(-170, 100),
                                  'lat' = rep(57, 100),
                                  'doy' = seq(min(fhs_egg$doy), 
                                              max(fhs_egg$doy), 
                                              length = 100),
                                  'year' = rep(2014, 100))
grid_extent_eggfhs2$pred1 <- predict(fhs_egg_model1, newdata = grid_extent_eggfhs2)
grid_extent_eggfhs2$pred2 <- predict(fhs_egg_model2, newdata = grid_extent_eggfhs2)
grid_extent_eggfhs2$pred <- grid_extent_eggfhs2$pred1 * grid_extent_eggfhs2$pred2
grid_extent_eggfhs2$se1 <- predict(fhs_egg_model1, newdata = grid_extent_eggfhs2, se = T)[[2]]
grid_extent_eggfhs2$se2 <- predict(fhs_egg_model2, newdata = grid_extent_eggfhs2, se = T)[[2]]
grid_extent_eggfhs2$se <- grid_extent_eggfhs2$se1 * grid_extent_eggfhs2$se2
grid_extent_eggfhs2$pred_up <- grid_extent_eggfhs2$pred + 1.96 * grid_extent_eggfhs2$se
grid_extent_eggfhs2$pred_lw <- grid_extent_eggfhs2$pred - 1.96 * grid_extent_eggfhs2$se

# Plot
windows(width = 8, height = 3.5)
par(mfrow = c(1, 2), 
    mai = c(0.8, 0.9, 0.5, 0.5))
map_phenology(fhs_egg, grid_extent_eggfhs, grid_extent_eggfhs2)
dev.copy(jpeg, here('results/flathead_hindcast', 'flathead_base_egg.jpg'), 
         height = 3.5, width = 8, res = 200, units = 'in')
dev.off()

#### Larvae ----
# Negative Binomial
fhs_larvae_basenb <- gam(count ~ factor(year) + 
                           s(lon, lat) + 
                           s(doy) +
                           s(roms_salinity) +
                           s(roms_temperature),
                         data = fhs_larvae,
                         family = nb(),
                         offset = log(volume_filtered))
summary(fhs_larvae_basenb)
# R2: 0.246
# Deviance: 76.8%

fhs_larvae_basenb$family$getTheta(TRUE) # 0.3615238

windows()
par(mfrow = c(2, 2))
plot(fhs_larvae_basenb, select = 2, main = "DOY")
plot(fhs_larvae_basenb, select = 3, main = "Salinity")
plot(fhs_larvae_basenb, select = 4, main = "Temperature")
dev.copy(jpeg, here('results/flathead_hindcast', 'flathead_larvae_basenb.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(fhs_larvae_basenb)

# quasiPoisson
fhs_larvae_basep <- gam(count ~ factor(year) + 
                          s(lon, lat) + 
                          s(doy) +
                          s(roms_salinity) +
                          s(roms_temperature),
                        data = fhs_larvae,
                        family = quasipoisson(link = "log"),
                        offset = log(volume_filtered))
summary(fhs_larvae_basep)
# R2: 0.493
# Deviance: 75%

windows()
par(mfrow = c(2, 2))
plot(fhs_larvae_basep, select = 2, main = "DOY")
plot(fhs_larvae_basep, select = 3, main = "Salinity")
plot(fhs_larvae_basep, select = 4, main = "Temperature")
dev.copy(jpeg, here('results/flathead_hindcast', 'flathead_larvae_basep.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(fhs_larvae_basep)

# Tweedie
fhs_larvae_baset <- gam(larvalcatchper10m2 ~ factor(year) + 
                          s(lon, lat) + 
                          s(doy) +
                          s(roms_salinity) +
                          s(roms_temperature),
                        data = fhs_larvae,
                        family = tw(link = 'log'),
                        method = 'REML')
summary(fhs_larvae_baset)
# R2: 0.293
# Deviance: 74.3%

windows()
par(mfrow = c(2, 2))
plot(fhs_larvae_baset, select = 2, main = "DOY")
plot(fhs_larvae_baset, select = 3, main = "Salinity")
plot(fhs_larvae_baset, select = 4, main = "Temperature")
dev.copy(jpeg, here('results/flathead_hindcast', 'flathead_larvae_baset.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(fhs_larvae_baset)

# Zero-inflated Poisson (1 stage and 2 stage)
fhs_larvae_zip <- gam(count ~ s(year) +  # cannot use year as a factor
                        s(lon, lat) + 
                        s(doy) +
                        s(roms_salinity) +
                        s(roms_temperature),
                      data = fhs_larvae,
                      family = ziP(),
                      offset = log(volume_filtered))
summary(fhs_larvae_zip)
# Deviance: 58.7% 

windows()
par(mfrow = c(2, 2))
plot(fhs_larvae_zip, select = 3, main = "DOY")
plot(fhs_larvae_zip, select = 4, main = "Salinity")
plot(fhs_larvae_zip, select = 5, main = "Temperature")
dev.copy(jpeg, here('results/flathead_hindcast', 'flathead_larvae_zip.jpg'), 
         height = 5, width = 5, units = 'in', res = 200)
dev.off()

par(mfrow = c(2, 2))
gam.check(fhs_larvae_zip)

fhs_larvae_ziplss <- gam(list(count ~ s(year) +
                                s(lon, lat) +
                                s(doy) +
                                s(roms_salinity) +
                                s(roms_temperature),
                              ~ s(year) +
                                s(lon, lat) +
                                s(doy)),
                         offset = log(volume_filtered),
                         data = fhs_larvae,
                         family = ziplss()) #ziplss
summary(fhs_larvae_ziplss)
# Deviance: 57.5%

windows()
par(mfrow = c(2, 2))
plot(fhs_larvae_ziplss, select = 3, main = "DOY")
plot(fhs_larvae_ziplss, select = 4, main = "Salinity")
plot(fhs_larvae_ziplss, select = 5, main = "Temperature")
dev.copy(jpeg, here('results/flathead_hindcast', 'flathead_larvae_ziplss.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(fhs_larvae_ziplss)

# Two-part binomial and gaussian
# binomial
fhs_larvae$presence <- 1 * (fhs_larvae$count > 0)
fhs_larvae_gam1 <- gam(presence ~ factor(year) +
                         s(doy) +
                         s(lon, lat),
                       data = fhs_larvae,
                       family = "binomial")
summary(fhs_larvae_gam1)
# R2: 0.514
# Deviance: 49.4%

windows()
plot(fhs_larvae_gam1, select = 1, main = "DOY")
dev.copy(jpeg, here('results/flathead_hindcast', 'flathead_larvae_gam1.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(fhs_larvae_gam1)

# gaussian
fhs_larvae_gam2 <- gam(log(larvalcatchper10m2 + 1) ~ factor(year) +
                         s(doy) +
                         s(lon, lat),
                       data = fhs_larvae[fhs_larvae$larvalcatchper10m2 > 0, ])
summary(fhs_larvae_gam2)
# R2: 0.383
# Deviance: 41.7% 

windows()
par(mfrow = c(2, 2))
plot(fhs_larvae_gam2, select = 1, main = "DOY")
plot(fhs_larvae_gam2, select = 3, main = "Salinity")
plot(fhs_larvae_gam2, select = 4, main = "Temperature")
dev.copy(jpeg, here('results/flathead_hindcast', 'flathead_larvae_gam2.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(fhs_larvae_gam2)

# Plot best model
# Plot results of average geography and phenology along with decrease of MSE
# Prediction grid
latd = seq(min(fhs_larvae$lat), max(fhs_larvae$lat), length.out = nlat)
lond = seq(min(fhs_larvae$lon), max(fhs_larvae$lon), length.out = nlon)
grid_extent_larvaefhs <- expand.grid(lond, latd)
names(grid_extent_larvaefhs) <- c('lon', 'lat')

# Calculate distance of each grid point to closest 'positive observation'
grid_extent_larvaefhs$dist <- NA
for (k in 1:nrow(grid_extent_larvaefhs)) {
  dist <- distance_function(grid_extent_larvaefhs$lat[k],
                            grid_extent_larvaefhs$lon[k],
                            fhs_larvae$lat,
                            fhs_larvae$lon)
  grid_extent_larvaefhs$dist[k] <- min(dist)
}

# Assign a within sample year and doy to the grid data
fhs_larvae_model1 <- fhs_larvae_gam1
fhs_larvae_model2 <- fhs_larvae_gam2
grid_extent_larvaefhs$year <- 2014
grid_extent_larvaefhs$doy <- median(fhs_larvae$doy)
grid_extent_larvaefhs$pred1 <- predict(fhs_larvae_model1, newdata = grid_extent_larvaefhs)
grid_extent_larvaefhs$pred2 <- predict(fhs_larvae_model2, newdata = grid_extent_larvaefhs)
grid_extent_larvaefhs$pred <- grid_extent_larvaefhs$pred1 * grid_extent_larvaefhs$pred2
grid_extent_larvaefhs$pred[grid_extent_larvaefhs$dist > 30000] <- NA

# Plot phenology
grid_extent_larvaefhs2 <- data.frame('lon' = rep(-170, 100),
                                     'lat' = rep(57, 100),
                                     'doy' = seq(min(fhs_larvae$doy), max(fhs_larvae$doy), length = 100),
                                     'year' = rep(2014, 100))
grid_extent_larvaefhs2$pred1 <- predict(fhs_larvae_model1, newdata = grid_extent_larvaefhs2)
grid_extent_larvaefhs2$pred2 <- predict(fhs_larvae_model2, newdata = grid_extent_larvaefhs2)
grid_extent_larvaefhs2$pred <- grid_extent_larvaefhs2$pred1 * grid_extent_larvaefhs2$pred2
grid_extent_larvaefhs2$se1 <- predict(fhs_larvae_model1, newdata = grid_extent_larvaefhs2, se = T)[[2]]
grid_extent_larvaefhs2$se2 <- predict(fhs_larvae_model2, newdata = grid_extent_larvaefhs2, se = T)[[2]]
grid_extent_larvaefhs2$se <- grid_extent_larvaefhs2$se1 * grid_extent_larvaefhs2$se2 # need to change this, can't just be multiplied
grid_extent_larvaefhs2$pred_up <- grid_extent_larvaefhs2$pred + 1.96 * grid_extent_larvaefhs2$se
grid_extent_larvaefhs2$pred_lw <- grid_extent_larvaefhs2$pred - 1.96 * grid_extent_larvaefhs2$se

# Plot
windows(width = 8, height = 3.5)
par(mfrow = c(1, 2), 
    mai = c(0.8, 0.9, 0.5, 0.5))
map_phenology(fhs_larvae, grid_extent_larvaefhs, grid_extent_larvaefhs2)
dev.copy(jpeg, here('results/flathead_hindcast', 'flathead_base_larvae.jpg'), height = 3.5, width = 8, res = 200, units = 'in')
dev.off()

### Alaska Plaice ----
#### Eggs ----
hist(akp_egg$larvalcatchper10m2)

# Negative Binomial
akp_egg_basenb <- gam(count ~ factor(year) + 
                        s(lon, lat) + 
                        s(doy) +
                        s(roms_salinity) +
                        s(roms_temperature),
                      data = akp_egg,
                      family = nb(),
                      offset = log(volume_filtered))
summary(akp_egg_basenb)
# R2: -6.93
# Deviance: 85.2%

akp_egg_basenb$family$getTheta(TRUE) # 0.3602255

windows()
par(mfrow = c(2, 2))
plot(akp_egg_basenb, select = 2, main = "DOY")
plot(akp_egg_basenb, select = 3, main = "Salinity")
plot(akp_egg_basenb, select = 4, main = "Temperature")
dev.copy(jpeg, here('results/plaice_hindcast', 'plaice_egg_basenb.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(akp_egg_basenb)

# quasiPoisson
akp_egg_basep <- gam(count ~ factor(year) + 
                       s(lon, lat) + 
                       s(doy) +
                       s(roms_salinity) +
                       s(roms_temperature),
                     data = akp_egg,
                     family = quasipoisson(link = "log"),
                     offset = log(volume_filtered))
summary(akp_egg_basep)
# R2: 0.537
# Deviance: 83.8%

windows()
par(mfrow = c(2, 2))
plot(akp_egg_basep, select = 2, main = "DOY")
plot(akp_egg_basep, select = 3, main = "Salinity")
plot(akp_egg_basep, select = 4, main = "Temperature")
dev.copy(jpeg, here('results/plaice_hindcast', 'plaice_egg_basep.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(akp_egg_basep)

# Tweedie
akp_egg_baset <- gam(larvalcatchper10m2 ~ factor(year) + 
                       s(lon, lat) + 
                       s(doy) +
                       s(roms_salinity) +
                       s(roms_temperature),
                     data = akp_egg,
                     family = tw(link = 'log'),
                     method = 'REML',
                     offset = log(volume_filtered))
summary(akp_egg_baset)
# R2: -13.1
# Deviance: 85.7%

windows()
par(mfrow = c(2, 2))
plot(akp_egg_baset, select = 2, main = "DOY")
plot(akp_egg_baset, select = 3, main = "Salinity")
plot(akp_egg_baset, select = 4, main = "Temperature")
dev.copy(jpeg, here('results/plaice_hindcast', 'plaice_egg_baset.jpg'), 
         height = 5, width = 5, units = 'in', res = 200)
dev.off()

par(mfrow = c(2, 2))
gam.check(akp_egg_baset)

# Zero-inflated Poisson (1 stage and 2 stage)
akp_egg_zip <- gam(count ~ s(year) +  # cannot use year as a factor
                     s(lon, lat) + 
                     s(doy) +
                     s(roms_salinity) +
                     s(roms_temperature),
                   data = akp_egg,
                   family = ziP(),
                   offset = log(volume_filtered))
summary(akp_egg_zip)
# Deviance: 69.7%

windows()
par(mfrow = c(2, 2))
plot(akp_egg_zip, select = 3, main = "DOY")
plot(akp_egg_zip, select = 4, main = "Salinity")
plot(akp_egg_zip, select = 5, main = "Temperature")
dev.copy(jpeg, here('results/plaice_hindcast', 'plaice_egg_zip.jpg'), 
         height = 5, width = 5, units = 'in', res = 200)
dev.off()

par(mfrow = c(2, 2))
gam.check(akp_egg_zip)

akp_egg_ziplss <- gam(list(count ~ s(year) +
                             s(lon, lat) +
                             s(doy) +
                             s(roms_salinity) +
                             s(roms_temperature),
                           ~ s(year) +
                             s(lon, lat) +
                             s(doy)),
                      offset = log(volume_filtered),
                      data = akp_egg,
                      family = ziplss()) #ziplss
summary(akp_egg_ziplss)
# Deviance: 68.2%

windows()
par(mfrow = c(2, 2))
plot(akp_egg_ziplss, select = 3, main = "DOY")
plot(akp_egg_ziplss, select = 4, main = "Salinity")
plot(akp_egg_ziplss, select = 5, main = "Temperature")
dev.copy(jpeg, here('results/plaice_hindcast', 'plaice_egg_ziplss.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(akp_egg_ziplss)

# Two-part binomial and gaussian
# binomial
akp_egg$presence <- 1 * (akp_egg$count > 0)
akp_egg_gam1 <- gam(presence ~ factor(year) +
                      s(doy) +
                      s(lon, lat),
                    data = akp_egg,
                    family = "binomial",
                    offset = log(volume_filtered))
summary(akp_egg_gam1)
# R2: 0.594
# Deviance: 61.3%

windows()
plot(akp_egg_gam1, select = 1, main = "DOY")
dev.copy(jpeg, here('results/plaice_hindcast', 'plaice_egg_gam1.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(akp_egg_gam1)

# gaussian
akp_egg_gam2 <- gam(log(larvalcatchper10m2 + 1) ~ factor(year) +
                      s(doy) +
                      s(lon, lat),
                    offset = log(volume_filtered),
                    data = akp_egg[akp_egg$larvalcatchper10m2 > 0, ])
summary(akp_egg_gam2)
# R2: 0.278
# Deviance: 46.3% 

windows()
par(mfrow = c(2, 2))
plot(akp_egg_gam2, select = 1, main = "DOY")
plot(akp_egg_gam2, select = 3, main = "Salinity")
plot(akp_egg_gam2, select = 4, main = "Temperature")
dev.copy(jpeg, here('results/plaice_hindcast', 'plaice_egg_gam2.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(akp_egg_gam2)


# Plot best model
# Plot results of average geography and phenology along with decrease of MSE
# Prediction grid
nlat = 80
nlon = 120
latd = seq(min(akp_egg$lat), max(akp_egg$lat), length.out = nlat)
lond = seq(min(akp_egg$lon), max(akp_egg$lon), length.out = nlon)
grid_extent_eggakp <- expand.grid(lond, latd)
names(grid_extent_eggakp) <- c('lon', 'lat')

# Calculate distance of each grid point to closest 'positive observation'
grid_extent_eggakp$dist <- NA
for (k in 1:nrow(grid_extent_eggakp)) {
  dist <- distance_function(grid_extent_eggakp$lat[k],
                            grid_extent_eggakp$lon[k],
                            akp_egg$lat,
                            akp_egg$lon)
  grid_extent_eggakp$dist[k] <- min(dist)
}

# Assign a within sample year and doy to the grid data
akp_egg_model1 <- akp_egg_gam1
akp_egg_model2 <- akp_egg_gam2
grid_extent_eggakp$year <- 2014
grid_extent_eggakp$doy <- median(akp_egg$doy)
grid_extent_eggakp$pred1 <- predict(akp_egg_model1, newdata = grid_extent_eggakp)
grid_extent_eggakp$pred2 <- predict(akp_egg_model2, newdata = grid_extent_eggakp)
grid_extent_eggakp$pred <- grid_extent_eggakp$pred1 * grid_extent_eggakp$pred2
grid_extent_eggakp$pred[grid_extent_eggakp$dist > 30000] <- NA

# Plot phenology
grid_extent_eggakp2 <- data.frame('lon' = rep(-170, 100),
                                  'lat' = rep(57, 100),
                                  'doy' = seq(min(akp_egg$doy), 
                                              max(akp_egg$doy), 
                                              length = 100),
                                  'year' = rep(2014, 100))
grid_extent_eggakp2$pred1 <- predict(akp_egg_model1, newdata = grid_extent_eggakp2)
grid_extent_eggakp2$pred2 <- predict(akp_egg_model2, newdata = grid_extent_eggakp2)
grid_extent_eggakp2$pred <- grid_extent_eggakp2$pred1 * grid_extent_eggakp2$pred2
grid_extent_eggakp2$se1 <- predict(akp_egg_model1, newdata = grid_extent_eggakp2, se = T)[[2]]
grid_extent_eggakp2$se2 <- predict(akp_egg_model2, newdata = grid_extent_eggakp2, se = T)[[2]]
grid_extent_eggakp2$se <- grid_extent_eggakp2$se1 * grid_extent_eggakp2$se2
grid_extent_eggakp2$pred_up <- grid_extent_eggakp2$pred + 1.96 * grid_extent_eggakp2$se
grid_extent_eggakp2$pred_lw <- grid_extent_eggakp2$pred - 1.96 * grid_extent_eggakp2$se

# Plot
windows(width = 8, height = 3.5)
par(mfrow = c(1, 2), 
    mai = c(0.8, 0.9, 0.5, 0.5))
map_phenology(akp_egg, grid_extent_eggakp, grid_extent_eggakp2)
dev.copy(jpeg, here('results/plaice_hindcast', 'plaice_base_egg.jpg'), 
         height = 3.5, width = 8, res = 200, units = 'in')
dev.off()

#### Larvae ----
# Negative Binomial
akp_larvae_basenb <- gam(count ~ factor(year) + 
                           s(lon, lat) + 
                           s(doy) +
                           s(roms_salinity) +
                           s(roms_temperature),
                         data = akp_larvae,
                         family = nb(),
                         offset = log(volume_filtered))
summary(akp_larvae_basenb)
# R2: 0.246
# Deviance: 87.5%

akp_larvae_basenb$family$getTheta(TRUE) # 0.3615238

windows()
par(mfrow = c(2, 2))
plot(akp_larvae_basenb, select = 2, main = "DOY")
plot(akp_larvae_basenb, select = 3, main = "Salinity")
plot(akp_larvae_basenb, select = 4, main = "Temperature")
dev.copy(jpeg, here('results/plaice_hindcast', 'plaice_larvae_basenb.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(akp_larvae_basenb)

# quasiPoisson
akp_larvae_basep <- gam(count ~ factor(year) + 
                          s(lon, lat) + 
                          s(doy) +
                          s(roms_salinity) +
                          s(roms_temperature),
                        data = akp_larvae,
                        family = quasipoisson(link = "log"),
                        offset = log(volume_filtered))
summary(akp_larvae_basep)
# R2: 0.659
# Deviance: 86.9%

windows()
par(mfrow = c(2, 2))
plot(akp_larvae_basep, select = 2, main = "DOY")
plot(akp_larvae_basep, select = 3, main = "Salinity")
plot(akp_larvae_basep, select = 4, main = "Temperature")
dev.copy(jpeg, here('results/plaice_hindcast', 'plaice_larvae_basep.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(akp_larvae_basep)

# Tweedie
akp_larvae_baset <- gam(larvalcatchper10m2 ~ factor(year) + 
                          s(lon, lat) + 
                          s(doy) +
                          s(roms_salinity) +
                          s(roms_temperature),
                        data = akp_larvae,
                        family = tw(link = 'log'),
                        method = 'REML')
summary(akp_larvae_baset)
# R2: 0.519
# Deviance: 82.4%

windows()
par(mfrow = c(2, 2))
plot(akp_larvae_baset, select = 2, main = "DOY")
plot(akp_larvae_baset, select = 3, main = "Salinity")
plot(akp_larvae_baset, select = 4, main = "Temperature")
dev.copy(jpeg, here('results/plaice_hindcast', 'plaice_larvae_baset.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(akp_larvae_baset)

# Zero-inflated Poisson (1 stage and 2 stage)
akp_larvae_zip <- gam(count ~ s(year) +  # cannot use year as a factor
                        s(lon, lat) + 
                        s(doy) +
                        s(roms_salinity) +
                        s(roms_temperature),
                      data = akp_larvae,
                      family = ziP(),
                      offset = log(volume_filtered))
summary(akp_larvae_zip)
# Deviance: 78.6% 

windows()
par(mfrow = c(2, 2))
plot(akp_larvae_zip, select = 3, main = "DOY")
plot(akp_larvae_zip, select = 4, main = "Salinity")
plot(akp_larvae_zip, select = 5, main = "Temperature")
dev.copy(jpeg, here('results/plaice_hindcast', 'plaice_larvae_zip.jpg'), 
         height = 5, width = 5, units = 'in', res = 200)
dev.off()

par(mfrow = c(2, 2))
gam.check(akp_larvae_zip)

akp_larvae_ziplss <- gam(list(count ~ s(year) +
                                s(lon, lat) +
                                s(doy) +
                                s(roms_salinity) +
                                s(roms_temperature),
                              ~ s(year) +
                                s(lon, lat) +
                                s(doy)),
                         offset = log(volume_filtered),
                         data = akp_larvae,
                         family = ziplss()) #ziplss
summary(akp_larvae_ziplss)
# Deviance: 71%

windows()
par(mfrow = c(2, 2))
plot(akp_larvae_ziplss, select = 3, main = "DOY")
plot(akp_larvae_ziplss, select = 4, main = "Salinity")
plot(akp_larvae_ziplss, select = 5, main = "Temperature")
dev.copy(jpeg, here('results/plaice_hindcast', 'plaice_larvae_ziplss.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(akp_larvae_ziplss)

# Two-part binomial and gaussian
# binomial
akp_larvae$presence <- 1 * (akp_larvae$count > 0)
akp_larvae_gam1 <- gam(presence ~ factor(year) +
                         s(doy) +
                         s(lon, lat),
                       data = akp_larvae,
                       family = "binomial")
summary(akp_larvae_gam1)
# R2: 0.432
# Deviance: 49.4%

windows()
plot(akp_larvae_gam1, select = 1, main = "DOY")
dev.copy(jpeg, here('results/plaice_hindcast', 'plaice_larvae_gam1.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(akp_larvae_gam1)

# gaussian
akp_larvae_gam2 <- gam(log(larvalcatchper10m2 + 1) ~ factor(year) +
                         s(doy) +
                         s(lon, lat),
                       data = akp_larvae[akp_larvae$larvalcatchper10m2 > 0, ])
summary(akp_larvae_gam2)
# R2: 0.512
# Deviance: 58.6% 

windows()
par(mfrow = c(2, 2))
plot(akp_larvae_gam2, select = 1, main = "DOY")
plot(akp_larvae_gam2, select = 3, main = "Salinity")
plot(akp_larvae_gam2, select = 4, main = "Temperature")
dev.copy(jpeg, here('results/plaice_hindcast', 'plaice_larvae_gam2.jpg'), 
         height = 5, width = 5, units = 'in', res = 200 )
dev.off()

par(mfrow = c(2, 2))
gam.check(akp_larvae_gam2)

# Plot best model
# Plot results of average geography and phenology along with decrease of MSE
# Prediction grid
latd = seq(min(akp_larvae$lat), max(akp_larvae$lat), length.out = nlat)
lond = seq(min(akp_larvae$lon), max(akp_larvae$lon), length.out = nlon)
grid_extent_larvaeakp <- expand.grid(lond, latd)
names(grid_extent_larvaeakp) <- c('lon', 'lat')

# Calculate distance of each grid point to closest 'positive observation'
grid_extent_larvaeakp$dist <- NA
for (k in 1:nrow(grid_extent_larvaeakp)) {
  dist <- distance_function(grid_extent_larvaeakp$lat[k],
                            grid_extent_larvaeakp$lon[k],
                            akp_larvae$lat,
                            akp_larvae$lon)
  grid_extent_larvaeakp$dist[k] <- min(dist)
}

# Assign a within sample year and doy to the grid data
akp_larvae_model1 <- akp_larvae_gam1
akp_larvae_model2 <- akp_larvae_gam2
grid_extent_larvaeakp$year <- 2014
grid_extent_larvaeakp$doy <- median(akp_larvae$doy)
grid_extent_larvaeakp$pred1 <- predict(akp_larvae_model1, newdata = grid_extent_larvaeakp)
grid_extent_larvaeakp$pred2 <- predict(akp_larvae_model2, newdata = grid_extent_larvaeakp)
grid_extent_larvaeakp$pred <- grid_extent_larvaeakp$pred1 * grid_extent_larvaeakp$pred2
grid_extent_larvaeakp$pred[grid_extent_larvaeakp$dist > 30000] <- NA

# Plot phenology
grid_extent_larvaeakp2 <- data.frame('lon' = rep(-170, 100),
                                     'lat' = rep(57, 100),
                                     'doy' = seq(min(akp_larvae$doy), max(akp_larvae$doy), length = 100),
                                     'year' = rep(2014, 100))
grid_extent_larvaeakp2$pred1 <- predict(akp_larvae_model1, newdata = grid_extent_larvaeakp2)
grid_extent_larvaeakp2$pred2 <- predict(akp_larvae_model2, newdata = grid_extent_larvaeakp2)
grid_extent_larvaeakp2$pred <- grid_extent_larvaeakp2$pred1 * grid_extent_larvaeakp2$pred2
grid_extent_larvaeakp2$se1 <- predict(akp_larvae_model1, newdata = grid_extent_larvaeakp2, se = T)[[2]]
grid_extent_larvaeakp2$se2 <- predict(akp_larvae_model2, newdata = grid_extent_larvaeakp2, se = T)[[2]]
grid_extent_larvaeakp2$se <- grid_extent_larvaeakp2$se1 * grid_extent_larvaeakp2$se2 # need to change this, can't just be multiplied
grid_extent_larvaeakp2$pred_up <- grid_extent_larvaeakp2$pred + 1.96 * grid_extent_larvaeakp2$se
grid_extent_larvaeakp2$pred_lw <- grid_extent_larvaeakp2$pred - 1.96 * grid_extent_larvaeakp2$se

# Plot
windows(width = 8, height = 3.5)
par(mfrow = c(1, 2), 
    mai = c(0.8, 0.9, 0.5, 0.5))
map_phenology(akp_larvae, grid_extent_larvaeakp, grid_extent_larvaeakp2)
dev.copy(jpeg, here('results/plaice_hindcast', 'plaice_base_larvae.jpg'), height = 3.5, width = 8, res = 200, units = 'in')
dev.off()