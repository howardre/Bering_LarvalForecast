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
library(ggplot2)
source(here('code/functions', 'vis_gam_COLORS.R'))

### Load fish data ----
yfs_egg <- readRDS(here('data', 'yfs_egg.rds'))
yfs_larvae <- readRDS(here('data', 'yfs_larvae.rds'))
akp_egg <- readRDS(here('data', 'akp_egg.rds'))
akp_larvae <- readRDS(here('data', 'akp_larvae.rds'))
fhs_egg <- readRDS(here('data', 'fhs_egg.rds'))
fhs_larvae <- readRDS(here('data', 'fhs_larvae.rds'))
pk_egg <- readRDS(here('data', 'pk_egg.rds'))
pk_larvae <- readRDS(here('data', 'pk_larvae.rds'))

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

### Remove outliers ----
yfs_egg <- filter(yfs_egg, larvalcatchper10m2 < 330)
yfs_larvae <- filter(yfs_larvae, larvalcatchper10m2 < 30000,
                     roms_salinity > 25)
akp_egg <- filter(akp_egg, larvalcatchper10m2 < 4000)
akp_larvae <- filter(akp_larvae, larvalcatchper10m2 < 500,
                     roms_salinity > 29)
fhs_egg <- filter(fhs_egg, larvalcatchper10m2 < 1400,
                  roms_salinity > 29)
fhs_larvae <- filter(fhs_larvae, larvalcatchper10m2 < 9000)
pk_egg <- filter(pk_egg, larvalcatchper10m2 < 350000,
                 roms_salinity > 29)
pk_larvae <- filter(pk_larvae, larvalcatchper10m2 < 100000,
                    roms_salinity > 29)

# map of data
# Create bathymetry dataset
BS_bathy <- getNOAA.bathy(lon1= -170, 
                          lon2= -156, 
                          lat1= 60.1, 
                          lat2= 54.2, 
                          resolution=1)
blues <- c("lightsteelblue4", 
           "lightsteelblue3", 
           "lightsteelblue2", 
           "lightsteelblue1")
greys <- c(grey(0.6), 
           grey(0.93), 
           grey(0.99))
plot.bathy(
  BS_bathy,
  image = T,
  axes = T,
  lwd = 0.03,
  land = T,
  n = 0,
  bpal = list(c(0, max(BS_bathy), greys), c(min(BS_bathy), 0, blues)),
  ylim = c(54.2, 60.1),
  xlim = c(-170,-156),
  ylab = "Latitude °N",
  xlab = "Longitude °W",
  cex.lab = 2,
  cex.main = 2.5,
  cex.axis = 1.8)
points(yfs_egg$lon[yfs_egg$larvalcatchper10m2 == 0],
       yfs_egg$lat[yfs_egg$larvalcatchper10m2 == 0],
       pch = 4,
       col = 'darkgray',
       cex = 1)
points(yfs_egg$lon[yfs_egg$larvalcatchper10m2 > 0],
       yfs_egg$lat[yfs_egg$larvalcatchper10m2 > 0],
       pch = 16,
       col = 'black',
       cex = 1.3)
# doesn't seem to have anything weird going on
# was fixed with the distance function during cleaning, don't need to do this for all species


### Yellowfin sole ----
#### Eggs ----
hist(yfs_egg$larvalcatchper10m2)
yfs_egg_gam1 <- gam(larvalcatchper10m2 ~ factor(year) +
                      s(lon, lat) +
                      s(doy) +
                      s(roms_temperature) +
                      s(roms_salinity),
                    data = yfs_egg)
summary(yfs_egg_gam1)
gam.check(yfs_egg_gam1)
plot(yfs_egg_gam1)

yfs_egg$counts <- as.integer(yfs_egg$larvalcatchper10m2)
yfs_egg_gam2 <- gam(counts ~ s(year, k = 5) +
                      s(lon, lat) +
                      s(doy, k = 5) +
                      s(roms_temperature, k = 5) +
                      s(roms_salinity, k = 5),
                    data = yfs_egg,
                    family = ziP()) # ZIP family
summary(yfs_egg_gam2)
par(mfrow = c(2, 2))
gam.check(yfs_egg_gam2)
# not great unconstrained
# possible overfitting

yfs_egg_gam3 <- gam(list(counts ~ s(year, k = 5) +
                           s(lon, lat) +
                           s(doy, k = 5) +
                           s(roms_temperature, k = 5) +
                           s(roms_salinity, k = 5),
                         ~ s(year, k = 5) +
                           s(lon, lat) +
                           s(doy, k = 5)),
                    data = yfs_egg,
                    family = ziplss()) #ziplss
summary(yfs_egg_gam3)
par(mfrow = c(2, 2))
gam.check(yfs_egg_gam3)

#### Larvae ----
hist(yfs_larvae$larvalcatchper10m2)
yfs_larvae$counts <- as.integer(yfs_larvae$larvalcatchper10m2)
yfs_larvae_gam1 <- gam(larvalcatchper10m2 ~ factor(year) +
                      s(lon, lat) +
                      s(doy) +
                      s(roms_temperature) +
                      s(roms_salinity),
                    data = yfs_larvae)
summary(yfs_larvae_gam1)
gam.check(yfs_larvae_gam1)
plot(yfs_larvae_gam1)

yfs_larvae_gam2 <- gam(counts ~ s(year, k = 5) +
                         s(lon, lat) +
                         s(doy, k = 5) +
                         s(roms_temperature, k = 5) +
                         s(roms_salinity, k = 5),
                       data = yfs_larvae,
                       family = ziP())
summary(yfs_larvae_gam2)
gam.check(yfs_larvae_gam2)
plot(yfs_larvae_gam2)

yfs_larvae_gam3 <- gam(list(counts ~ s(year, k = 5) +
                           s(lon, lat, k = 5) +
                           s(doy, k = 5) +
                           s(roms_temperature, k = 5) +
                           s(roms_salinity, k = 5),
                         ~ s(year, k = 5) +
                           s(lon, lat, k = 5) +
                           s(doy, k = 5)),
                    data = yfs_larvae,
                    family = ziplss(),
                    method = "REML") #ziplss
summary(yfs_larvae_gam3)
par(mfrow = c(2, 2))
gam.check(yfs_larvae_gam3)

# Tweedie
yfs_larvae_gam4 <- gam(counts ~ s(year, k = 5) +
                      s(lon, lat) +
                      s(doy, k = 5) +
                      s(roms_temperature, k = 5) +
                      s(roms_salinity, k = 5),
                    data = yfs_larvae,
                    family = tw(),
                    method = "REML")
summary(yfs_larvae_gam4)
gam.check(yfs_larvae_gam4)
plot(yfs_larvae_gam4)
# concerned about