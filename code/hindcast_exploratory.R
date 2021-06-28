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
yfs_egg_basenb <- gam(count ~ offset(log(volume_filtered)) + 
                      factor(year) + 
                      s(lon, lat) + 
                      s(doy, k = 5),
    data = yfs_egg,
    family = nb())
summary(yfs_egg_basenb)

yfs_egg_basenb$family$getTheta(TRUE) # 0.1145417

plot(yfs_egg_basenb, select = 2, main = "Neg Binomial")
par(mfrow = c(2, 2))
gam.check(yfs_egg_basenb)

# quasiPoisson
yfs_egg_basep <- gam(count ~ offset(log(volume_filtered)) + 
                       factor(year) + 
                       s(lon, lat) + 
                       s(doy),
                     data = yfs_egg,
                     family = quasipoisson(link = "log"))
summary(yfs_egg_basep)

plot(yfs_egg_basep, select = 2, main = "Poisson")

par(mfrow = c(2, 2))
gam.check(yfs_egg_basep)

# Tweedie
yfs_egg_baset <- gam(larvalcatchper10m2 ~ factor(year) + 
                       s(lon, lat) + 
                       s(doy),
                     data = yfs_egg,
                     family = tw(link = 'log'),
                     method = 'REML')
summary(yfs_egg_baset)

par(mfrow = c(2, 2))
gam.check(yfs_egg_baset)

plot(yfs_egg_baset, select = 2, main = "Tweedie")

# Zero-inflated Poisson (1 stage and 2 stage)
yfs_egg_zip <- gam(count ~ offset(log(volume_filtered)) + 
                       s(year) +  # cannot use year as a factor
                       s(lon, lat) + 
                       s(doy),
                     data = yfs_egg,
                     family = ziP())
summary(yfs_egg_basezip)

par(mfrow = c(2, 2))
gam.check(yfs_egg_zip)

plot(yfs_egg_zip, select = 3, main = "Zero-inflated Poisson (1-stage)")

yfs_egg_baseziplss <- gam(list(count ~ s(year) +
                           s(lon, lat) +
                           s(doy),
                         ~ s(year) +
                           s(lon, lat) +
                           s(doy)),
                    data = yfs_egg,
                    family = ziplss()) #ziplss
summary(yfs_egg_baseziplss)
par(mfrow = c(2, 2))
gam.check(yfs_egg_baseziplss)

# Zero-inflated negative binomial
# Haven't found a good option for two-stage with different models
# Using the gamlss package
yfs_egg_zinb1 <- gamlss(count ~ cs(year)  +
                         cs(doy) +
                         cs(lon, lat),
                       ~ year +
                         doy +
                         lon, lat,
                        data = na.omit(yfs_egg),
                 family = NBI(), k = 5) # won't allow smooth terms for both models
summary(yfs_egg_zinb1)

# Using the pscl package
yfs_egg_zinb2 <- zeroinfl(count ~ year +
                               doy +
                               roms_temperature |
                               year +
                               doy ,
                             data = yfs_egg)
summary(yfs_egg_zinb2)



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
                    family = ziplss()) #ziplss
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

pk_egg_basenb$family$getTheta(TRUE) # 0.2818324

par(mfrow = c(2, 2))
plot(pk_egg_basenb, select = 2, main = "DOY")
plot(pk_egg_basenb, select = 3, main = "Salinity")
plot(pk_egg_basenb, select = 4, main = "Temperature")

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

par(mfrow = c(2, 2))
plot(pk_egg_basep, select = 2, main = "DOY")
plot(pk_egg_basep, select = 3, main = "Salinity")
plot(pk_egg_basep, select = 4, main = "Temperature")

par(mfrow = c(2, 2))
gam.check(pk_egg_basep)

# Tweedie
pk_egg_baset <- gam(larvalcatchper10m2 ~ factor(year) + 
                      s(lon, lat) + 
                      s(doy),
                     data = pk_egg,
                     family = tw(link = 'log'),
                     method = 'REML')
summary(pk_egg_baset)

par(mfrow = c(2, 2))
plot(pk_egg_baset, select = 2, main = "DOY")
plot(pk_egg_baset, select = 3, main = "Salinity")
plot(pk_egg_baset, select = 4, main = "Temperature")

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
par(mfrow = c(2, 2))
gam.check(pk_egg_ziplss)

# Zero-inflated negative binomial
# Haven't found a good option for two-stage with different models
# Using the gamlss package
library(gamlss)
pk_egg_zinb1 <- gamlss(count ~ cs(year)  +
                         cs(doy) +
                         cs(lon, lat) +
                         cs(roms_salinity) +
                         cs(roms_temperature),
                       ~ year +
                         doy +
                         lon, lat,
                        data = na.omit(pk_egg),
                        family = NBI(), k = 5) # won't allow smooth terms for both models
summary(pk_egg_zinb1)

# Using the brms package modified to frequentist
# library(brms) # takes too long
# pk_egg_zinb2 <- brm(count ~ s(year, k = 4) +
#                       s(doy, k = 4) +
#                       s(lon, lat),
#                     data = pk_egg,
#                     family = zero_inflated_negbinomial(),
#                     chains = 4,
#                     cores = 4,
#                     control = list(adapt_delta = 0.999))
# summary(pk_egg_zinb2, WAIC = FALSE)
# plot(marginal_effects(pk_egg_zinb2))

# Using the VGAM package (likely the best option)


# Two-part binomial and gaussian
# binomial
pk_egg$presence <- 1 * (pk_egg$count > 0)
pk_egg_gam1 <- gam(presence ~ factor(year) +
                     s(doy) +
                     s(lon, lat),
                   data = pk_egg,
                   family = "binomial")
summary(pk_egg_gam1)
par(mfrow = c(2, 2))
gam.check(pk_egg_gam1)

# gaussian
pk_egg_gam2 <- gam(log(larvalcatchper10m2 + 1) ~ factor(year) +
                     s(doy) +
                     s(lon, lat) +
                     s(roms_temperature) +
                     s(roms_salinity),
                   data = pk_egg[pk_egg$larvalcatchper10m2 > 0, ])
summary(pk_egg_gam2)
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
pk_egg_model <- pk_egg_baset
grid_extent_eggpk$year <- 2014
grid_extent_eggpk$doy <- median(pk_egg$doy)
grid_extent_eggpk$pred <- predict(pk_egg_model, newdata = grid_extent_eggpk)
grid_extent_eggpk$pred[grid_extent_eggpk$dist > 30000] <- NA

# Plot
windows(width = 8, height = 3.5)
par(mfrow = c(1, 2), 
    mai = c(0.8, 0.9, 0.5, 0.5))
image(lond,
      latd,
      t(matrix(grid_extent_eggpk$pred,
               nrow = length(latd),
               ncol = length(lond),
               byrow = T)),
      col = viridis(100, option = "F", direction = 1),
      ylab = "Latitude",
      xlab = "Longitude",
      xlim = c(-176.5, -156.5),
      ylim = c(52, 62),
      main = 'Distribution',
      cex.main = 1.2,
      cex.lab = 1.1,
      cex.axis = 1.1)
# contour(unique(bathy_dat$lon), # need to change to marmap contours
#         sort(unique(bathy_dat$lat)),
#         bathy_mat,
#         levels = -c(50, 200),
#         labcex = 1,
#         col = 'black',
#         add = T)
symbols(pk_egg$lon[pk_egg$larvalcatchper10m2 > 0],
        pk_egg$lat[pk_egg$larvalcatchper10m2 > 0],
        circles = log(pk_egg$larvalcatchper10m2 + 1)[pk_egg$larvalcatchper10m2 > 0],
        inches = 0.1,
        bg = alpha('grey', 0.1),
        fg = alpha('black', 0.05),
        add = T)
points(pk_egg$lon[pk_egg$larvalcatchper10m2 == 0], 
       pk_egg$lat[pk_egg$larvalcatchper10m2 == 0], 
       pch = '')
map("worldHires",
    fill = T,
    col = "wheat4",
    add = T)
# mtext('Egg density ln(n/10m2)', 1, line = 4.0, cex = 1.2)

# Plot phenology
grid_extent_eggpk2 <- data.frame('lon' = rep(-170, 100),
                           'lat' = rep(57, 100),
                           'doy' = seq(min(pk_egg$doy), max(pk_egg$doy), length = 100),
                           'year' = rep(2014, 100))
grid_extent_eggpk2$pred <- predict(pk_egg_model, newdata = grid_extent_eggpk2)
grid_extent_eggpk2$se <- predict(pk_egg_model, newdata = grid_extent_eggpk2, se = T)[[2]]
grid_extent_eggpk2$pred.up <- grid_extent_eggpk2$pred + 1.96 * grid_extent_eggpk2$se
grid_extent_eggpk2$pred.lw <- grid_extent_eggpk2$pred - 1.96 * grid_extent_eggpk2$se
plot(grid_extent_eggpk2$doy,
     grid_extent_eggpk2$pred,
     main = 'Phenology',
     type = 'l',
     ylab = 'Egg density ln(n/10m2)',
     xlab = 'Day of the year',
     cex.lab = 1.1,
     cex.axis = 1.1,
     cex.main = 1.2,
     xlim = c(60, 215),
     ylim = range(c(grid_extent_eggpk2$pred.up, grid_extent_eggpk2$pred.lw)),
     col = 'blue',
     lwd = 2)
polygon(c(grid_extent_eggpk2$doy, rev(grid_extent_eggpk2$doy)),
        c(grid_extent_eggpk2$pred.lw, rev(grid_extent_eggpk2$pred.up)),
        col = alpha('gray', 0.6),
        lty = 0)
dev.copy(jpeg, here('results', 'pollok_base_egg.jpg'), 
         height = 3.5, width = 8, res = 200, units = 'in')
dev.off()

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

pk_larvae_basenb$family$getTheta(TRUE) # 0.3223544

windows()
par(mfrow = c(2, 2))
plot(pk_larvae_basenb, select = 2, main = "DOY")
plot(pk_larvae_basenb, select = 3, main = "Salinity")
plot(pk_larvae_basenb, select = 4, main = "Temperature")
dev.copy(jpeg, here('results', 'pollok_larvae_basenb.jpg'), 
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

par(mfrow = c(2, 2))
plot(pk_larvae_basep, select = 2, main = "DOY")
plot(pk_larvae_basep, select = 3, main = "Salinity")
plot(pk_larvae_basep, select = 4, main = "Temperature")

par(mfrow = c(2, 2))
gam.check(pk_larvae_basep)

# Tweedie
pk_larvae_baset <- gam(larvalcatchper10m2 ~ factor(year) + 
                      s(lon, lat) + 
                      s(doy),
                    data = pk_larvae,
                    family = tw(link = 'log'),
                    method = 'REML')
summary(pk_larvae_baset)

par(mfrow = c(2, 2))
plot(pk_larvae_baset, select = 2, main = "DOY")
plot(pk_larvae_baset, select = 3, main = "Salinity")
plot(pk_larvae_baset, select = 4, main = "Temperature")

par(mfrow = c(2, 2))
gam.check(pk_larvae_baset)

plot(pk_larvae_baset, select = 2, main = "Tweedie")

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
par(mfrow = c(2, 2))
gam.check(pk_larvae_ziplss)

# Plot best model
# Plot results of average geography and phenology along with decrease of MSE
# Prediction grid
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
pk_larvae_model <- pk_larvae_baset
grid_extent_larvaepk$year <- 2014
grid_extent_larvaepk$doy <- median(pk_larvae$doy)
grid_extent_larvaepk$pred <- predict(pk_larvae_model, newdata = grid_extent_larvaepk)
grid_extent_larvaepk$pred[grid_extent_larvaepk$dist > 30000] <- NA

# Plot
windows(width = 8, height = 3.5)
par(mfrow = c(1, 2), 
    mai = c(0.8, 0.9, 0.5, 0.5))
image(lond,
      latd,
      t(matrix(grid_extent_larvaepk$pred,
               nrow = length(latd),
               ncol = length(lond),
               byrow = T)),
      col = viridis(100, option = "F", direction = 1),
      ylab = "Latitude",
      xlab = "Longitude",
      xlim = c(-176.5, -156.5),
      ylim = c(52, 62),
      main = 'Distribution',
      cex.main = 1.2,
      cex.lab = 1.1,
      cex.axis = 1.1)
# contour(unique(bathy_dat$lon), # need to change to marmap contours
#         sort(unique(bathy_dat$lat)),
#         bathy_mat,
#         levels = -c(50, 200),
#         labcex = 1,
#         col = 'black',
#         add = T)
symbols(pk_larvae$lon[pk_larvae$larvalcatchper10m2 > 0],
        pk_larvae$lat[pk_larvae$larvalcatchper10m2 > 0],
        circles = log(pk_larvae$larvalcatchper10m2 + 1)[pk_larvae$larvalcatchper10m2 > 0],
        inches = 0.1,
        bg = alpha('grey', 0.1),
        fg = alpha('black', 0.05),
        add = T)
points(pk_larvae$lon[pk_larvae$larvalcatchper10m2 == 0], 
       pk_larvae$lat[pk_larvae$larvalcatchper10m2 == 0], 
       pch = '')
map("worldHires",
    fill = T,
    col = "wheat4",
    add = T)
# mtext('larvae density ln(n/10m2)', 1, line = 4.0, cex = 1.2)

# Plot phenology
grid_extent_larvaepk2 <- data.frame('lon' = rep(-170, 100),
                           'lat' = rep(57, 100),
                           'doy' = seq(min(pk_larvae$doy), max(pk_larvae$doy), length = 100),
                           'year' = rep(2014, 100))
grid_extent_larvaepk2$pred <- predict(pk_larvae_model, newdata = grid_extent_larvaepk2)
grid_extent_larvaepk2$se <- predict(pk_larvae_model, newdata = grid_extent_larvaepk2, se = T)[[2]]
grid_extent_larvaepk2$pred.up <- grid_extent_larvaepk2$pred + 1.96 * grid_extent_larvaepk2$se
grid_extent_larvaepk2$pred.lw <- grid_extent_larvaepk2$pred - 1.96 * grid_extent_larvaepk2$se
plot(grid_extent_larvaepk2$doy,
     grid_extent_larvaepk2$pred,
     main = 'Phenology',
     type = 'l',
     ylab = 'Larval density ln(n/10m2)',
     xlab = 'Day of the year',
     cex.lab = 1.1,
     cex.axis = 1.1,
     cex.main = 1.2,
     xlim = c(60, 215),
     ylim = range(c(grid_extent_larvaepk2$pred.up, grid_extent_larvaepk2$pred.lw)),
     col = 'blue',
     lwd = 2)
polygon(c(grid_extent_larvaepk2$doy, rev(grid_extent_larvaepk2$doy)),
        c(grid_extent_larvaepk2$pred.lw, rev(grid_extent_larvaepk2$pred.up)),
        col = alpha('gray', 0.6),
        lty = 0)
dev.copy(jpeg, here('results', 'pollok_base_larvae.jpg'), height = 3.5, width = 8, res = 200, units = 'in')
dev.off()