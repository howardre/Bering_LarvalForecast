# Libraries
library(marmap)
library(here)
library(dplyr)


# Load fish data
clean_data <- function(file, roms){
  df <- as.data.frame(filter(readRDS(here('data', file))))
  temperature <- as.data.frame(readRDS(here('data', roms)))
  df$mean_temp <- temperature$mean[match(df$year, temperature$year)]
  df$catch <- df$larvalcatchper10m2 + 1
  df <- na.omit(df)
  df$presence <- 1 * (df$count > 0)
  return(df)
}

pk_egg <- clean_data('pk_egg.rds', 'roms_temps.rds')

# Load bathymetry
BS_bathy <- getNOAA.bathy(lon1 = -176.5, lon2 = -156.5, 
                          lat1 = 51, lat2 = 65.5, 
                          resolution = 1)
blues <- c("lightsteelblue4", "lightsteelblue3", "lightsteelblue2", "lightsteelblue1")
greys <- c(grey(0.6), grey(0.93), grey(0.99))



windows(width = 17,
        height = 14,)
par(mfrow = c(1, 1),
    family = 'serif',
    mar = c(4, 5, 3, .2) + .15)
plot.bathy(
  BS_bathy,
  image = T,
  axes = T,
  lwd = 0.03,
  land = T,
  n = 0,
  bpal = list(c(0,
                max(BS_bathy),
                greys),
              c(min(BS_bathy),
                0,
                blues)),
  xlim = c(-176.5,-156.8),
  ylim = c(51, 65.5),
  ylab = "Latitude °N",
  xlab = "Longitude °W",
  main = "",
  cex.lab = 2,
  cex.main = 2.5,
  cex.axis = 1.8)
points(pk_egg$lon[pk_egg$presence == 0],
       pk_egg$lat[pk_egg$presence == 0],
       pch = 4,
       col = 'red',
       cex = 1.1)
points(pk_egg$lon[pk_egg$presence == 1],
       pk_egg$lat[pk_egg$presence == 1],
       pch = 18,
       col = 'black',
       cex = 1)


