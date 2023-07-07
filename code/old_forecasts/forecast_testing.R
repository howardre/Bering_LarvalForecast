### Libraries, functions, and data ----
library(maps)
library(maptools)
library(mapdata)
library(spacetime)
library(fields)
library(here)
library(tidyverse)
library(lubridate)
library(date)
library(rgdal)
library(RColorBrewer)
library(mgcv)
library(RANN)
library(scales)
library(magick)
source(here('code/functions', 'distance_function.R'))

# New method to use bias corrected ROMS output
get_preds <- function(data, the_year, doy,
                      the_month, proj, temp_output,
                      salt_output, formula){
  nlat = 40
  nlon = 60
  latd = seq(min(data$lat), max(data$lat), length.out = nlat)
  lond = seq(min(data$lon), max(data$lon), length.out = nlon)
  
  grid_extent <- expand.grid(lond, latd)
  names(grid_extent) <- c('lon', 'lat')
  
  # Calculate distance of each grid point to closest 'positive observation'
  grid_extent$dist <- NA
  
  for (k in 1:nrow(grid_extent)) {
    dist <- distance_function(grid_extent$lat[k],
                              grid_extent$lon[k],
                              data$lat,
                              data$lon)
    grid_extent$dist[k] <- min(dist)
  }
  
  # Assign a within sample year and doy to the grid data
  grid_extent$year <- the_year
  grid_extent$doy <- rep(doy, length(grid_extent))
  grid_extent$month <- the_month
  
  temp_output$lon <- -temp_output$lon
  salt_output$lon <- -salt_output$lon
  
  # Use RANN package to match nearest temperature value based on lat, lon, month, year
  bc_temps <- temp_output %>% filter(month == the_month & year == the_year & projection == proj)
  bc_salts <- salt_output %>% filter(month == the_month & year == the_year & projection == proj)
  
  grid_extent[, c(7, 8)] <- as.data.frame(RANN::nn2(bc_temps[, c('lat', 'lon')],
                                                    grid_extent[, c('lat', 'lon')],
                                                    k = 1))
  grid_extent$roms_temperature <- bc_temps[c(grid_extent$nn.idx), 10] # Match nearest temp
  grid_extent <- grid_extent[-c(7:8)] # remove extra columns before predicting
  grid_extent[, c(8, 9)] <- as.data.frame(RANN::nn2(bc_salts[, c('lat', 'lon')],
                                                    grid_extent[, c('lat', 'lon')],
                                                    k = 1))
  grid_extent$roms_salinity <- bc_salts[c(grid_extent$nn.idx), 10]
  grid_extent <- grid_extent[-c(8, 9)] 
  
  # Calculate mean temperature
  temp_filtered <- temp_output %>% filter(lon >= -170 & lon <= -165, 
                                          lat >= 56 & lat <= 58,
                                          month >= 2 & month <= 4,
                                          year == the_year)
  mean <- mean(temp_filtered$bc, na.rm = T)
  grid_extent$mean_temp <- mean
  
  # Parameterized model
  gam <- formula
  
  grid_extent$pred <- exp(predict(gam,
                                  newdata = grid_extent,
                                  type = "link",
                                  exclude = "s(year)"))
  grid_extent$pred[grid_extent$dist > 30000] <- NA
  return(grid_extent)
 }

pred_loop <- function(range, data, doy, month, 
                      proj, temp_output, salt_output,
                      the_formula){
  grids <- list()
  for(j in range) {
    grid <- get_preds(data, j, doy, month,
                      proj, temp_output, salt_output,
                      the_formula)
    grids[[paste("year", j, sep = "")]] <- grid
  }
  return(grids)
}


grid_predict <- function(grid, title){
  nlat = 40
  nlon = 60
  latd = seq(min(grid$lat), max(grid$lat), length.out = nlat)
  lond = seq(min(grid$lon), max(grid$lon), length.out = nlon)
  my_color = colorRampPalette(rev(c("#FFFFCC", "#FBF2A8", "#F9E585",
                                    "#F5D363", "#EFBA55", "#EAA352",
                                    "#E68C51", "#E0754F", "#D75C4D",
                                    "#BB4A48", "#994240", "#763931", 
                                    "#542D20", "#352311", "#191900")))
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
        col = my_color(100), 
        ylab = "Latitude",
        xlab = "Longitude",
        xlim = c(-176.5, -156.5),
        ylim = c(52, 62),
        zlim = c(min(grid$pred, na.rm = T), 
                 max(grid$pred, na.rm = T)),
        main = title,
        cex.main = 1.2,
        cex.lab = 1.1,
        cex.axis = 1.1)
  maps::map("worldHires",
            fill = T,
            col = "wheat4",
            add = T)
  image.plot(legend.only = T,
             col = my_color(100),
             legend.shrink = 0.2,
             smallplot = c(.79, .82, .20, .37),
             legend.cex = 0.8,
             axis.args = list(cex.axis = 0.8),
             legend.width = 0.5,
             legend.mar = 6,
             zlim = c(min(grid$pred, na.rm = T), 
                      max(grid$pred, na.rm = T)),
             legend.args = list("Scaled Abundance",
                                side = 2, cex = 1))
}



# Functions
# Loads the data, adds ROMS temperatures, sets up catch for GAM
load_data <- function(file, data, temps){
  data <- as.data.frame(readRDS(here('data', file)))
  data$mean_temp <- temps$mean[match(data$year, temps$year)]
  data$catch <- data$larvalcatchper10m2 + 1
  return(data)
}

# Create GAM formulas
formula_pheno <- function(data){
  gam(catch ~ s(year, bs = 're', k = 9) +
        s(doy, k = 9) +
        s(lon, lat) +
        s(roms_temperature, k = 9) +
        s(roms_salinity, k = 9) +
        s(doy, by = mean_temp, k = 9), # phenology
      data = data,
      family = tw(link = "log"),
      method = 'REML')
}

formula_geog <- function(data){
  gam(catch ~ s(year, bs = 're', k = 9) +
        s(doy, k = 9) +
        s(lon, lat) +
        s(roms_temperature, k = 9) +
        s(roms_salinity, k = 9) +
        s(lat, lon, by = mean_temp, k = 9), # geography
      data = data,
      family = tw(link = "log"),
      method = 'REML')
}

base_dir <- getwd()

# Load ROMS temperature means and forecast
roms_temps <- readRDS(here('data', 'roms_temps.rds'))

### Pollock Eggs --------------------------------------------------------------------------------------------------------------------------
pk_egg <- load_data('pk_egg.rds', pk_egg, roms_temps)
pkegg_formula <- formula_geog(pk_egg)

#### Forecast and average into 3 time periods ---------------------------------------------------------------------------------------------
##### CESM 126 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
cesm_temps1 <- readRDS(here('data', 'cesm_forecast_temp1.rds'))
cesm_salts1 <- readRDS(here('data', 'cesm_forecast_salt1.rds'))

test1 <- pred_loop(2015:2039, pk_egg, 134,
                   5, 'ssp126', cesm_temps1, 
                   cesm_salts1, pkegg_formula)

# Plot each individual year
# Multiple years show strange patterns
grid_predict(test1$year2015,
             "Forecasted Distribution 2015 \n CESM SSP1-2.6")
grid_predict(test1$year2016,
             "Forecasted Distribution 2016 \n CESM SSP1-2.6")
grid_predict(test1$year2017,
             "Forecasted Distribution 2017 \n CESM SSP1-2.6")
grid_predict(test1$year2018,
             "Forecasted Distribution 2018 \n CESM SSP1-2.6")
grid_predict(test1$year2019,
             "Forecasted Distribution 2019 \n CESM SSP1-2.6")
grid_predict(test1$year2020,
             "Forecasted Distribution 2020 \n CESM SSP1-2.6")
grid_predict(test1$year2021,
             "Forecasted Distribution 2021 \n CESM SSP1-2.6")
grid_predict(test1$year2022,
             "Forecasted Distribution 2022 \n CESM SSP1-2.6")
grid_predict(test1$year2023,
             "Forecasted Distribution 2023 \n CESM SSP1-2.6")
grid_predict(test1$year2024,
             "Forecasted Distribution 2024 \n CESM SSP1-2.6")
grid_predict(test1$year2025,
             "Forecasted Distribution 2025 \n CESM SSP1-2.6")
grid_predict(test1$year2026,
             "Forecasted Distribution 2026 \n CESM SSP1-2.6")
grid_predict(test1$year2027,
             "Forecasted Distribution 2027 \n CESM SSP1-2.6")
grid_predict(test1$year2028,
             "Forecasted Distribution 2028 \n CESM SSP1-2.6")
grid_predict(test1$year2029,
             "Forecasted Distribution 2029 \n CESM SSP1-2.6")
grid_predict(test1$year2030,
             "Forecasted Distribution 2030 \n CESM SSP1-2.6")
grid_predict(test1$year2031,
             "Forecasted Distribution 2031 \n CESM SSP1-2.6")
grid_predict(test1$year2032,
             "Forecasted Distribution 2032 \n CESM SSP1-2.6")
grid_predict(test1$year2033,
             "Forecasted Distribution 2033 \n CESM SSP1-2.6")
grid_predict(test1$year2034,
             "Forecasted Distribution 2034 \n CESM SSP1-2.6")
grid_predict(test1$year2035,
             "Forecasted Distribution 2035 \n CESM SSP1-2.6")
grid_predict(test1$year2036,
             "Forecasted Distribution 2036 \n CESM SSP1-2.6")
grid_predict(test1$year2037,
             "Forecasted Distribution 2037 \n CESM SSP1-2.6")
grid_predict(test1$year2038,
             "Forecasted Distribution 2038 \n CESM SSP1-2.6")
grid_predict(test1$year2039,
             "Forecasted Distribution 2039 \n CESM SSP1-2.6")





# Test distance function
nlat = 60
nlon = 80
latd = seq(min(pk_egg$lat), max(pk_egg$lat), length.out = nlat)
lond = seq(min(pk_egg$lon), max(pk_egg$lon), length.out = nlon)

grid_extent <- expand.grid(lond, latd)
names(grid_extent) <- c('lon', 'lat')

# Calculate distance of each grid point to closest 'positive observation'
grid_extent$dist <- NA

for (k in 1:nrow(grid_extent)) {
  dist <- distance_function(grid_extent$lat[k],
                            grid_extent$lon[k],
                            pk_egg$lat,
                            pk_egg$lon)
  grid_extent$dist[k] <- min(dist)
}

# Assign a within sample year and 134 to the grid pk_egg
grid_extent$year <- 2021
grid_extent$doy <- rep(134, length(grid_extent))
grid_extent$month <- 5

cesm_temps1 <- cesm_temps1 %>% 
  mutate(lat = lat,
         lon = case_when(lon >= 180 ~ lon - 360,
                         lon < 180 ~ lon))
cesm_salts1 <- cesm_salts1 %>% 
  mutate(lat = lat,
         lon = case_when(lon >= 180 ~ lon - 360,
                         lon < 180 ~ lon))

# Use RANN package to match nearest temperature value based on lat, lon, month, year
bc_temps <- cesm_temps1 %>% filter(month == 5 & year == 2021 & projection == 'ssp126')
bc_salts <- cesm_salts1 %>% filter(month == 5 & year == 2021 & projection == 'ssp126')

grid_extent[, c(7, 8)] <- as.data.frame(RANN::nn2(bc_temps[, c('lon', 'lat')],
                                                  grid_extent[, c('lon', 'lat')],
                                                  k = 1))
grid_extent$roms_temperature <- bc_temps[c(grid_extent$nn.idx), 10] # Match nearest temp
grid_extent <- grid_extent[-c(7:8)] # remove extra columns before predicting
grid_extent[, c(8, 9)] <- as.data.frame(RANN::nn2(bc_salts[, c('lon', 'lat')],
                                                  grid_extent[, c('lon', 'lat')],
                                                  k = 1))
grid_extent$roms_salinity <- bc_salts[c(grid_extent$nn.idx), 10]
grid_extent <- grid_extent[-c(8, 9)] 

# Try making coarser resolution grid
#put points into a 'SpatialPointsDataFrame' 'sp' object
# Rotate the grid
bc_temps <- cesm_temps1 %>% filter(month == 5 & year == 2021 & projection == 'ssp126')

bc_temps <- bc_temps %>% 
  mutate(lat = lat,
         lon = case_when(lon >= 180 ~ lon - 360,
                         lon < 180 ~ lon))

coords <- cbind(x = bc_temps[['lon']], y = bc_temps[['lat']])
sPDF <- SpatialPointsDataFrame(coords, data = bc_temps)

#set number of rows & columns in the grid
nrows <- 150
ncols <- 150

#setting extents from the data
xmn <- min(bc_temps[['lon']]) 
ymn <- min(bc_temps[['lat']]) 
xmx <- max(bc_temps[['lon']]) 
ymx <- max(bc_temps[['lat']]) 

#create a grid
blankRaster <- raster(nrows = nrows, ncols = ncols, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx)
#adding data into raster to avoid 'no data' error
blankRaster[] <- 1:ncell(blankRaster)

#calc mean (or other function) of points per cell 
rasterMeanPoints <- rasterize(x = sPDF, y = blankRaster, field = 'bc', fun = mean)

coords <- as.matrix(coordinates(rasterMeanPoints))
temp_output <- data.frame(lon = coords[, 1], lat = coords[, 2],
                          as.data.frame(rasterMeanPoints))

#plot to get an idea whether it's doing the right thing
plot(rasterMeanPoints)

grid_extent[, c(7, 8)] <- as.data.frame(RANN::nn2(temp_output[, c('lon', 'lat')],
                                                  grid_extent[, c('lon', 'lat')],
                                                  k = 1))
grid_extent$roms_temperature <- temp_output[c(grid_extent$nn.idx), 3] # Match nearest temp
grid_extent <- grid_extent[-c(7, 8)] # remove extra columns before predicting
grid_extent[, c(8, 9)] <- as.data.frame(RANN::nn2(bc_salts[, c('lon', 'lat')],
                                                  grid_extent[, c('lon', 'lat')],
                                                  k = 1))
grid_extent$roms_salinity <- bc_salts[c(grid_extent$nn.idx), 10]
grid_extent <- grid_extent[-c(8, 9)] 




# Calculate mean temperature
temp_filtered <- cesm_temps1 %>% filter(lon >= -170 & lon <= -165, 
                                        lat >= 56 & lat <= 58,
                                        month >= 2 & month <= 4,
                                        year == 2021)
mean <- mean(temp_filtered$bc, na.rm = T)
grid_extent$mean_temp <- mean

# Parameterized model
gam <- pkegg_formula

grid_extent$pred <- exp(predict(gam,
                                newdata = grid_extent,
                                type = "link",
                                exclude = "s(year)"))
grid_extent$pred[grid_extent$dist > 30000] <- NA


### Plot temperatures
nlat = 60
nlon = 80
latd = seq(min(grid_extent$lat), max(grid_extent$lat), length.out = nlat)
lond = seq(min(grid_extent$lon), max(grid_extent$lon), length.out = nlon)
my_color = colorRampPalette(rev(c("#FFFFCC", "#FBF2A8", "#F9E585",
                                  "#F5D363", "#EFBA55", "#EAA352",
                                  "#E68C51", "#E0754F", "#D75C4D",
                                  "#BB4A48", "#994240", "#763931", 
                                  "#542D20", "#352311", "#191900")))
image(lond,
      latd,
      t(matrix(grid_extent$roms_temperature,
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
      t(matrix(grid_extent$roms_temperature,
               nrow = length(latd),
               ncol = length(lond),
               byrow = T)),
      col = my_color(100), 
      ylab = "Latitude",
      xlab = "Longitude",
      xlim = c(-176.5, -156.5),
      ylim = c(52, 62),
      zlim = c(min(grid_extent$roms_temperature, na.rm = T), 
               max(grid_extent$roms_temperature, na.rm = T)),
      main = "ROMS Temperature",
      cex.main = 1.2,
      cex.lab = 1.1,
      cex.axis = 1.1)
maps::map("worldHires",
          fill = T,
          col = "wheat4",
          add = T)
image.plot(legend.only = T,
           col = my_color(100),
           legend.shrink = 0.2,
           smallplot = c(.79, .82, .20, .37),
           legend.cex = 0.8,
           axis.args = list(cex.axis = 0.8),
           legend.width = 0.5,
           legend.mar = 6,
           zlim = c(min(grid_extent$roms_temperature, na.rm = T), 
                    max(grid_extent$roms_temperature, na.rm = T)),
           legend.args = list("Temperature",
                              side = 2, cex = 1))


terms <- predict(pkegg_formula,
                 newdata = 
                   type = "terms") # largest by far is DOY