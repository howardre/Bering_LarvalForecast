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
source(here('code/functions', 'distance_function.R'))


# Load ROMS temperature means and forecast
roms_temps <- readRDS(here('data', 'roms_temps.rds'))

# Load fish data
yfs_egg <- as.data.frame(filter((readRDS(here('data', 'yfs_egg.rds'))),
                                lat >= 52 & lat <= 62,
                                lon >= -176.5 & lon <= -156.5))
yfs_egg$mean_temp <- roms_temps$mean[match(yfs_egg$year, roms_temps$year)]
yfs_egg$catch <- yfs_egg$larvalcatchper10m2 + 1

yfs_larvae <- as.data.frame(filter(readRDS(here('data', 'yfs_larvae.rds')),
                                   lat >= 52 & lat <= 62,
                                   lon >= -176.5 & lon <= -156.5))
yfs_larvae$mean_temp <- roms_temps$mean[match(yfs_larvae$year, roms_temps$year)]
yfs_larvae$catch <- yfs_larvae$larvalcatchper10m2 + 1

# Match ROMS output function
varid_match <- function(data, model_output1, model_output2, list){
  data$roms_date <- NA
  data$roms_temperature <- NA
  data$roms_salinity <- NA
  for (i in 1:nrow(data)) {
    idx_time <- order(abs(model_output1[[list]][[3]] - data$date[i]))[1]
    data$roms_date[i] <- model_output1[[list]][[3]][idx_time]
    idx_grid <- order(distance_function(data$lat[i],
                                        data$lon[i],
                                        c(model_output1[[list]][[2]]),
                                        c(model_output1[[list]][[1]])))[1]
    data$roms_temperature[i] <- c(model_output1[[list]][[4]][, , idx_time])[idx_grid]
    data$roms_salinity[i] <- c(model_output2[[list]][[4]][, , idx_time])[idx_grid]
  }
  return(data)
}

# Function to make maps
grid_predict <- function(grid, title){
  nlat = 80
  nlon = 120
  latd = seq(min(grid$lat), max(grid$lat), length.out = nlat)
  lond = seq(min(grid$lon), max(grid$lon), length.out = nlon)
  my_color = colorRampPalette(rev(c("#FFFFCC", "#FBF2A8", "#F9E585",
                                    "#F5D363", "#EFBA55", "#EAA352",
                                    "#E68C51", "#E0754F", "#D75C4D",
                                    "#BB4A48", "#994240", "#763931", 
                                    "#542D20", "#352311", "#191900")))
  image(lond,
        latd,
        t(matrix(grid$avg_pred,
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
        t(matrix(grid$avg_pred,
                 nrow = length(latd),
                 ncol = length(lond),
                 byrow = T)),
        col = my_color(100), 
        ylab = "Latitude",
        xlab = "Longitude",
        xlim = c(-176.5, -156.5),
        ylim = c(52, 62),
        zlim = c(min(grid$avg_pred, na.rm = T), 
                 max(grid$avg_pred, na.rm = T)),
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
             zlim = c(min(grid$avg_pred, na.rm = T), 
                      max(grid$avg_pred, na.rm = T)),
             legend.args = list("Avg. Predicted \n Occurrence",
                                side = 2, cex = 1))
}

egg_formula <- gam(catch ~ s(year) +
                     s(doy, k = 8) +
                     s(lon, lat) +
                     s(roms_temperature, k = 6) +
                     s(roms_salinity, k = 6),
                   data = yfs_egg,
                   family = tw(link = 'log'),
                   method = 'REML')

larval_formula <- gam(catch ~ s(year) +
                        s(doy, k = 8) +
                        s(lon, lat) +
                        s(roms_temperature, k = 6) +
                        s(roms_salinity, k = 6) +
                        s(lat, lon, by = mean_temp, k = 6),
                      data = yfs_larvae,
                      family = tw(link = 'log'),
                      method = 'REML')

### Bias correction testing ----
# Average by month for hindcast and historical found in each fish dataset
# Need to match the forecast to the fish dataset by lat, lon, year, month
# Then subtract the baseline (historical) from the forecast to get the deltas
# Finally add the deltas to the hindcast to get final values



# Function to get predictions
get_preds <- function(data, year, date, doy, 
                      start_date, end_date,
                      temp_output, salt_output,
                      list, formula){
  # Prediction grid
  nlat = 80
  nlon = 120
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
  grid_extent$year <- year
  grid_extent$date <- rep(as.Date(date),
                          length(grid_extent))
  grid_extent$doy <- rep(doy, length(grid_extent))
  
  # Attach ROMs forecast
  grid_extent <- varid_match(grid_extent, temp_output, salt_output, list)
  
  # Calculate mean temperature
  time_index <- temp_output[[list]][[3]] >= start_date & temp_output[[list]][[3]] <= end_date
  temp_array <- temp_output[[list]][[4]][, , time_index]
  temp_data <- as.data.frame(cbind(lon = as.vector(temp_output[[list]][[1]]), 
                                   lat = as.vector(temp_output[[list]][[2]]), 
                                   temp = as.vector(temp_array)))
  temp_filtered <- temp_data %>% filter(lon >= -170 & lon <= -165, lat >= 56 & lat <= 58)
  mean <- mean(temp_filtered$temp, na.rm = T)
  
  grid_extent$mean_temp <- mean
  
  # Parameterized model
  gam <- formula
  
  # Predict on forecasted output
  grid_extent$pred <- predict(gam,
                              newdata = grid_extent,
                              type = "response")
  grid_extent$pred[grid_extent$dist > 30000] <- NA
  
  return(grid_extent)
  
}

# Function to loop through years
pred_loop <- function(range, data, doy, 
                      temp_output, salt_output,
                      list, formula){
  grids <- list()
  for(j in range) {
    date1 <- paste(j, "-05-10", sep = "")
    date2 <- paste(j, "-02-01", sep = "")
    date3 <- paste(j, "-04-30", sep = "")
    grid <- get_preds(data, j, date1, doy,
                      date2, date3,
                      temp_output, salt_output,
                      list, formula)
    grids[[paste("year", j, sep = "")]] <- grid
  }
  return(grids)
}


### Yellowfin Eggs --------------------------------------------------------------------------------------------------------------------------
#### Forecast and average into 3 time periods ---------------------------------------------------------------------------------------------
##### CESM 126 ----------------------------------------------------------------------------------------------------------------------------
temps_cesm_ssp126 <- readRDS(here('data', 'temps_cesm_ssp126.rds'))
salts_cesm_ssp126 <- readRDS(here('data', 'salts_cesm_ssp126.rds'))

## 2015 - 2039
grids_yfsegg1 <- pred_loop(2015:2019, yfs_egg, 130, 
                           temps_cesm_ssp126,
                           salts_cesm_ssp126, 1,
                           egg_formula)
grids_yfsegg2 <- pred_loop(2020:2024, yfs_egg, 130, 
                           temps_cesm_ssp126,
                           salts_cesm_ssp126, 2,
                           egg_formula)
grids_yfsegg3 <- pred_loop(2025:2029, yfs_egg, 130,
                           temps_cesm_ssp126,
                           salts_cesm_ssp126, 3,
                           egg_formula)
grids_yfsegg4 <- pred_loop(2030:2034, yfs_egg, 130, 
                           temps_cesm_ssp126,
                           salts_cesm_ssp126, 4,
                           egg_formula)
grids_yfsegg5 <- pred_loop(2035:2039, yfs_egg, 130, 
                           temps_cesm_ssp126,
                           salts_cesm_ssp126, 5,
                           egg_formula)

# Combine into one data frame
df_yfsegg1 <- list(grids_yfsegg1[[1]], grids_yfsegg1[[2]], grids_yfsegg1[[3]], 
                   grids_yfsegg1[[4]], grids_yfsegg1[[5]], grids_yfsegg2[[1]], 
                   grids_yfsegg2[[2]], grids_yfsegg2[[3]], grids_yfsegg2[[4]],
                   grids_yfsegg2[[5]], grids_yfsegg3[[1]], grids_yfsegg3[[2]], 
                   grids_yfsegg3[[3]], grids_yfsegg3[[4]], grids_yfsegg3[[5]],
                   grids_yfsegg4[[1]], grids_yfsegg4[[2]], grids_yfsegg4[[3]], 
                   grids_yfsegg4[[4]], grids_yfsegg4[[5]], grids_yfsegg5[[1]], 
                   grids_yfsegg5[[2]], grids_yfsegg5[[3]], grids_yfsegg5[[4]],
                   grids_yfsegg5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_yfsegg1), fixed = T)
df_yfsegg_avg1_cesm126 <- data.frame(lat = df_yfsegg1$lat, 
                                     lon = df_yfsegg1$lon, 
                                     dist = df_yfsegg1$dist,
                                     avg_pred = rowSums(df_yfsegg1[, x])/25)
saveRDS(df_yfsegg_avg1_cesm126, file = here("data", "df_yfsegg_avg1_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfsegg_avg1_cesm126, "Forecasted Distribution 2015 - 2039 \n CESM SSP126")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_egg_cesm_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_yfsegg6 <- pred_loop(2040:2044, yfs_egg, 130, 
                           temps_cesm_ssp126,
                           salts_cesm_ssp126, 6,
                           egg_formula)
grids_yfsegg7 <- pred_loop(2045:2049, yfs_egg, 130, 
                           temps_cesm_ssp126,
                           salts_cesm_ssp126, 7,
                           egg_formula)
grids_yfsegg8 <- pred_loop(2050:2054, yfs_egg, 130, 
                           temps_cesm_ssp126,
                           salts_cesm_ssp126, 8,
                           egg_formula)
grids_yfsegg9 <- pred_loop(2055:2059, yfs_egg, 130, 
                           temps_cesm_ssp126,
                           salts_cesm_ssp126, 9,
                           egg_formula)
grids_yfsegg10 <- pred_loop(2060:2064, yfs_egg, 130, 
                            temps_cesm_ssp126,
                            salts_cesm_ssp126, 10,
                            egg_formula)
grids_yfsegg11 <- pred_loop(2065:2069, yfs_egg, 130, 
                            temps_cesm_ssp126,
                            salts_cesm_ssp126, 11,
                            egg_formula)

# Combine into one data frame
df_yfsegg2 <- list(grids_yfsegg6[[1]], grids_yfsegg6[[2]], grids_yfsegg6[[3]], 
                   grids_yfsegg6[[4]], grids_yfsegg6[[5]], grids_yfsegg7[[1]], 
                   grids_yfsegg7[[2]], grids_yfsegg7[[3]], grids_yfsegg7[[4]],
                   grids_yfsegg7[[5]], grids_yfsegg8[[1]], grids_yfsegg8[[2]], 
                   grids_yfsegg8[[3]], grids_yfsegg8[[4]], grids_yfsegg8[[5]],
                   grids_yfsegg9[[1]], grids_yfsegg9[[2]], grids_yfsegg9[[3]], 
                   grids_yfsegg9[[4]], grids_yfsegg9[[5]], grids_yfsegg10[[1]], 
                   grids_yfsegg10[[2]], grids_yfsegg10[[3]], grids_yfsegg10[[4]],
                   grids_yfsegg11[[5]], grids_yfsegg11[[1]], grids_yfsegg11[[2]],
                   grids_yfsegg11[[3]], grids_yfsegg11[[4]], grids_yfsegg11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_yfsegg2), fixed = T)
df_yfsegg_avg2_cesm126 <- data.frame(lat = df_yfsegg2$lat, 
                                     lon = df_yfsegg2$lon, 
                                     dist = df_yfsegg2$dist,
                                     avg_pred = rowSums(df_yfsegg2[, x])/30)
saveRDS(df_yfsegg_avg2_cesm126, file = here("data", "df_yfsegg_avg2_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfsegg_avg2_cesm126, "Forecasted Distribution 2040 - 2069 \n CESM SSP126")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_egg_cesm_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_yfsegg12 <- pred_loop(2070:2074, yfs_egg, 130,
                            temps_cesm_ssp126,
                            salts_cesm_ssp126, 12,
                            egg_formula)
grids_yfsegg13 <- pred_loop(2075:2079, yfs_egg, 130,
                            temps_cesm_ssp126,
                            salts_cesm_ssp126, 13,
                            egg_formula)
grids_yfsegg14 <- pred_loop(2080:2084, yfs_egg, 130,
                            temps_cesm_ssp126,
                            salts_cesm_ssp126, 14,
                            egg_formula)
grids_yfsegg15 <- pred_loop(2085:2089, yfs_egg, 130,
                            temps_cesm_ssp126,
                            salts_cesm_ssp126, 15,
                            egg_formula)
grids_yfsegg16 <- pred_loop(2090:2094, yfs_egg, 130,
                            temps_cesm_ssp126,
                            salts_cesm_ssp126, 16,
                            egg_formula)
grids_yfsegg17 <- pred_loop(2095:2099, yfs_egg, 130, 
                            temps_cesm_ssp126,
                            salts_cesm_ssp126, 17,
                            egg_formula)

# Combine into one data frame
df_yfsegg3 <- list(grids_yfsegg12[[1]], grids_yfsegg12[[2]], grids_yfsegg12[[3]], 
                   grids_yfsegg12[[4]], grids_yfsegg12[[5]], grids_yfsegg13[[1]], 
                   grids_yfsegg13[[2]], grids_yfsegg13[[3]], grids_yfsegg13[[4]],
                   grids_yfsegg13[[5]], grids_yfsegg14[[1]], grids_yfsegg14[[2]], 
                   grids_yfsegg14[[3]], grids_yfsegg14[[4]], grids_yfsegg14[[5]],
                   grids_yfsegg15[[1]], grids_yfsegg15[[2]], grids_yfsegg15[[3]], 
                   grids_yfsegg15[[4]], grids_yfsegg15[[5]], grids_yfsegg16[[1]], 
                   grids_yfsegg16[[2]], grids_yfsegg16[[3]], grids_yfsegg16[[4]],
                   grids_yfsegg17[[5]], grids_yfsegg17[[1]], grids_yfsegg17[[2]],
                   grids_yfsegg17[[3]], grids_yfsegg17[[4]], grids_yfsegg17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_yfsegg3), fixed = T)
df_yfsegg_avg3_cesm126 <- data.frame(lat = df_yfsegg3$lat, 
                                     lon = df_yfsegg3$lon, 
                                     dist = df_yfsegg3$dist,
                                     avg_pred = rowSums(df_yfsegg3[, x])/30)
saveRDS(df_yfsegg_avg3_cesm126, file = here("data", "df_yfsegg_avg3_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfsegg_avg3_cesm126, "Forecasted Distribution 2070 - 2099 \n CESM SSP126")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_egg_cesm_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

##### CESM 585 -----------------------------------------------------------------------------------------------------------------
temps_cesm_ssp585 <- readRDS(here('data', 'temps_cesm_ssp585.rds'))
salts_cesm_ssp585 <- readRDS(here('data', 'salts_cesm_ssp585.rds'))

## 2015 - 2039
grids_yfsegg1 <- pred_loop(2015:2019, yfs_egg, 130, 
                           temps_cesm_ssp585,
                           salts_cesm_ssp585, 1,
                           egg_formula)
grids_yfsegg2 <- pred_loop(2020:2024, yfs_egg, 130, 
                           temps_cesm_ssp585,
                           salts_cesm_ssp585, 2,
                           egg_formula)
grids_yfsegg3 <- pred_loop(2025:2029, yfs_egg, 130,
                           temps_cesm_ssp585,
                           salts_cesm_ssp585, 3,
                           egg_formula)
grids_yfsegg4 <- pred_loop(2030:2034, yfs_egg, 130, 
                           temps_cesm_ssp585,
                           salts_cesm_ssp585, 4,
                           egg_formula)
grids_yfsegg5 <- pred_loop(2035:2039, yfs_egg, 130, 
                           temps_cesm_ssp585,
                           salts_cesm_ssp585, 5,
                           egg_formula)

# Combine into one data frame
df_yfsegg4 <- list(grids_yfsegg1[[1]], grids_yfsegg1[[2]], grids_yfsegg1[[3]], 
                   grids_yfsegg1[[4]], grids_yfsegg1[[5]], grids_yfsegg2[[1]], 
                   grids_yfsegg2[[2]], grids_yfsegg2[[3]], grids_yfsegg2[[4]],
                   grids_yfsegg2[[5]], grids_yfsegg3[[1]], grids_yfsegg3[[2]], 
                   grids_yfsegg3[[3]], grids_yfsegg3[[4]], grids_yfsegg3[[5]],
                   grids_yfsegg4[[1]], grids_yfsegg4[[2]], grids_yfsegg4[[3]], 
                   grids_yfsegg4[[4]], grids_yfsegg4[[5]], grids_yfsegg5[[1]], 
                   grids_yfsegg5[[2]], grids_yfsegg5[[3]], grids_yfsegg5[[4]],
                   grids_yfsegg5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_yfsegg4), fixed = T)
df_yfsegg_avg4_cesm585 <- data.frame(lat = df_yfsegg4$lat, 
                                     lon = df_yfsegg4$lon, 
                                     dist = df_yfsegg4$dist,
                                     avg_pred = rowSums(df_yfsegg4[, x])/25)
saveRDS(df_yfsegg_avg4_cesm585, file = here("data", "df_yfsegg_avg4_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfsegg_avg4_cesm585, "Forecasted Distribution 2015 - 2039 \n CESM SSP585")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_egg_cesm_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_yfsegg6 <- pred_loop(2040:2044, yfs_egg, 130, 
                           temps_cesm_ssp585,
                           salts_cesm_ssp585, 6,
                           egg_formula)
grids_yfsegg7 <- pred_loop(2045:2049, yfs_egg, 130, 
                           temps_cesm_ssp585,
                           salts_cesm_ssp585, 7,
                           egg_formula)
grids_yfsegg8 <- pred_loop(2050:2054, yfs_egg, 130, 
                           temps_cesm_ssp585,
                           salts_cesm_ssp585, 8,
                           egg_formula)
grids_yfsegg9 <- pred_loop(2055:2059, yfs_egg, 130, 
                           temps_cesm_ssp585,
                           salts_cesm_ssp585, 9,
                           egg_formula)
grids_yfsegg10 <- pred_loop(2060:2064, yfs_egg, 130, 
                            temps_cesm_ssp585,
                            salts_cesm_ssp585, 10,
                            egg_formula)
grids_yfsegg11 <- pred_loop(2065:2069, yfs_egg, 130, 
                            temps_cesm_ssp585,
                            salts_cesm_ssp585, 11,
                            egg_formula)

# Combine into one data frame
df_yfsegg5 <- list(grids_yfsegg6[[1]], grids_yfsegg6[[2]], grids_yfsegg6[[3]], 
                   grids_yfsegg6[[4]], grids_yfsegg6[[5]], grids_yfsegg7[[1]], 
                   grids_yfsegg7[[2]], grids_yfsegg7[[3]], grids_yfsegg7[[4]],
                   grids_yfsegg7[[5]], grids_yfsegg8[[1]], grids_yfsegg8[[2]], 
                   grids_yfsegg8[[3]], grids_yfsegg8[[4]], grids_yfsegg8[[5]],
                   grids_yfsegg9[[1]], grids_yfsegg9[[2]], grids_yfsegg9[[3]], 
                   grids_yfsegg9[[4]], grids_yfsegg9[[5]], grids_yfsegg10[[1]], 
                   grids_yfsegg10[[2]], grids_yfsegg10[[3]], grids_yfsegg10[[4]],
                   grids_yfsegg11[[5]], grids_yfsegg11[[1]], grids_yfsegg11[[2]],
                   grids_yfsegg11[[3]], grids_yfsegg11[[4]], grids_yfsegg11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_yfsegg5), fixed = T)
df_yfsegg_avg5_cesm585 <- data.frame(lat = df_yfsegg5$lat, 
                                     lon = df_yfsegg5$lon, 
                                     dist = df_yfsegg5$dist,
                                     avg_pred = rowSums(df_yfsegg5[, x])/30)
saveRDS(df_yfsegg_avg5_cesm585, file = here("data", "df_yfsegg_avg5_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfsegg_avg5_cesm585, "Forecasted Distribution 2040 - 2069 \n CESM SSP585")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_egg_cesm_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_yfsegg12 <- pred_loop(2070:2074, yfs_egg, 130,
                            temps_cesm_ssp585,
                            salts_cesm_ssp585, 12,
                            egg_formula)
grids_yfsegg13 <- pred_loop(2075:2079, yfs_egg, 130,
                            temps_cesm_ssp585,
                            salts_cesm_ssp585, 13,
                            egg_formula)
grids_yfsegg14 <- pred_loop(2080:2084, yfs_egg, 130,
                            temps_cesm_ssp585,
                            salts_cesm_ssp585, 14,
                            egg_formula)
grids_yfsegg15 <- pred_loop(2085:2089, yfs_egg, 130,
                            temps_cesm_ssp585,
                            salts_cesm_ssp585, 15,
                            egg_formula)
grids_yfsegg16 <- pred_loop(2090:2094, yfs_egg, 130,
                            temps_cesm_ssp585,
                            salts_cesm_ssp585, 16,
                            egg_formula)
grids_yfsegg17 <- pred_loop(2095:2099, yfs_egg, 130, 
                            temps_cesm_ssp585,
                            salts_cesm_ssp585, 17,
                            egg_formula)

# Combine into one data frame
df_yfsegg6 <- list(grids_yfsegg12[[1]], grids_yfsegg12[[2]], grids_yfsegg12[[3]], 
                   grids_yfsegg12[[4]], grids_yfsegg12[[5]], grids_yfsegg13[[1]], 
                   grids_yfsegg13[[2]], grids_yfsegg13[[3]], grids_yfsegg13[[4]],
                   grids_yfsegg13[[5]], grids_yfsegg14[[1]], grids_yfsegg14[[2]], 
                   grids_yfsegg14[[3]], grids_yfsegg14[[4]], grids_yfsegg14[[5]],
                   grids_yfsegg15[[1]], grids_yfsegg15[[2]], grids_yfsegg15[[3]], 
                   grids_yfsegg15[[4]], grids_yfsegg15[[5]], grids_yfsegg16[[1]], 
                   grids_yfsegg16[[2]], grids_yfsegg16[[3]], grids_yfsegg16[[4]],
                   grids_yfsegg17[[5]], grids_yfsegg17[[1]], grids_yfsegg17[[2]],
                   grids_yfsegg17[[3]], grids_yfsegg17[[4]], grids_yfsegg17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_yfsegg6), fixed = T)
df_yfsegg_avg6_cesm585 <- data.frame(lat = df_yfsegg6$lat, 
                                     lon = df_yfsegg6$lon, 
                                     dist = df_yfsegg6$dist,
                                     avg_pred = rowSums(df_yfsegg6[, x])/30)
saveRDS(df_yfsegg_avg6_cesm585, file = here("data", "df_yfsegg_avg6_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfsegg_avg6_cesm585, "Forecasted Distribution 2070 - 2099 \n CESM SSP585")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_egg_cesm_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### GFDL 126 ----------------------------------------------------------------------------------------------------------------------------
temps_gfdl_ssp126 <- readRDS(here('data', 'temps_gfdl_ssp126.rds'))
salts_gfdl_ssp126 <- readRDS(here('data', 'salts_gfdl_ssp126.rds'))

## 2015 - 2039
grids_yfsegg1 <- pred_loop(2015:2019, yfs_egg, 130, 
                           temps_gfdl_ssp126,
                           salts_gfdl_ssp126, 1,
                           egg_formula)
grids_yfsegg2 <- pred_loop(2020:2024, yfs_egg, 130, 
                           temps_gfdl_ssp126,
                           salts_gfdl_ssp126, 2,
                           egg_formula)
grids_yfsegg3 <- pred_loop(2025:2029, yfs_egg, 130,
                           temps_gfdl_ssp126,
                           salts_gfdl_ssp126, 3,
                           egg_formula)
grids_yfsegg4 <- pred_loop(2030:2034, yfs_egg, 130, 
                           temps_gfdl_ssp126,
                           salts_gfdl_ssp126, 4,
                           egg_formula)
grids_yfsegg5 <- pred_loop(2035:2039, yfs_egg, 130, 
                           temps_gfdl_ssp126,
                           salts_gfdl_ssp126, 5,
                           egg_formula)

# Combine into one data frame
df_yfsegg1 <- list(grids_yfsegg1[[1]], grids_yfsegg1[[2]], grids_yfsegg1[[3]], 
                   grids_yfsegg1[[4]], grids_yfsegg1[[5]], grids_yfsegg2[[1]], 
                   grids_yfsegg2[[2]], grids_yfsegg2[[3]], grids_yfsegg2[[4]],
                   grids_yfsegg2[[5]], grids_yfsegg3[[1]], grids_yfsegg3[[2]], 
                   grids_yfsegg3[[3]], grids_yfsegg3[[4]], grids_yfsegg3[[5]],
                   grids_yfsegg4[[1]], grids_yfsegg4[[2]], grids_yfsegg4[[3]], 
                   grids_yfsegg4[[4]], grids_yfsegg4[[5]], grids_yfsegg5[[1]], 
                   grids_yfsegg5[[2]], grids_yfsegg5[[3]], grids_yfsegg5[[4]],
                   grids_yfsegg5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_yfsegg1), fixed = T)
df_yfsegg_avg1_gfdl126 <- data.frame(lat = df_yfsegg1$lat, 
                                     lon = df_yfsegg1$lon, 
                                     dist = df_yfsegg1$dist,
                                     avg_pred = rowSums(df_yfsegg1[, x])/25)
saveRDS(df_yfsegg_avg1_gfdl126, file = here("data", "df_yfsegg_avg1_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfsegg_avg1_gfdl126, "Forecasted Distribution 2015 - 2039 \n GFDL SSP126")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_egg_gfdl_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_yfsegg6 <- pred_loop(2040:2044, yfs_egg, 130, 
                           temps_gfdl_ssp126,
                           salts_gfdl_ssp126, 6,
                           egg_formula)
grids_yfsegg7 <- pred_loop(2045:2049, yfs_egg, 130, 
                           temps_gfdl_ssp126,
                           salts_gfdl_ssp126, 7,
                           egg_formula)
grids_yfsegg8 <- pred_loop(2050:2054, yfs_egg, 130, 
                           temps_gfdl_ssp126,
                           salts_gfdl_ssp126, 8,
                           egg_formula)
grids_yfsegg9 <- pred_loop(2055:2059, yfs_egg, 130, 
                           temps_gfdl_ssp126,
                           salts_gfdl_ssp126, 9,
                           egg_formula)
grids_yfsegg10 <- pred_loop(2060:2064, yfs_egg, 130, 
                            temps_gfdl_ssp126,
                            salts_gfdl_ssp126, 10,
                            egg_formula)
grids_yfsegg11 <- pred_loop(2065:2069, yfs_egg, 130, 
                            temps_gfdl_ssp126,
                            salts_gfdl_ssp126, 11,
                            egg_formula)

# Combine into one data frame
df_yfsegg2 <- list(grids_yfsegg6[[1]], grids_yfsegg6[[2]], grids_yfsegg6[[3]], 
                   grids_yfsegg6[[4]], grids_yfsegg6[[5]], grids_yfsegg7[[1]], 
                   grids_yfsegg7[[2]], grids_yfsegg7[[3]], grids_yfsegg7[[4]],
                   grids_yfsegg7[[5]], grids_yfsegg8[[1]], grids_yfsegg8[[2]], 
                   grids_yfsegg8[[3]], grids_yfsegg8[[4]], grids_yfsegg8[[5]],
                   grids_yfsegg9[[1]], grids_yfsegg9[[2]], grids_yfsegg9[[3]], 
                   grids_yfsegg9[[4]], grids_yfsegg9[[5]], grids_yfsegg10[[1]], 
                   grids_yfsegg10[[2]], grids_yfsegg10[[3]], grids_yfsegg10[[4]],
                   grids_yfsegg11[[5]], grids_yfsegg11[[1]], grids_yfsegg11[[2]],
                   grids_yfsegg11[[3]], grids_yfsegg11[[4]], grids_yfsegg11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_yfsegg2), fixed = T)
df_yfsegg_avg2_gfdl126 <- data.frame(lat = df_yfsegg2$lat, 
                                     lon = df_yfsegg2$lon, 
                                     dist = df_yfsegg2$dist,
                                     avg_pred = rowSums(df_yfsegg2[, x])/30)
saveRDS(df_yfsegg_avg2_gfdl126, file = here("data", "df_yfsegg_avg2_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfsegg_avg2_gfdl126, "Forecasted Distribution 2040 - 2069 \n GFDL SSP126")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_egg_gfdl_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_yfsegg12 <- pred_loop(2070:2074, yfs_egg, 130,
                            temps_gfdl_ssp126,
                            salts_gfdl_ssp126, 12,
                            egg_formula)
grids_yfsegg13 <- pred_loop(2075:2079, yfs_egg, 130,
                            temps_gfdl_ssp126,
                            salts_gfdl_ssp126, 13,
                            egg_formula)
grids_yfsegg14 <- pred_loop(2080:2084, yfs_egg, 130,
                            temps_gfdl_ssp126,
                            salts_gfdl_ssp126, 14,
                            egg_formula)
grids_yfsegg15 <- pred_loop(2085:2089, yfs_egg, 130,
                            temps_gfdl_ssp126,
                            salts_gfdl_ssp126, 15,
                            egg_formula)
grids_yfsegg16 <- pred_loop(2090:2094, yfs_egg, 130,
                            temps_gfdl_ssp126,
                            salts_gfdl_ssp126, 16,
                            egg_formula)
grids_yfsegg17 <- pred_loop(2095:2099, yfs_egg, 130, 
                            temps_gfdl_ssp126,
                            salts_gfdl_ssp126, 17,
                            egg_formula)

# Combine into one data frame
df_yfsegg3 <- list(grids_yfsegg12[[1]], grids_yfsegg12[[2]], grids_yfsegg12[[3]], 
                   grids_yfsegg12[[4]], grids_yfsegg12[[5]], grids_yfsegg13[[1]], 
                   grids_yfsegg13[[2]], grids_yfsegg13[[3]], grids_yfsegg13[[4]],
                   grids_yfsegg13[[5]], grids_yfsegg14[[1]], grids_yfsegg14[[2]], 
                   grids_yfsegg14[[3]], grids_yfsegg14[[4]], grids_yfsegg14[[5]],
                   grids_yfsegg15[[1]], grids_yfsegg15[[2]], grids_yfsegg15[[3]], 
                   grids_yfsegg15[[4]], grids_yfsegg15[[5]], grids_yfsegg16[[1]], 
                   grids_yfsegg16[[2]], grids_yfsegg16[[3]], grids_yfsegg16[[4]],
                   grids_yfsegg17[[5]], grids_yfsegg17[[1]], grids_yfsegg17[[2]],
                   grids_yfsegg17[[3]], grids_yfsegg17[[4]], grids_yfsegg17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_yfsegg3), fixed = T)
df_yfsegg_avg3_gfdl126 <- data.frame(lat = df_yfsegg3$lat, 
                                     lon = df_yfsegg3$lon, 
                                     dist = df_yfsegg3$dist,
                                     avg_pred = rowSums(df_yfsegg3[, x])/30)
saveRDS(df_yfsegg_avg3_gfdl126, file = here("data", "df_yfsegg_avg3_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfsegg_avg3_gfdl126, "Forecasted Distribution 2070 - 2099 \n GFDL SSP126")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_egg_gfdl_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

##### GFDL 585 -----------------------------------------------------------------------------------------------------------------
temps_gfdl_ssp585 <- readRDS(here('data', 'temps_gfdl_ssp585.rds'))
salts_gfdl_ssp585 <- readRDS(here('data', 'salts_gfdl_ssp585.rds'))

## 2015 - 2039
grids_yfsegg1 <- pred_loop(2015:2019, yfs_egg, 130, 
                           temps_gfdl_ssp585,
                           salts_gfdl_ssp585, 1,
                           egg_formula)
grids_yfsegg2 <- pred_loop(2020:2024, yfs_egg, 130, 
                           temps_gfdl_ssp585,
                           salts_gfdl_ssp585, 2,
                           egg_formula)
grids_yfsegg3 <- pred_loop(2025:2029, yfs_egg, 130,
                           temps_gfdl_ssp585,
                           salts_gfdl_ssp585, 3,
                           egg_formula)
grids_yfsegg4 <- pred_loop(2030:2034, yfs_egg, 130, 
                           temps_gfdl_ssp585,
                           salts_gfdl_ssp585, 4,
                           egg_formula)
grids_yfsegg5 <- pred_loop(2035:2039, yfs_egg, 130, 
                           temps_gfdl_ssp585,
                           salts_gfdl_ssp585, 5,
                           egg_formula)

# Combine into one data frame
df_yfsegg4 <- list(grids_yfsegg1[[1]], grids_yfsegg1[[2]], grids_yfsegg1[[3]], 
                   grids_yfsegg1[[4]], grids_yfsegg1[[5]], grids_yfsegg2[[1]], 
                   grids_yfsegg2[[2]], grids_yfsegg2[[3]], grids_yfsegg2[[4]],
                   grids_yfsegg2[[5]], grids_yfsegg3[[1]], grids_yfsegg3[[2]], 
                   grids_yfsegg3[[3]], grids_yfsegg3[[4]], grids_yfsegg3[[5]],
                   grids_yfsegg4[[1]], grids_yfsegg4[[2]], grids_yfsegg4[[3]], 
                   grids_yfsegg4[[4]], grids_yfsegg4[[5]], grids_yfsegg5[[1]], 
                   grids_yfsegg5[[2]], grids_yfsegg5[[3]], grids_yfsegg5[[4]],
                   grids_yfsegg5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_yfsegg4), fixed = T)
df_yfsegg_avg4_gfdl585 <- data.frame(lat = df_yfsegg4$lat, 
                                     lon = df_yfsegg4$lon, 
                                     dist = df_yfsegg4$dist,
                                     avg_pred = rowSums(df_yfsegg4[, x])/25)
saveRDS(df_yfsegg_avg4_gfdl585, file = here("data", "df_yfsegg_avg4_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfsegg_avg4_gfdl585, "Forecasted Distribution 2015 - 2039 \n GFDL SSP585")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_egg_gfdl_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_yfsegg6 <- pred_loop(2040:2044, yfs_egg, 130, 
                           temps_gfdl_ssp585,
                           salts_gfdl_ssp585, 6,
                           egg_formula)
grids_yfsegg7 <- pred_loop(2045:2049, yfs_egg, 130, 
                           temps_gfdl_ssp585,
                           salts_gfdl_ssp585, 7,
                           egg_formula)
grids_yfsegg8 <- pred_loop(2050:2054, yfs_egg, 130, 
                           temps_gfdl_ssp585,
                           salts_gfdl_ssp585, 8,
                           egg_formula)
grids_yfsegg9 <- pred_loop(2055:2059, yfs_egg, 130, 
                           temps_gfdl_ssp585,
                           salts_gfdl_ssp585, 9,
                           egg_formula)
grids_yfsegg10 <- pred_loop(2060:2064, yfs_egg, 130, 
                            temps_gfdl_ssp585,
                            salts_gfdl_ssp585, 10,
                            egg_formula)
grids_yfsegg11 <- pred_loop(2065:2069, yfs_egg, 130, 
                            temps_gfdl_ssp585,
                            salts_gfdl_ssp585, 11,
                            egg_formula)

# Combine into one data frame
df_yfsegg5 <- list(grids_yfsegg6[[1]], grids_yfsegg6[[2]], grids_yfsegg6[[3]], 
                   grids_yfsegg6[[4]], grids_yfsegg6[[5]], grids_yfsegg7[[1]], 
                   grids_yfsegg7[[2]], grids_yfsegg7[[3]], grids_yfsegg7[[4]],
                   grids_yfsegg7[[5]], grids_yfsegg8[[1]], grids_yfsegg8[[2]], 
                   grids_yfsegg8[[3]], grids_yfsegg8[[4]], grids_yfsegg8[[5]],
                   grids_yfsegg9[[1]], grids_yfsegg9[[2]], grids_yfsegg9[[3]], 
                   grids_yfsegg9[[4]], grids_yfsegg9[[5]], grids_yfsegg10[[1]], 
                   grids_yfsegg10[[2]], grids_yfsegg10[[3]], grids_yfsegg10[[4]],
                   grids_yfsegg11[[5]], grids_yfsegg11[[1]], grids_yfsegg11[[2]],
                   grids_yfsegg11[[3]], grids_yfsegg11[[4]], grids_yfsegg11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_yfsegg5), fixed = T)
df_yfsegg_avg5_gfdl585 <- data.frame(lat = df_yfsegg5$lat, 
                                     lon = df_yfsegg5$lon, 
                                     dist = df_yfsegg5$dist,
                                     avg_pred = rowSums(df_yfsegg5[, x])/30)
saveRDS(df_yfsegg_avg5_gfdl585, file = here("data", "df_yfsegg_avg5_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfsegg_avg5_gfdl585, "Forecasted Distribution 2040 - 2069 \n GFDL SSP585")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_egg_gfdl_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_yfsegg12 <- pred_loop(2070:2074, yfs_egg, 130,
                            temps_gfdl_ssp585,
                            salts_gfdl_ssp585, 12,
                            egg_formula)
grids_yfsegg13 <- pred_loop(2075:2079, yfs_egg, 130,
                            temps_gfdl_ssp585,
                            salts_gfdl_ssp585, 13,
                            egg_formula)
grids_yfsegg14 <- pred_loop(2080:2084, yfs_egg, 130,
                            temps_gfdl_ssp585,
                            salts_gfdl_ssp585, 14,
                            egg_formula)
grids_yfsegg15 <- pred_loop(2085:2089, yfs_egg, 130,
                            temps_gfdl_ssp585,
                            salts_gfdl_ssp585, 15,
                            egg_formula)
grids_yfsegg16 <- pred_loop(2090:2094, yfs_egg, 130,
                            temps_gfdl_ssp585,
                            salts_gfdl_ssp585, 16,
                            egg_formula)
grids_yfsegg17 <- pred_loop(2095:2099, yfs_egg, 130, 
                            temps_gfdl_ssp585,
                            salts_gfdl_ssp585, 17,
                            egg_formula)

# Combine into one data frame
df_yfsegg6 <- list(grids_yfsegg12[[1]], grids_yfsegg12[[2]], grids_yfsegg12[[3]], 
                   grids_yfsegg12[[4]], grids_yfsegg12[[5]], grids_yfsegg13[[1]], 
                   grids_yfsegg13[[2]], grids_yfsegg13[[3]], grids_yfsegg13[[4]],
                   grids_yfsegg13[[5]], grids_yfsegg14[[1]], grids_yfsegg14[[2]], 
                   grids_yfsegg14[[3]], grids_yfsegg14[[4]], grids_yfsegg14[[5]],
                   grids_yfsegg15[[1]], grids_yfsegg15[[2]], grids_yfsegg15[[3]], 
                   grids_yfsegg15[[4]], grids_yfsegg15[[5]], grids_yfsegg16[[1]], 
                   grids_yfsegg16[[2]], grids_yfsegg16[[3]], grids_yfsegg16[[4]],
                   grids_yfsegg17[[5]], grids_yfsegg17[[1]], grids_yfsegg17[[2]],
                   grids_yfsegg17[[3]], grids_yfsegg17[[4]], grids_yfsegg17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_yfsegg6), fixed = T)
df_yfsegg_avg6_gfdl585 <- data.frame(lat = df_yfsegg6$lat, 
                                     lon = df_yfsegg6$lon, 
                                     dist = df_yfsegg6$dist,
                                     avg_pred = rowSums(df_yfsegg6[, x])/30)
saveRDS(df_yfsegg_avg6_gfdl585, file = here("data", "df_yfsegg_avg6_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfsegg_avg6_gfdl585, "Forecasted Distribution 2070 - 2099 \n GFDL SSP585")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_egg_gfdl_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### MIROC 126 ----------------------------------------------------------------------------------------------------------------------------
temps_miroc_ssp126 <- readRDS(here('data', 'temps_miroc_ssp126.rds'))
salts_miroc_ssp126 <- readRDS(here('data', 'salts_miroc_ssp126.rds'))

## 2015 - 2039
grids_yfsegg1 <- pred_loop(2015:2019, yfs_egg, 130, 
                           temps_miroc_ssp126,
                           salts_miroc_ssp126, 1,
                           egg_formula)
grids_yfsegg2 <- pred_loop(2020:2024, yfs_egg, 130, 
                           temps_miroc_ssp126,
                           salts_miroc_ssp126, 2,
                           egg_formula)
grids_yfsegg3 <- pred_loop(2025:2029, yfs_egg, 130,
                           temps_miroc_ssp126,
                           salts_miroc_ssp126, 3,
                           egg_formula)
grids_yfsegg4 <- pred_loop(2030:2034, yfs_egg, 130, 
                           temps_miroc_ssp126,
                           salts_miroc_ssp126, 4,
                           egg_formula)
grids_yfsegg5 <- pred_loop(2035:2039, yfs_egg, 130, 
                           temps_miroc_ssp126,
                           salts_miroc_ssp126, 5,
                           egg_formula)

# Combine into one data frame
df_yfsegg1 <- list(grids_yfsegg1[[1]], grids_yfsegg1[[2]], grids_yfsegg1[[3]], 
                   grids_yfsegg1[[4]], grids_yfsegg1[[5]], grids_yfsegg2[[1]], 
                   grids_yfsegg2[[2]], grids_yfsegg2[[3]], grids_yfsegg2[[4]],
                   grids_yfsegg2[[5]], grids_yfsegg3[[1]], grids_yfsegg3[[2]], 
                   grids_yfsegg3[[3]], grids_yfsegg3[[4]], grids_yfsegg3[[5]],
                   grids_yfsegg4[[1]], grids_yfsegg4[[2]], grids_yfsegg4[[3]], 
                   grids_yfsegg4[[4]], grids_yfsegg4[[5]], grids_yfsegg5[[1]], 
                   grids_yfsegg5[[2]], grids_yfsegg5[[3]], grids_yfsegg5[[4]],
                   grids_yfsegg5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_yfsegg1), fixed = T)
df_yfsegg_avg1_miroc126 <- data.frame(lat = df_yfsegg1$lat, 
                                      lon = df_yfsegg1$lon, 
                                      dist = df_yfsegg1$dist,
                                      avg_pred = rowSums(df_yfsegg1[, x])/25)
saveRDS(df_yfsegg_avg1_miroc126, file = here("data", "df_yfsegg_avg1_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfsegg_avg1_miroc126, "Forecasted Distribution 2015 - 2039 \n MIROC SSP126")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_egg_miroc_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_yfsegg6 <- pred_loop(2040:2044, yfs_egg, 130, 
                           temps_miroc_ssp126,
                           salts_miroc_ssp126, 6,
                           egg_formula)
grids_yfsegg7 <- pred_loop(2045:2049, yfs_egg, 130, 
                           temps_miroc_ssp126,
                           salts_miroc_ssp126, 7,
                           egg_formula)
grids_yfsegg8 <- pred_loop(2050:2054, yfs_egg, 130, 
                           temps_miroc_ssp126,
                           salts_miroc_ssp126, 8,
                           egg_formula)
grids_yfsegg9 <- pred_loop(2055:2059, yfs_egg, 130, 
                           temps_miroc_ssp126,
                           salts_miroc_ssp126, 9,
                           egg_formula)
grids_yfsegg10 <- pred_loop(2060:2064, yfs_egg, 130, 
                            temps_miroc_ssp126,
                            salts_miroc_ssp126, 10,
                            egg_formula)
grids_yfsegg11 <- pred_loop(2065:2069, yfs_egg, 130, 
                            temps_miroc_ssp126,
                            salts_miroc_ssp126, 11,
                            egg_formula)

# Combine into one data frame
df_yfsegg2 <- list(grids_yfsegg6[[1]], grids_yfsegg6[[2]], grids_yfsegg6[[3]], 
                   grids_yfsegg6[[4]], grids_yfsegg6[[5]], grids_yfsegg7[[1]], 
                   grids_yfsegg7[[2]], grids_yfsegg7[[3]], grids_yfsegg7[[4]],
                   grids_yfsegg7[[5]], grids_yfsegg8[[1]], grids_yfsegg8[[2]], 
                   grids_yfsegg8[[3]], grids_yfsegg8[[4]], grids_yfsegg8[[5]],
                   grids_yfsegg9[[1]], grids_yfsegg9[[2]], grids_yfsegg9[[3]], 
                   grids_yfsegg9[[4]], grids_yfsegg9[[5]], grids_yfsegg10[[1]], 
                   grids_yfsegg10[[2]], grids_yfsegg10[[3]], grids_yfsegg10[[4]],
                   grids_yfsegg11[[5]], grids_yfsegg11[[1]], grids_yfsegg11[[2]],
                   grids_yfsegg11[[3]], grids_yfsegg11[[4]], grids_yfsegg11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_yfsegg2), fixed = T)
df_yfsegg_avg2_miroc126 <- data.frame(lat = df_yfsegg2$lat, 
                                      lon = df_yfsegg2$lon, 
                                      dist = df_yfsegg2$dist,
                                      avg_pred = rowSums(df_yfsegg2[, x])/30)
saveRDS(df_yfsegg_avg2_miroc126, file = here("data", "df_yfsegg_avg2_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfsegg_avg2_miroc126, "Forecasted Distribution 2040 - 2069 \n MIROC SSP126")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_egg_miroc_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_yfsegg12 <- pred_loop(2070:2074, yfs_egg, 130,
                            temps_miroc_ssp126,
                            salts_miroc_ssp126, 12,
                            egg_formula)
grids_yfsegg13 <- pred_loop(2075:2079, yfs_egg, 130,
                            temps_miroc_ssp126,
                            salts_miroc_ssp126, 13,
                            egg_formula)
grids_yfsegg14 <- pred_loop(2080:2084, yfs_egg, 130,
                            temps_miroc_ssp126,
                            salts_miroc_ssp126, 14,
                            egg_formula)
grids_yfsegg15 <- pred_loop(2085:2089, yfs_egg, 130,
                            temps_miroc_ssp126,
                            salts_miroc_ssp126, 15,
                            egg_formula)
grids_yfsegg16 <- pred_loop(2090:2094, yfs_egg, 130,
                            temps_miroc_ssp126,
                            salts_miroc_ssp126, 16,
                            egg_formula)
grids_yfsegg17 <- pred_loop(2095:2099, yfs_egg, 130, 
                            temps_miroc_ssp126,
                            salts_miroc_ssp126, 17,
                            egg_formula)

# Combine into one data frame
df_yfsegg3 <- list(grids_yfsegg12[[1]], grids_yfsegg12[[2]], grids_yfsegg12[[3]], 
                   grids_yfsegg12[[4]], grids_yfsegg12[[5]], grids_yfsegg13[[1]], 
                   grids_yfsegg13[[2]], grids_yfsegg13[[3]], grids_yfsegg13[[4]],
                   grids_yfsegg13[[5]], grids_yfsegg14[[1]], grids_yfsegg14[[2]], 
                   grids_yfsegg14[[3]], grids_yfsegg14[[4]], grids_yfsegg14[[5]],
                   grids_yfsegg15[[1]], grids_yfsegg15[[2]], grids_yfsegg15[[3]], 
                   grids_yfsegg15[[4]], grids_yfsegg15[[5]], grids_yfsegg16[[1]], 
                   grids_yfsegg16[[2]], grids_yfsegg16[[3]], grids_yfsegg16[[4]],
                   grids_yfsegg17[[5]], grids_yfsegg17[[1]], grids_yfsegg17[[2]],
                   grids_yfsegg17[[3]], grids_yfsegg17[[4]], grids_yfsegg17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_yfsegg3), fixed = T)
df_yfsegg_avg3_miroc126 <- data.frame(lat = df_yfsegg3$lat, 
                                      lon = df_yfsegg3$lon, 
                                      dist = df_yfsegg3$dist,
                                      avg_pred = rowSums(df_yfsegg3[, x])/30)
saveRDS(df_yfsegg_avg3_miroc126, file = here("data", "df_yfsegg_avg3_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfsegg_avg3_miroc126, "Forecasted Distribution 2070 - 2099 \n MIROC SSP126")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_egg_miroc_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

##### MIROC 585 -----------------------------------------------------------------------------------------------------------------
temps_miroc_ssp585 <- readRDS(here('data', 'temps_miroc_ssp585.rds'))
salts_miroc_ssp585 <- readRDS(here('data', 'salts_miroc_ssp585.rds'))

## 2015 - 2039
grids_yfsegg1 <- pred_loop(2015:2019, yfs_egg, 130, 
                           temps_miroc_ssp585,
                           salts_miroc_ssp585, 1,
                           egg_formula)
grids_yfsegg2 <- pred_loop(2020:2024, yfs_egg, 130, 
                           temps_miroc_ssp585,
                           salts_miroc_ssp585, 2,
                           egg_formula)
grids_yfsegg3 <- pred_loop(2025:2029, yfs_egg, 130,
                           temps_miroc_ssp585,
                           salts_miroc_ssp585, 3,
                           egg_formula)
grids_yfsegg4 <- pred_loop(2030:2034, yfs_egg, 130, 
                           temps_miroc_ssp585,
                           salts_miroc_ssp585, 4,
                           egg_formula)
grids_yfsegg5 <- pred_loop(2035:2039, yfs_egg, 130, 
                           temps_miroc_ssp585,
                           salts_miroc_ssp585, 5,
                           egg_formula)

# Combine into one data frame
df_yfsegg4 <- list(grids_yfsegg1[[1]], grids_yfsegg1[[2]], grids_yfsegg1[[3]], 
                   grids_yfsegg1[[4]], grids_yfsegg1[[5]], grids_yfsegg2[[1]], 
                   grids_yfsegg2[[2]], grids_yfsegg2[[3]], grids_yfsegg2[[4]],
                   grids_yfsegg2[[5]], grids_yfsegg3[[1]], grids_yfsegg3[[2]], 
                   grids_yfsegg3[[3]], grids_yfsegg3[[4]], grids_yfsegg3[[5]],
                   grids_yfsegg4[[1]], grids_yfsegg4[[2]], grids_yfsegg4[[3]], 
                   grids_yfsegg4[[4]], grids_yfsegg4[[5]], grids_yfsegg5[[1]], 
                   grids_yfsegg5[[2]], grids_yfsegg5[[3]], grids_yfsegg5[[4]],
                   grids_yfsegg5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_yfsegg4), fixed = T)
df_yfsegg_avg4_miroc585 <- data.frame(lat = df_yfsegg4$lat, 
                                      lon = df_yfsegg4$lon, 
                                      dist = df_yfsegg4$dist,
                                      avg_pred = rowSums(df_yfsegg4[, x])/25)
saveRDS(df_yfsegg_avg4_miroc585, file = here("data", "df_yfsegg_avg4_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfsegg_avg4_miroc585, "Forecasted Distribution 2015 - 2039 \n MIROC SSP585")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_egg_miroc_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_yfsegg6 <- pred_loop(2040:2044, yfs_egg, 130, 
                           temps_miroc_ssp585,
                           salts_miroc_ssp585, 6,
                           egg_formula)
grids_yfsegg7 <- pred_loop(2045:2049, yfs_egg, 130, 
                           temps_miroc_ssp585,
                           salts_miroc_ssp585, 7,
                           egg_formula)
grids_yfsegg8 <- pred_loop(2050:2054, yfs_egg, 130, 
                           temps_miroc_ssp585,
                           salts_miroc_ssp585, 8,
                           egg_formula)
grids_yfsegg9 <- pred_loop(2055:2059, yfs_egg, 130, 
                           temps_miroc_ssp585,
                           salts_miroc_ssp585, 9,
                           egg_formula)
grids_yfsegg10 <- pred_loop(2060:2064, yfs_egg, 130, 
                            temps_miroc_ssp585,
                            salts_miroc_ssp585, 10,
                            egg_formula)
grids_yfsegg11 <- pred_loop(2065:2069, yfs_egg, 130, 
                            temps_miroc_ssp585,
                            salts_miroc_ssp585, 11,
                            egg_formula)

# Combine into one data frame
df_yfsegg5 <- list(grids_yfsegg6[[1]], grids_yfsegg6[[2]], grids_yfsegg6[[3]], 
                   grids_yfsegg6[[4]], grids_yfsegg6[[5]], grids_yfsegg7[[1]], 
                   grids_yfsegg7[[2]], grids_yfsegg7[[3]], grids_yfsegg7[[4]],
                   grids_yfsegg7[[5]], grids_yfsegg8[[1]], grids_yfsegg8[[2]], 
                   grids_yfsegg8[[3]], grids_yfsegg8[[4]], grids_yfsegg8[[5]],
                   grids_yfsegg9[[1]], grids_yfsegg9[[2]], grids_yfsegg9[[3]], 
                   grids_yfsegg9[[4]], grids_yfsegg9[[5]], grids_yfsegg10[[1]], 
                   grids_yfsegg10[[2]], grids_yfsegg10[[3]], grids_yfsegg10[[4]],
                   grids_yfsegg11[[5]], grids_yfsegg11[[1]], grids_yfsegg11[[2]],
                   grids_yfsegg11[[3]], grids_yfsegg11[[4]], grids_yfsegg11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_yfsegg5), fixed = T)
df_yfsegg_avg5_miroc585 <- data.frame(lat = df_yfsegg5$lat, 
                                      lon = df_yfsegg5$lon, 
                                      dist = df_yfsegg5$dist,
                                      avg_pred = rowSums(df_yfsegg5[, x])/30)
saveRDS(df_yfsegg_avg5_miroc585, file = here("data", "df_yfsegg_avg5_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfsegg_avg5_miroc585, "Forecasted Distribution 2040 - 2069 \n MIROC SSP585")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_egg_miroc_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_yfsegg12 <- pred_loop(2070:2074, yfs_egg, 130,
                            temps_miroc_ssp585,
                            salts_miroc_ssp585, 12,
                            egg_formula)
grids_yfsegg13 <- pred_loop(2075:2079, yfs_egg, 130,
                            temps_miroc_ssp585,
                            salts_miroc_ssp585, 13,
                            egg_formula)
grids_yfsegg14 <- pred_loop(2080:2084, yfs_egg, 130,
                            temps_miroc_ssp585,
                            salts_miroc_ssp585, 14,
                            egg_formula)
grids_yfsegg15 <- pred_loop(2085:2089, yfs_egg, 130,
                            temps_miroc_ssp585,
                            salts_miroc_ssp585, 15,
                            egg_formula)
grids_yfsegg16 <- pred_loop(2090:2094, yfs_egg, 130,
                            temps_miroc_ssp585,
                            salts_miroc_ssp585, 16,
                            egg_formula)
grids_yfsegg17 <- pred_loop(2095:2099, yfs_egg, 130, 
                            temps_miroc_ssp585,
                            salts_miroc_ssp585, 17,
                            egg_formula)

# Combine into one data frame
df_yfsegg6 <- list(grids_yfsegg12[[1]], grids_yfsegg12[[2]], grids_yfsegg12[[3]], 
                   grids_yfsegg12[[4]], grids_yfsegg12[[5]], grids_yfsegg13[[1]], 
                   grids_yfsegg13[[2]], grids_yfsegg13[[3]], grids_yfsegg13[[4]],
                   grids_yfsegg13[[5]], grids_yfsegg14[[1]], grids_yfsegg14[[2]], 
                   grids_yfsegg14[[3]], grids_yfsegg14[[4]], grids_yfsegg14[[5]],
                   grids_yfsegg15[[1]], grids_yfsegg15[[2]], grids_yfsegg15[[3]], 
                   grids_yfsegg15[[4]], grids_yfsegg15[[5]], grids_yfsegg16[[1]], 
                   grids_yfsegg16[[2]], grids_yfsegg16[[3]], grids_yfsegg16[[4]],
                   grids_yfsegg17[[5]], grids_yfsegg17[[1]], grids_yfsegg17[[2]],
                   grids_yfsegg17[[3]], grids_yfsegg17[[4]], grids_yfsegg17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_yfsegg6), fixed = T)
df_yfsegg_avg6_miroc585 <- data.frame(lat = df_yfsegg6$lat, 
                                      lon = df_yfsegg6$lon, 
                                      dist = df_yfsegg6$dist,
                                      avg_pred = rowSums(df_yfsegg6[, x])/30)
saveRDS(df_yfsegg_avg6_miroc585, file = here("data", "df_yfsegg_avg6_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfsegg_avg6_miroc585, "Forecasted Distribution 2070 - 2099 \n MIROC SSP585")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_egg_miroc_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()



### Yellowfin Larvae --------------------------------------------------------------------------------------------------------------------------
#### Forecast and average into 3 time periods ---------------------------------------------------------------------------------------------
##### CESM 126 ----------------------------------------------------------------------------------------------------------------------------
temps_cesm_ssp126 <- readRDS(here('data', 'temps_cesm_ssp126.rds'))
salts_cesm_ssp126 <- readRDS(here('data', 'salts_cesm_ssp126.rds'))

## 2015 - 2039
grids_yfslarvae1 <- pred_loop(2015:2019, yfs_larvae, 130, 
                              temps_cesm_ssp126,
                              salts_cesm_ssp126, 1,
                              larval_formula)
grids_yfslarvae2 <- pred_loop(2020:2024, yfs_larvae, 130, 
                              temps_cesm_ssp126,
                              salts_cesm_ssp126, 2,
                              larval_formula)
grids_yfslarvae3 <- pred_loop(2025:2029, yfs_larvae, 130,
                              temps_cesm_ssp126,
                              salts_cesm_ssp126, 3,
                              larval_formula)
grids_yfslarvae4 <- pred_loop(2030:2034, yfs_larvae, 130, 
                              temps_cesm_ssp126,
                              salts_cesm_ssp126, 4,
                              larval_formula)
grids_yfslarvae5 <- pred_loop(2035:2039, yfs_larvae, 130, 
                              temps_cesm_ssp126,
                              salts_cesm_ssp126, 5,
                              larval_formula)

# Combine into one data frame
df_yfslarvae1 <- list(grids_yfslarvae1[[1]], grids_yfslarvae1[[2]], grids_yfslarvae1[[3]], 
                      grids_yfslarvae1[[4]], grids_yfslarvae1[[5]], grids_yfslarvae2[[1]], 
                      grids_yfslarvae2[[2]], grids_yfslarvae2[[3]], grids_yfslarvae2[[4]],
                      grids_yfslarvae2[[5]], grids_yfslarvae3[[1]], grids_yfslarvae3[[2]], 
                      grids_yfslarvae3[[3]], grids_yfslarvae3[[4]], grids_yfslarvae3[[5]],
                      grids_yfslarvae4[[1]], grids_yfslarvae4[[2]], grids_yfslarvae4[[3]], 
                      grids_yfslarvae4[[4]], grids_yfslarvae4[[5]], grids_yfslarvae5[[1]], 
                      grids_yfslarvae5[[2]], grids_yfslarvae5[[3]], grids_yfslarvae5[[4]],
                      grids_yfslarvae5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_yfslarvae1), fixed = T)
df_yfslarvae_avg1_cesm126 <- data.frame(lat = df_yfslarvae1$lat, 
                                        lon = df_yfslarvae1$lon, 
                                        dist = df_yfslarvae1$dist,
                                        avg_pred = rowSums(df_yfslarvae1[, x])/25)
saveRDS(df_yfslarvae_avg1_cesm126, file = here("data", "df_yfslarvae_avg1_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae_avg1_cesm126, "Forecasted Distribution 2015 - 2039 \n CESM SSP126")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_cesm_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_yfslarvae6 <- pred_loop(2040:2044, yfs_larvae, 130, 
                              temps_cesm_ssp126,
                              salts_cesm_ssp126, 6,
                              larval_formula)
grids_yfslarvae7 <- pred_loop(2045:2049, yfs_larvae, 130, 
                              temps_cesm_ssp126,
                              salts_cesm_ssp126, 7,
                              larval_formula)
grids_yfslarvae8 <- pred_loop(2050:2054, yfs_larvae, 130, 
                              temps_cesm_ssp126,
                              salts_cesm_ssp126, 8,
                              larval_formula)
grids_yfslarvae9 <- pred_loop(2055:2059, yfs_larvae, 130, 
                              temps_cesm_ssp126,
                              salts_cesm_ssp126, 9,
                              larval_formula)
grids_yfslarvae10 <- pred_loop(2060:2064, yfs_larvae, 130, 
                               temps_cesm_ssp126,
                               salts_cesm_ssp126, 10,
                               larval_formula)
grids_yfslarvae11 <- pred_loop(2065:2069, yfs_larvae, 130, 
                               temps_cesm_ssp126,
                               salts_cesm_ssp126, 11,
                               larval_formula)

# Combine into one data frame
df_yfslarvae2 <- list(grids_yfslarvae6[[1]], grids_yfslarvae6[[2]], grids_yfslarvae6[[3]], 
                      grids_yfslarvae6[[4]], grids_yfslarvae6[[5]], grids_yfslarvae7[[1]], 
                      grids_yfslarvae7[[2]], grids_yfslarvae7[[3]], grids_yfslarvae7[[4]],
                      grids_yfslarvae7[[5]], grids_yfslarvae8[[1]], grids_yfslarvae8[[2]], 
                      grids_yfslarvae8[[3]], grids_yfslarvae8[[4]], grids_yfslarvae8[[5]],
                      grids_yfslarvae9[[1]], grids_yfslarvae9[[2]], grids_yfslarvae9[[3]], 
                      grids_yfslarvae9[[4]], grids_yfslarvae9[[5]], grids_yfslarvae10[[1]], 
                      grids_yfslarvae10[[2]], grids_yfslarvae10[[3]], grids_yfslarvae10[[4]],
                      grids_yfslarvae11[[5]], grids_yfslarvae11[[1]], grids_yfslarvae11[[2]],
                      grids_yfslarvae11[[3]], grids_yfslarvae11[[4]], grids_yfslarvae11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_yfslarvae2), fixed = T)
df_yfslarvae_avg2_cesm126 <- data.frame(lat = df_yfslarvae2$lat, 
                                        lon = df_yfslarvae2$lon, 
                                        dist = df_yfslarvae2$dist,
                                        avg_pred = rowSums(df_yfslarvae2[, x])/30)
saveRDS(df_yfslarvae_avg2_cesm126, file = here("data", "df_yfslarvae_avg2_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae_avg2_cesm126, "Forecasted Distribution 2040 - 2069 \n CESM SSP126")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_cesm_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_yfslarvae12 <- pred_loop(2070:2074, yfs_larvae, 130,
                               temps_cesm_ssp126,
                               salts_cesm_ssp126, 12,
                               larval_formula)
grids_yfslarvae13 <- pred_loop(2075:2079, yfs_larvae, 130,
                               temps_cesm_ssp126,
                               salts_cesm_ssp126, 13,
                               larval_formula)
grids_yfslarvae14 <- pred_loop(2080:2084, yfs_larvae, 130,
                               temps_cesm_ssp126,
                               salts_cesm_ssp126, 14,
                               larval_formula)
grids_yfslarvae15 <- pred_loop(2085:2089, yfs_larvae, 130,
                               temps_cesm_ssp126,
                               salts_cesm_ssp126, 15,
                               larval_formula)
grids_yfslarvae16 <- pred_loop(2090:2094, yfs_larvae, 130,
                               temps_cesm_ssp126,
                               salts_cesm_ssp126, 16,
                               larval_formula)
grids_yfslarvae17 <- pred_loop(2095:2099, yfs_larvae, 130, 
                               temps_cesm_ssp126,
                               salts_cesm_ssp126, 17,
                               larval_formula)

# Combine into one data frame
df_yfslarvae3 <- list(grids_yfslarvae12[[1]], grids_yfslarvae12[[2]], grids_yfslarvae12[[3]], 
                      grids_yfslarvae12[[4]], grids_yfslarvae12[[5]], grids_yfslarvae13[[1]], 
                      grids_yfslarvae13[[2]], grids_yfslarvae13[[3]], grids_yfslarvae13[[4]],
                      grids_yfslarvae13[[5]], grids_yfslarvae14[[1]], grids_yfslarvae14[[2]], 
                      grids_yfslarvae14[[3]], grids_yfslarvae14[[4]], grids_yfslarvae14[[5]],
                      grids_yfslarvae15[[1]], grids_yfslarvae15[[2]], grids_yfslarvae15[[3]], 
                      grids_yfslarvae15[[4]], grids_yfslarvae15[[5]], grids_yfslarvae16[[1]], 
                      grids_yfslarvae16[[2]], grids_yfslarvae16[[3]], grids_yfslarvae16[[4]],
                      grids_yfslarvae17[[5]], grids_yfslarvae17[[1]], grids_yfslarvae17[[2]],
                      grids_yfslarvae17[[3]], grids_yfslarvae17[[4]], grids_yfslarvae17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_yfslarvae3), fixed = T)
df_yfslarvae_avg3_cesm126 <- data.frame(lat = df_yfslarvae3$lat, 
                                        lon = df_yfslarvae3$lon, 
                                        dist = df_yfslarvae3$dist,
                                        avg_pred = rowSums(df_yfslarvae3[, x])/30)
saveRDS(df_yfslarvae_avg3_cesm126, file = here("data", "df_yfslarvae_avg3_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae_avg3_cesm126, "Forecasted Distribution 2070 - 2099 \n CESM SSP126")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_cesm_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

##### CESM 585 -----------------------------------------------------------------------------------------------------------------
temps_cesm_ssp585 <- readRDS(here('data', 'temps_cesm_ssp585.rds'))
salts_cesm_ssp585 <- readRDS(here('data', 'salts_cesm_ssp585.rds'))

## 2015 - 2039
grids_yfslarvae1 <- pred_loop(2015:2019, yfs_larvae, 130, 
                              temps_cesm_ssp585,
                              salts_cesm_ssp585, 1,
                              larval_formula)
grids_yfslarvae2 <- pred_loop(2020:2024, yfs_larvae, 130, 
                              temps_cesm_ssp585,
                              salts_cesm_ssp585, 2,
                              larval_formula)
grids_yfslarvae3 <- pred_loop(2025:2029, yfs_larvae, 130,
                              temps_cesm_ssp585,
                              salts_cesm_ssp585, 3,
                              larval_formula)
grids_yfslarvae4 <- pred_loop(2030:2034, yfs_larvae, 130, 
                              temps_cesm_ssp585,
                              salts_cesm_ssp585, 4,
                              larval_formula)
grids_yfslarvae5 <- pred_loop(2035:2039, yfs_larvae, 130, 
                              temps_cesm_ssp585,
                              salts_cesm_ssp585, 5,
                              larval_formula)

# Combine into one data frame
df_yfslarvae4 <- list(grids_yfslarvae1[[1]], grids_yfslarvae1[[2]], grids_yfslarvae1[[3]], 
                      grids_yfslarvae1[[4]], grids_yfslarvae1[[5]], grids_yfslarvae2[[1]], 
                      grids_yfslarvae2[[2]], grids_yfslarvae2[[3]], grids_yfslarvae2[[4]],
                      grids_yfslarvae2[[5]], grids_yfslarvae3[[1]], grids_yfslarvae3[[2]], 
                      grids_yfslarvae3[[3]], grids_yfslarvae3[[4]], grids_yfslarvae3[[5]],
                      grids_yfslarvae4[[1]], grids_yfslarvae4[[2]], grids_yfslarvae4[[3]], 
                      grids_yfslarvae4[[4]], grids_yfslarvae4[[5]], grids_yfslarvae5[[1]], 
                      grids_yfslarvae5[[2]], grids_yfslarvae5[[3]], grids_yfslarvae5[[4]],
                      grids_yfslarvae5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_yfslarvae4), fixed = T)
df_yfslarvae_avg4_cesm585 <- data.frame(lat = df_yfslarvae4$lat, 
                                        lon = df_yfslarvae4$lon, 
                                        dist = df_yfslarvae4$dist,
                                        avg_pred = rowSums(df_yfslarvae4[, x])/25)
saveRDS(df_yfslarvae_avg4_cesm585, file = here("data", "df_yfslarvae_avg4_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae_avg4_cesm585, "Forecasted Distribution 2015 - 2039 \n CESM SSP585")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_cesm_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_yfslarvae6 <- pred_loop(2040:2044, yfs_larvae, 130, 
                              temps_cesm_ssp585,
                              salts_cesm_ssp585, 6,
                              larval_formula)
grids_yfslarvae7 <- pred_loop(2045:2049, yfs_larvae, 130, 
                              temps_cesm_ssp585,
                              salts_cesm_ssp585, 7,
                              larval_formula)
grids_yfslarvae8 <- pred_loop(2050:2054, yfs_larvae, 130, 
                              temps_cesm_ssp585,
                              salts_cesm_ssp585, 8,
                              larval_formula)
grids_yfslarvae9 <- pred_loop(2055:2059, yfs_larvae, 130, 
                              temps_cesm_ssp585,
                              salts_cesm_ssp585, 9,
                              larval_formula)
grids_yfslarvae10 <- pred_loop(2060:2064, yfs_larvae, 130, 
                               temps_cesm_ssp585,
                               salts_cesm_ssp585, 10,
                               larval_formula)
grids_yfslarvae11 <- pred_loop(2065:2069, yfs_larvae, 130, 
                               temps_cesm_ssp585,
                               salts_cesm_ssp585, 11,
                               larval_formula)

# Combine into one data frame
df_yfslarvae5 <- list(grids_yfslarvae6[[1]], grids_yfslarvae6[[2]], grids_yfslarvae6[[3]], 
                      grids_yfslarvae6[[4]], grids_yfslarvae6[[5]], grids_yfslarvae7[[1]], 
                      grids_yfslarvae7[[2]], grids_yfslarvae7[[3]], grids_yfslarvae7[[4]],
                      grids_yfslarvae7[[5]], grids_yfslarvae8[[1]], grids_yfslarvae8[[2]], 
                      grids_yfslarvae8[[3]], grids_yfslarvae8[[4]], grids_yfslarvae8[[5]],
                      grids_yfslarvae9[[1]], grids_yfslarvae9[[2]], grids_yfslarvae9[[3]], 
                      grids_yfslarvae9[[4]], grids_yfslarvae9[[5]], grids_yfslarvae10[[1]], 
                      grids_yfslarvae10[[2]], grids_yfslarvae10[[3]], grids_yfslarvae10[[4]],
                      grids_yfslarvae11[[5]], grids_yfslarvae11[[1]], grids_yfslarvae11[[2]],
                      grids_yfslarvae11[[3]], grids_yfslarvae11[[4]], grids_yfslarvae11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_yfslarvae5), fixed = T)
df_yfslarvae_avg5_cesm585 <- data.frame(lat = df_yfslarvae5$lat, 
                                        lon = df_yfslarvae5$lon, 
                                        dist = df_yfslarvae5$dist,
                                        avg_pred = rowSums(df_yfslarvae5[, x])/30)
saveRDS(df_yfslarvae_avg5_cesm585, file = here("data", "df_yfslarvae_avg5_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae_avg5_cesm585, "Forecasted Distribution 2040 - 2069 \n CESM SSP585")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_cesm_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_yfslarvae12 <- pred_loop(2070:2074, yfs_larvae, 130,
                               temps_cesm_ssp585,
                               salts_cesm_ssp585, 12,
                               larval_formula)
grids_yfslarvae13 <- pred_loop(2075:2079, yfs_larvae, 130,
                               temps_cesm_ssp585,
                               salts_cesm_ssp585, 13,
                               larval_formula)
grids_yfslarvae14 <- pred_loop(2080:2084, yfs_larvae, 130,
                               temps_cesm_ssp585,
                               salts_cesm_ssp585, 14,
                               larval_formula)
grids_yfslarvae15 <- pred_loop(2085:2089, yfs_larvae, 130,
                               temps_cesm_ssp585,
                               salts_cesm_ssp585, 15,
                               larval_formula)
grids_yfslarvae16 <- pred_loop(2090:2094, yfs_larvae, 130,
                               temps_cesm_ssp585,
                               salts_cesm_ssp585, 16,
                               larval_formula)
grids_yfslarvae17 <- pred_loop(2095:2099, yfs_larvae, 130, 
                               temps_cesm_ssp585,
                               salts_cesm_ssp585, 17,
                               larval_formula)

# Combine into one data frame
df_yfslarvae6 <- list(grids_yfslarvae12[[1]], grids_yfslarvae12[[2]], grids_yfslarvae12[[3]], 
                      grids_yfslarvae12[[4]], grids_yfslarvae12[[5]], grids_yfslarvae13[[1]], 
                      grids_yfslarvae13[[2]], grids_yfslarvae13[[3]], grids_yfslarvae13[[4]],
                      grids_yfslarvae13[[5]], grids_yfslarvae14[[1]], grids_yfslarvae14[[2]], 
                      grids_yfslarvae14[[3]], grids_yfslarvae14[[4]], grids_yfslarvae14[[5]],
                      grids_yfslarvae15[[1]], grids_yfslarvae15[[2]], grids_yfslarvae15[[3]], 
                      grids_yfslarvae15[[4]], grids_yfslarvae15[[5]], grids_yfslarvae16[[1]], 
                      grids_yfslarvae16[[2]], grids_yfslarvae16[[3]], grids_yfslarvae16[[4]],
                      grids_yfslarvae17[[5]], grids_yfslarvae17[[1]], grids_yfslarvae17[[2]],
                      grids_yfslarvae17[[3]], grids_yfslarvae17[[4]], grids_yfslarvae17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_yfslarvae6), fixed = T)
df_yfslarvae_avg6_cesm585 <- data.frame(lat = df_yfslarvae6$lat, 
                                        lon = df_yfslarvae6$lon, 
                                        dist = df_yfslarvae6$dist,
                                        avg_pred = rowSums(df_yfslarvae6[, x])/30)
saveRDS(df_yfslarvae_avg6_cesm585, file = here("data", "df_yfslarvae_avg6_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae_avg6_cesm585, "Forecasted Distribution 2070 - 2099 \n CESM SSP585")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_cesm_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### GFDL 126 ----------------------------------------------------------------------------------------------------------------------------
temps_gfdl_ssp126 <- readRDS(here('data', 'temps_gfdl_ssp126.rds'))
salts_gfdl_ssp126 <- readRDS(here('data', 'salts_gfdl_ssp126.rds'))

## 2015 - 2039
grids_yfslarvae1 <- pred_loop(2015:2019, yfs_larvae, 130, 
                              temps_gfdl_ssp126,
                              salts_gfdl_ssp126, 1,
                              larval_formula)
grids_yfslarvae2 <- pred_loop(2020:2024, yfs_larvae, 130, 
                              temps_gfdl_ssp126,
                              salts_gfdl_ssp126, 2,
                              larval_formula)
grids_yfslarvae3 <- pred_loop(2025:2029, yfs_larvae, 130,
                              temps_gfdl_ssp126,
                              salts_gfdl_ssp126, 3,
                              larval_formula)
grids_yfslarvae4 <- pred_loop(2030:2034, yfs_larvae, 130, 
                              temps_gfdl_ssp126,
                              salts_gfdl_ssp126, 4,
                              larval_formula)
grids_yfslarvae5 <- pred_loop(2035:2039, yfs_larvae, 130, 
                              temps_gfdl_ssp126,
                              salts_gfdl_ssp126, 5,
                              larval_formula)

# Combine into one data frame
df_yfslarvae1 <- list(grids_yfslarvae1[[1]], grids_yfslarvae1[[2]], grids_yfslarvae1[[3]], 
                      grids_yfslarvae1[[4]], grids_yfslarvae1[[5]], grids_yfslarvae2[[1]], 
                      grids_yfslarvae2[[2]], grids_yfslarvae2[[3]], grids_yfslarvae2[[4]],
                      grids_yfslarvae2[[5]], grids_yfslarvae3[[1]], grids_yfslarvae3[[2]], 
                      grids_yfslarvae3[[3]], grids_yfslarvae3[[4]], grids_yfslarvae3[[5]],
                      grids_yfslarvae4[[1]], grids_yfslarvae4[[2]], grids_yfslarvae4[[3]], 
                      grids_yfslarvae4[[4]], grids_yfslarvae4[[5]], grids_yfslarvae5[[1]], 
                      grids_yfslarvae5[[2]], grids_yfslarvae5[[3]], grids_yfslarvae5[[4]],
                      grids_yfslarvae5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_yfslarvae1), fixed = T)
df_yfslarvae_avg1_gfdl126 <- data.frame(lat = df_yfslarvae1$lat, 
                                        lon = df_yfslarvae1$lon, 
                                        dist = df_yfslarvae1$dist,
                                        avg_pred = rowSums(df_yfslarvae1[, x])/25)
saveRDS(df_yfslarvae_avg1_gfdl126, file = here("data", "df_yfslarvae_avg1_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae_avg1_gfdl126, "Forecasted Distribution 2015 - 2039 \n GFDL SSP126")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_gfdl_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_yfslarvae6 <- pred_loop(2040:2044, yfs_larvae, 130, 
                              temps_gfdl_ssp126,
                              salts_gfdl_ssp126, 6,
                              larval_formula)
grids_yfslarvae7 <- pred_loop(2045:2049, yfs_larvae, 130, 
                              temps_gfdl_ssp126,
                              salts_gfdl_ssp126, 7,
                              larval_formula)
grids_yfslarvae8 <- pred_loop(2050:2054, yfs_larvae, 130, 
                              temps_gfdl_ssp126,
                              salts_gfdl_ssp126, 8,
                              larval_formula)
grids_yfslarvae9 <- pred_loop(2055:2059, yfs_larvae, 130, 
                              temps_gfdl_ssp126,
                              salts_gfdl_ssp126, 9,
                              larval_formula)
grids_yfslarvae10 <- pred_loop(2060:2064, yfs_larvae, 130, 
                               temps_gfdl_ssp126,
                               salts_gfdl_ssp126, 10,
                               larval_formula)
grids_yfslarvae11 <- pred_loop(2065:2069, yfs_larvae, 130, 
                               temps_gfdl_ssp126,
                               salts_gfdl_ssp126, 11,
                               larval_formula)

# Combine into one data frame
df_yfslarvae2 <- list(grids_yfslarvae6[[1]], grids_yfslarvae6[[2]], grids_yfslarvae6[[3]], 
                      grids_yfslarvae6[[4]], grids_yfslarvae6[[5]], grids_yfslarvae7[[1]], 
                      grids_yfslarvae7[[2]], grids_yfslarvae7[[3]], grids_yfslarvae7[[4]],
                      grids_yfslarvae7[[5]], grids_yfslarvae8[[1]], grids_yfslarvae8[[2]], 
                      grids_yfslarvae8[[3]], grids_yfslarvae8[[4]], grids_yfslarvae8[[5]],
                      grids_yfslarvae9[[1]], grids_yfslarvae9[[2]], grids_yfslarvae9[[3]], 
                      grids_yfslarvae9[[4]], grids_yfslarvae9[[5]], grids_yfslarvae10[[1]], 
                      grids_yfslarvae10[[2]], grids_yfslarvae10[[3]], grids_yfslarvae10[[4]],
                      grids_yfslarvae11[[5]], grids_yfslarvae11[[1]], grids_yfslarvae11[[2]],
                      grids_yfslarvae11[[3]], grids_yfslarvae11[[4]], grids_yfslarvae11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_yfslarvae2), fixed = T)
df_yfslarvae_avg2_gfdl126 <- data.frame(lat = df_yfslarvae2$lat, 
                                        lon = df_yfslarvae2$lon, 
                                        dist = df_yfslarvae2$dist,
                                        avg_pred = rowSums(df_yfslarvae2[, x])/30)
saveRDS(df_yfslarvae_avg2_gfdl126, file = here("data", "df_yfslarvae_avg2_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae_avg2_gfdl126, "Forecasted Distribution 2040 - 2069 \n GFDL SSP126")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_gfdl_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_yfslarvae12 <- pred_loop(2070:2074, yfs_larvae, 130,
                               temps_gfdl_ssp126,
                               salts_gfdl_ssp126, 12,
                               larval_formula)
grids_yfslarvae13 <- pred_loop(2075:2079, yfs_larvae, 130,
                               temps_gfdl_ssp126,
                               salts_gfdl_ssp126, 13,
                               larval_formula)
grids_yfslarvae14 <- pred_loop(2080:2084, yfs_larvae, 130,
                               temps_gfdl_ssp126,
                               salts_gfdl_ssp126, 14,
                               larval_formula)
grids_yfslarvae15 <- pred_loop(2085:2089, yfs_larvae, 130,
                               temps_gfdl_ssp126,
                               salts_gfdl_ssp126, 15,
                               larval_formula)
grids_yfslarvae16 <- pred_loop(2090:2094, yfs_larvae, 130,
                               temps_gfdl_ssp126,
                               salts_gfdl_ssp126, 16,
                               larval_formula)
grids_yfslarvae17 <- pred_loop(2095:2099, yfs_larvae, 130, 
                               temps_gfdl_ssp126,
                               salts_gfdl_ssp126, 17,
                               larval_formula)

# Combine into one data frame
df_yfslarvae3 <- list(grids_yfslarvae12[[1]], grids_yfslarvae12[[2]], grids_yfslarvae12[[3]], 
                      grids_yfslarvae12[[4]], grids_yfslarvae12[[5]], grids_yfslarvae13[[1]], 
                      grids_yfslarvae13[[2]], grids_yfslarvae13[[3]], grids_yfslarvae13[[4]],
                      grids_yfslarvae13[[5]], grids_yfslarvae14[[1]], grids_yfslarvae14[[2]], 
                      grids_yfslarvae14[[3]], grids_yfslarvae14[[4]], grids_yfslarvae14[[5]],
                      grids_yfslarvae15[[1]], grids_yfslarvae15[[2]], grids_yfslarvae15[[3]], 
                      grids_yfslarvae15[[4]], grids_yfslarvae15[[5]], grids_yfslarvae16[[1]], 
                      grids_yfslarvae16[[2]], grids_yfslarvae16[[3]], grids_yfslarvae16[[4]],
                      grids_yfslarvae17[[5]], grids_yfslarvae17[[1]], grids_yfslarvae17[[2]],
                      grids_yfslarvae17[[3]], grids_yfslarvae17[[4]], grids_yfslarvae17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_yfslarvae3), fixed = T)
df_yfslarvae_avg3_gfdl126 <- data.frame(lat = df_yfslarvae3$lat, 
                                        lon = df_yfslarvae3$lon, 
                                        dist = df_yfslarvae3$dist,
                                        avg_pred = rowSums(df_yfslarvae3[, x])/30)
saveRDS(df_yfslarvae_avg3_gfdl126, file = here("data", "df_yfslarvae_avg3_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae_avg3_gfdl126, "Forecasted Distribution 2070 - 2099 \n GFDL SSP126")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_gfdl_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

##### GFDL 585 -----------------------------------------------------------------------------------------------------------------
temps_gfdl_ssp585 <- readRDS(here('data', 'temps_gfdl_ssp585.rds'))
salts_gfdl_ssp585 <- readRDS(here('data', 'salts_gfdl_ssp585.rds'))

## 2015 - 2039
grids_yfslarvae1 <- pred_loop(2015:2019, yfs_larvae, 130, 
                              temps_gfdl_ssp585,
                              salts_gfdl_ssp585, 1,
                              larval_formula)
grids_yfslarvae2 <- pred_loop(2020:2024, yfs_larvae, 130, 
                              temps_gfdl_ssp585,
                              salts_gfdl_ssp585, 2,
                              larval_formula)
grids_yfslarvae3 <- pred_loop(2025:2029, yfs_larvae, 130,
                              temps_gfdl_ssp585,
                              salts_gfdl_ssp585, 3,
                              larval_formula)
grids_yfslarvae4 <- pred_loop(2030:2034, yfs_larvae, 130, 
                              temps_gfdl_ssp585,
                              salts_gfdl_ssp585, 4,
                              larval_formula)
grids_yfslarvae5 <- pred_loop(2035:2039, yfs_larvae, 130, 
                              temps_gfdl_ssp585,
                              salts_gfdl_ssp585, 5,
                              larval_formula)

# Combine into one data frame
df_yfslarvae4 <- list(grids_yfslarvae1[[1]], grids_yfslarvae1[[2]], grids_yfslarvae1[[3]], 
                      grids_yfslarvae1[[4]], grids_yfslarvae1[[5]], grids_yfslarvae2[[1]], 
                      grids_yfslarvae2[[2]], grids_yfslarvae2[[3]], grids_yfslarvae2[[4]],
                      grids_yfslarvae2[[5]], grids_yfslarvae3[[1]], grids_yfslarvae3[[2]], 
                      grids_yfslarvae3[[3]], grids_yfslarvae3[[4]], grids_yfslarvae3[[5]],
                      grids_yfslarvae4[[1]], grids_yfslarvae4[[2]], grids_yfslarvae4[[3]], 
                      grids_yfslarvae4[[4]], grids_yfslarvae4[[5]], grids_yfslarvae5[[1]], 
                      grids_yfslarvae5[[2]], grids_yfslarvae5[[3]], grids_yfslarvae5[[4]],
                      grids_yfslarvae5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_yfslarvae4), fixed = T)
df_yfslarvae_avg4_gfdl585 <- data.frame(lat = df_yfslarvae4$lat, 
                                        lon = df_yfslarvae4$lon, 
                                        dist = df_yfslarvae4$dist,
                                        avg_pred = rowSums(df_yfslarvae4[, x])/25)
saveRDS(df_yfslarvae_avg4_gfdl585, file = here("data", "df_yfslarvae_avg4_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae_avg4_gfdl585, "Forecasted Distribution 2015 - 2039 \n GFDL SSP585")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_gfdl_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_yfslarvae6 <- pred_loop(2040:2044, yfs_larvae, 130, 
                              temps_gfdl_ssp585,
                              salts_gfdl_ssp585, 6,
                              larval_formula)
grids_yfslarvae7 <- pred_loop(2045:2049, yfs_larvae, 130, 
                              temps_gfdl_ssp585,
                              salts_gfdl_ssp585, 7,
                              larval_formula)
grids_yfslarvae8 <- pred_loop(2050:2054, yfs_larvae, 130, 
                              temps_gfdl_ssp585,
                              salts_gfdl_ssp585, 8,
                              larval_formula)
grids_yfslarvae9 <- pred_loop(2055:2059, yfs_larvae, 130, 
                              temps_gfdl_ssp585,
                              salts_gfdl_ssp585, 9,
                              larval_formula)
grids_yfslarvae10 <- pred_loop(2060:2064, yfs_larvae, 130, 
                               temps_gfdl_ssp585,
                               salts_gfdl_ssp585, 10,
                               larval_formula)
grids_yfslarvae11 <- pred_loop(2065:2069, yfs_larvae, 130, 
                               temps_gfdl_ssp585,
                               salts_gfdl_ssp585, 11,
                               larval_formula)

# Combine into one data frame
df_yfslarvae5 <- list(grids_yfslarvae6[[1]], grids_yfslarvae6[[2]], grids_yfslarvae6[[3]], 
                      grids_yfslarvae6[[4]], grids_yfslarvae6[[5]], grids_yfslarvae7[[1]], 
                      grids_yfslarvae7[[2]], grids_yfslarvae7[[3]], grids_yfslarvae7[[4]],
                      grids_yfslarvae7[[5]], grids_yfslarvae8[[1]], grids_yfslarvae8[[2]], 
                      grids_yfslarvae8[[3]], grids_yfslarvae8[[4]], grids_yfslarvae8[[5]],
                      grids_yfslarvae9[[1]], grids_yfslarvae9[[2]], grids_yfslarvae9[[3]], 
                      grids_yfslarvae9[[4]], grids_yfslarvae9[[5]], grids_yfslarvae10[[1]], 
                      grids_yfslarvae10[[2]], grids_yfslarvae10[[3]], grids_yfslarvae10[[4]],
                      grids_yfslarvae11[[5]], grids_yfslarvae11[[1]], grids_yfslarvae11[[2]],
                      grids_yfslarvae11[[3]], grids_yfslarvae11[[4]], grids_yfslarvae11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_yfslarvae5), fixed = T)
df_yfslarvae_avg5_gfdl585 <- data.frame(lat = df_yfslarvae5$lat, 
                                        lon = df_yfslarvae5$lon, 
                                        dist = df_yfslarvae5$dist,
                                        avg_pred = rowSums(df_yfslarvae5[, x])/30)
saveRDS(df_yfslarvae_avg5_gfdl585, file = here("data", "df_yfslarvae_avg5_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae_avg5_gfdl585, "Forecasted Distribution 2040 - 2069 \n GFDL SSP585")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_gfdl_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_yfslarvae12 <- pred_loop(2070:2074, yfs_larvae, 130,
                               temps_gfdl_ssp585,
                               salts_gfdl_ssp585, 12,
                               larval_formula)
grids_yfslarvae13 <- pred_loop(2075:2079, yfs_larvae, 130,
                               temps_gfdl_ssp585,
                               salts_gfdl_ssp585, 13,
                               larval_formula)
grids_yfslarvae14 <- pred_loop(2080:2084, yfs_larvae, 130,
                               temps_gfdl_ssp585,
                               salts_gfdl_ssp585, 14,
                               larval_formula)
grids_yfslarvae15 <- pred_loop(2085:2089, yfs_larvae, 130,
                               temps_gfdl_ssp585,
                               salts_gfdl_ssp585, 15,
                               larval_formula)
grids_yfslarvae16 <- pred_loop(2090:2094, yfs_larvae, 130,
                               temps_gfdl_ssp585,
                               salts_gfdl_ssp585, 16,
                               larval_formula)
grids_yfslarvae17 <- pred_loop(2095:2099, yfs_larvae, 130, 
                               temps_gfdl_ssp585,
                               salts_gfdl_ssp585, 17,
                               larval_formula)

# Combine into one data frame
df_yfslarvae6 <- list(grids_yfslarvae12[[1]], grids_yfslarvae12[[2]], grids_yfslarvae12[[3]], 
                      grids_yfslarvae12[[4]], grids_yfslarvae12[[5]], grids_yfslarvae13[[1]], 
                      grids_yfslarvae13[[2]], grids_yfslarvae13[[3]], grids_yfslarvae13[[4]],
                      grids_yfslarvae13[[5]], grids_yfslarvae14[[1]], grids_yfslarvae14[[2]], 
                      grids_yfslarvae14[[3]], grids_yfslarvae14[[4]], grids_yfslarvae14[[5]],
                      grids_yfslarvae15[[1]], grids_yfslarvae15[[2]], grids_yfslarvae15[[3]], 
                      grids_yfslarvae15[[4]], grids_yfslarvae15[[5]], grids_yfslarvae16[[1]], 
                      grids_yfslarvae16[[2]], grids_yfslarvae16[[3]], grids_yfslarvae16[[4]],
                      grids_yfslarvae17[[5]], grids_yfslarvae17[[1]], grids_yfslarvae17[[2]],
                      grids_yfslarvae17[[3]], grids_yfslarvae17[[4]], grids_yfslarvae17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_yfslarvae6), fixed = T)
df_yfslarvae_avg6_gfdl585 <- data.frame(lat = df_yfslarvae6$lat, 
                                        lon = df_yfslarvae6$lon, 
                                        dist = df_yfslarvae6$dist,
                                        avg_pred = rowSums(df_yfslarvae6[, x])/30)
saveRDS(df_yfslarvae_avg6_gfdl585, file = here("data", "df_yfslarvae_avg6_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae_avg6_gfdl585, "Forecasted Distribution 2070 - 2099 \n GFDL SSP585")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_gfdl_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### MIROC 126 ----------------------------------------------------------------------------------------------------------------------------
temps_miroc_ssp126 <- readRDS(here('data', 'temps_miroc_ssp126.rds'))
salts_miroc_ssp126 <- readRDS(here('data', 'salts_miroc_ssp126.rds'))

## 2015 - 2039
grids_yfslarvae1 <- pred_loop(2015:2019, yfs_larvae, 130, 
                              temps_miroc_ssp126,
                              salts_miroc_ssp126, 1,
                              larval_formula)
grids_yfslarvae2 <- pred_loop(2020:2024, yfs_larvae, 130, 
                              temps_miroc_ssp126,
                              salts_miroc_ssp126, 2,
                              larval_formula)
grids_yfslarvae3 <- pred_loop(2025:2029, yfs_larvae, 130,
                              temps_miroc_ssp126,
                              salts_miroc_ssp126, 3,
                              larval_formula)
grids_yfslarvae4 <- pred_loop(2030:2034, yfs_larvae, 130, 
                              temps_miroc_ssp126,
                              salts_miroc_ssp126, 4,
                              larval_formula)
grids_yfslarvae5 <- pred_loop(2035:2039, yfs_larvae, 130, 
                              temps_miroc_ssp126,
                              salts_miroc_ssp126, 5,
                              larval_formula)

# Combine into one data frame
df_yfslarvae1 <- list(grids_yfslarvae1[[1]], grids_yfslarvae1[[2]], grids_yfslarvae1[[3]], 
                      grids_yfslarvae1[[4]], grids_yfslarvae1[[5]], grids_yfslarvae2[[1]], 
                      grids_yfslarvae2[[2]], grids_yfslarvae2[[3]], grids_yfslarvae2[[4]],
                      grids_yfslarvae2[[5]], grids_yfslarvae3[[1]], grids_yfslarvae3[[2]], 
                      grids_yfslarvae3[[3]], grids_yfslarvae3[[4]], grids_yfslarvae3[[5]],
                      grids_yfslarvae4[[1]], grids_yfslarvae4[[2]], grids_yfslarvae4[[3]], 
                      grids_yfslarvae4[[4]], grids_yfslarvae4[[5]], grids_yfslarvae5[[1]], 
                      grids_yfslarvae5[[2]], grids_yfslarvae5[[3]], grids_yfslarvae5[[4]],
                      grids_yfslarvae5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_yfslarvae1), fixed = T)
df_yfslarvae_avg1_miroc126 <- data.frame(lat = df_yfslarvae1$lat, 
                                         lon = df_yfslarvae1$lon, 
                                         dist = df_yfslarvae1$dist,
                                         avg_pred = rowSums(df_yfslarvae1[, x])/25)
saveRDS(df_yfslarvae_avg1_miroc126, file = here("data", "df_yfslarvae_avg1_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae_avg1_miroc126, "Forecasted Distribution 2015 - 2039 \n MIROC SSP126")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_miroc_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_yfslarvae6 <- pred_loop(2040:2044, yfs_larvae, 130, 
                              temps_miroc_ssp126,
                              salts_miroc_ssp126, 6,
                              larval_formula)
grids_yfslarvae7 <- pred_loop(2045:2049, yfs_larvae, 130, 
                              temps_miroc_ssp126,
                              salts_miroc_ssp126, 7,
                              larval_formula)
grids_yfslarvae8 <- pred_loop(2050:2054, yfs_larvae, 130, 
                              temps_miroc_ssp126,
                              salts_miroc_ssp126, 8,
                              larval_formula)
grids_yfslarvae9 <- pred_loop(2055:2059, yfs_larvae, 130, 
                              temps_miroc_ssp126,
                              salts_miroc_ssp126, 9,
                              larval_formula)
grids_yfslarvae10 <- pred_loop(2060:2064, yfs_larvae, 130, 
                               temps_miroc_ssp126,
                               salts_miroc_ssp126, 10,
                               larval_formula)
grids_yfslarvae11 <- pred_loop(2065:2069, yfs_larvae, 130, 
                               temps_miroc_ssp126,
                               salts_miroc_ssp126, 11,
                               larval_formula)

# Combine into one data frame
df_yfslarvae2 <- list(grids_yfslarvae6[[1]], grids_yfslarvae6[[2]], grids_yfslarvae6[[3]], 
                      grids_yfslarvae6[[4]], grids_yfslarvae6[[5]], grids_yfslarvae7[[1]], 
                      grids_yfslarvae7[[2]], grids_yfslarvae7[[3]], grids_yfslarvae7[[4]],
                      grids_yfslarvae7[[5]], grids_yfslarvae8[[1]], grids_yfslarvae8[[2]], 
                      grids_yfslarvae8[[3]], grids_yfslarvae8[[4]], grids_yfslarvae8[[5]],
                      grids_yfslarvae9[[1]], grids_yfslarvae9[[2]], grids_yfslarvae9[[3]], 
                      grids_yfslarvae9[[4]], grids_yfslarvae9[[5]], grids_yfslarvae10[[1]], 
                      grids_yfslarvae10[[2]], grids_yfslarvae10[[3]], grids_yfslarvae10[[4]],
                      grids_yfslarvae11[[5]], grids_yfslarvae11[[1]], grids_yfslarvae11[[2]],
                      grids_yfslarvae11[[3]], grids_yfslarvae11[[4]], grids_yfslarvae11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_yfslarvae2), fixed = T)
df_yfslarvae_avg2_miroc126 <- data.frame(lat = df_yfslarvae2$lat, 
                                         lon = df_yfslarvae2$lon, 
                                         dist = df_yfslarvae2$dist,
                                         avg_pred = rowSums(df_yfslarvae2[, x])/30)
saveRDS(df_yfslarvae_avg2_miroc126, file = here("data", "df_yfslarvae_avg2_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae_avg2_miroc126, "Forecasted Distribution 2040 - 2069 \n MIROC SSP126")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_miroc_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_yfslarvae12 <- pred_loop(2070:2074, yfs_larvae, 130,
                               temps_miroc_ssp126,
                               salts_miroc_ssp126, 12,
                               larval_formula)
grids_yfslarvae13 <- pred_loop(2075:2079, yfs_larvae, 130,
                               temps_miroc_ssp126,
                               salts_miroc_ssp126, 13,
                               larval_formula)
grids_yfslarvae14 <- pred_loop(2080:2084, yfs_larvae, 130,
                               temps_miroc_ssp126,
                               salts_miroc_ssp126, 14,
                               larval_formula)
grids_yfslarvae15 <- pred_loop(2085:2089, yfs_larvae, 130,
                               temps_miroc_ssp126,
                               salts_miroc_ssp126, 15,
                               larval_formula)
grids_yfslarvae16 <- pred_loop(2090:2094, yfs_larvae, 130,
                               temps_miroc_ssp126,
                               salts_miroc_ssp126, 16,
                               larval_formula)
grids_yfslarvae17 <- pred_loop(2095:2099, yfs_larvae, 130, 
                               temps_miroc_ssp126,
                               salts_miroc_ssp126, 17,
                               larval_formula)

# Combine into one data frame
df_yfslarvae3 <- list(grids_yfslarvae12[[1]], grids_yfslarvae12[[2]], grids_yfslarvae12[[3]], 
                      grids_yfslarvae12[[4]], grids_yfslarvae12[[5]], grids_yfslarvae13[[1]], 
                      grids_yfslarvae13[[2]], grids_yfslarvae13[[3]], grids_yfslarvae13[[4]],
                      grids_yfslarvae13[[5]], grids_yfslarvae14[[1]], grids_yfslarvae14[[2]], 
                      grids_yfslarvae14[[3]], grids_yfslarvae14[[4]], grids_yfslarvae14[[5]],
                      grids_yfslarvae15[[1]], grids_yfslarvae15[[2]], grids_yfslarvae15[[3]], 
                      grids_yfslarvae15[[4]], grids_yfslarvae15[[5]], grids_yfslarvae16[[1]], 
                      grids_yfslarvae16[[2]], grids_yfslarvae16[[3]], grids_yfslarvae16[[4]],
                      grids_yfslarvae17[[5]], grids_yfslarvae17[[1]], grids_yfslarvae17[[2]],
                      grids_yfslarvae17[[3]], grids_yfslarvae17[[4]], grids_yfslarvae17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_yfslarvae3), fixed = T)
df_yfslarvae_avg3_miroc126 <- data.frame(lat = df_yfslarvae3$lat, 
                                         lon = df_yfslarvae3$lon, 
                                         dist = df_yfslarvae3$dist,
                                         avg_pred = rowSums(df_yfslarvae3[, x])/30)
saveRDS(df_yfslarvae_avg3_miroc126, file = here("data", "df_yfslarvae_avg3_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae_avg3_miroc126, "Forecasted Distribution 2070 - 2099 \n MIROC SSP126")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_miroc_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

##### MIROC 585 -----------------------------------------------------------------------------------------------------------------
temps_miroc_ssp585 <- readRDS(here('data', 'temps_miroc_ssp585.rds'))
salts_miroc_ssp585 <- readRDS(here('data', 'salts_miroc_ssp585.rds'))

## 2015 - 2039
grids_yfslarvae1 <- pred_loop(2015:2019, yfs_larvae, 130, 
                              temps_miroc_ssp585,
                              salts_miroc_ssp585, 1,
                              larval_formula)
grids_yfslarvae2 <- pred_loop(2020:2024, yfs_larvae, 130, 
                              temps_miroc_ssp585,
                              salts_miroc_ssp585, 2,
                              larval_formula)
grids_yfslarvae3 <- pred_loop(2025:2029, yfs_larvae, 130,
                              temps_miroc_ssp585,
                              salts_miroc_ssp585, 3,
                              larval_formula)
grids_yfslarvae4 <- pred_loop(2030:2034, yfs_larvae, 130, 
                              temps_miroc_ssp585,
                              salts_miroc_ssp585, 4,
                              larval_formula)
grids_yfslarvae5 <- pred_loop(2035:2039, yfs_larvae, 130, 
                              temps_miroc_ssp585,
                              salts_miroc_ssp585, 5,
                              larval_formula)

# Combine into one data frame
df_yfslarvae4 <- list(grids_yfslarvae1[[1]], grids_yfslarvae1[[2]], grids_yfslarvae1[[3]], 
                      grids_yfslarvae1[[4]], grids_yfslarvae1[[5]], grids_yfslarvae2[[1]], 
                      grids_yfslarvae2[[2]], grids_yfslarvae2[[3]], grids_yfslarvae2[[4]],
                      grids_yfslarvae2[[5]], grids_yfslarvae3[[1]], grids_yfslarvae3[[2]], 
                      grids_yfslarvae3[[3]], grids_yfslarvae3[[4]], grids_yfslarvae3[[5]],
                      grids_yfslarvae4[[1]], grids_yfslarvae4[[2]], grids_yfslarvae4[[3]], 
                      grids_yfslarvae4[[4]], grids_yfslarvae4[[5]], grids_yfslarvae5[[1]], 
                      grids_yfslarvae5[[2]], grids_yfslarvae5[[3]], grids_yfslarvae5[[4]],
                      grids_yfslarvae5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_yfslarvae4), fixed = T)
df_yfslarvae_avg4_miroc585 <- data.frame(lat = df_yfslarvae4$lat, 
                                         lon = df_yfslarvae4$lon, 
                                         dist = df_yfslarvae4$dist,
                                         avg_pred = rowSums(df_yfslarvae4[, x])/25)
saveRDS(df_yfslarvae_avg4_miroc585, file = here("data", "df_yfslarvae_avg4_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae_avg4_miroc585, "Forecasted Distribution 2015 - 2039 \n MIROC SSP585")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_miroc_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_yfslarvae6 <- pred_loop(2040:2044, yfs_larvae, 130, 
                              temps_miroc_ssp585,
                              salts_miroc_ssp585, 6,
                              larval_formula)
grids_yfslarvae7 <- pred_loop(2045:2049, yfs_larvae, 130, 
                              temps_miroc_ssp585,
                              salts_miroc_ssp585, 7,
                              larval_formula)
grids_yfslarvae8 <- pred_loop(2050:2054, yfs_larvae, 130, 
                              temps_miroc_ssp585,
                              salts_miroc_ssp585, 8,
                              larval_formula)
grids_yfslarvae9 <- pred_loop(2055:2059, yfs_larvae, 130, 
                              temps_miroc_ssp585,
                              salts_miroc_ssp585, 9,
                              larval_formula)
grids_yfslarvae10 <- pred_loop(2060:2064, yfs_larvae, 130, 
                               temps_miroc_ssp585,
                               salts_miroc_ssp585, 10,
                               larval_formula)
grids_yfslarvae11 <- pred_loop(2065:2069, yfs_larvae, 130, 
                               temps_miroc_ssp585,
                               salts_miroc_ssp585, 11,
                               larval_formula)

# Combine into one data frame
df_yfslarvae5 <- list(grids_yfslarvae6[[1]], grids_yfslarvae6[[2]], grids_yfslarvae6[[3]], 
                      grids_yfslarvae6[[4]], grids_yfslarvae6[[5]], grids_yfslarvae7[[1]], 
                      grids_yfslarvae7[[2]], grids_yfslarvae7[[3]], grids_yfslarvae7[[4]],
                      grids_yfslarvae7[[5]], grids_yfslarvae8[[1]], grids_yfslarvae8[[2]], 
                      grids_yfslarvae8[[3]], grids_yfslarvae8[[4]], grids_yfslarvae8[[5]],
                      grids_yfslarvae9[[1]], grids_yfslarvae9[[2]], grids_yfslarvae9[[3]], 
                      grids_yfslarvae9[[4]], grids_yfslarvae9[[5]], grids_yfslarvae10[[1]], 
                      grids_yfslarvae10[[2]], grids_yfslarvae10[[3]], grids_yfslarvae10[[4]],
                      grids_yfslarvae11[[5]], grids_yfslarvae11[[1]], grids_yfslarvae11[[2]],
                      grids_yfslarvae11[[3]], grids_yfslarvae11[[4]], grids_yfslarvae11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_yfslarvae5), fixed = T)
df_yfslarvae_avg5_miroc585 <- data.frame(lat = df_yfslarvae5$lat, 
                                         lon = df_yfslarvae5$lon, 
                                         dist = df_yfslarvae5$dist,
                                         avg_pred = rowSums(df_yfslarvae5[, x])/30)
saveRDS(df_yfslarvae_avg5_miroc585, file = here("data", "df_yfslarvae_avg5_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae_avg5_miroc585, "Forecasted Distribution 2040 - 2069 \n MIROC SSP585")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_miroc_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_yfslarvae12 <- pred_loop(2070:2074, yfs_larvae, 130,
                               temps_miroc_ssp585,
                               salts_miroc_ssp585, 12,
                               larval_formula)
grids_yfslarvae13 <- pred_loop(2075:2079, yfs_larvae, 130,
                               temps_miroc_ssp585,
                               salts_miroc_ssp585, 13,
                               larval_formula)
grids_yfslarvae14 <- pred_loop(2080:2084, yfs_larvae, 130,
                               temps_miroc_ssp585,
                               salts_miroc_ssp585, 14,
                               larval_formula)
grids_yfslarvae15 <- pred_loop(2085:2089, yfs_larvae, 130,
                               temps_miroc_ssp585,
                               salts_miroc_ssp585, 15,
                               larval_formula)
grids_yfslarvae16 <- pred_loop(2090:2094, yfs_larvae, 130,
                               temps_miroc_ssp585,
                               salts_miroc_ssp585, 16,
                               larval_formula)
grids_yfslarvae17 <- pred_loop(2095:2099, yfs_larvae, 130, 
                               temps_miroc_ssp585,
                               salts_miroc_ssp585, 17,
                               larval_formula)

# Combine into one data frame
df_yfslarvae6 <- list(grids_yfslarvae12[[1]], grids_yfslarvae12[[2]], grids_yfslarvae12[[3]], 
                      grids_yfslarvae12[[4]], grids_yfslarvae12[[5]], grids_yfslarvae13[[1]], 
                      grids_yfslarvae13[[2]], grids_yfslarvae13[[3]], grids_yfslarvae13[[4]],
                      grids_yfslarvae13[[5]], grids_yfslarvae14[[1]], grids_yfslarvae14[[2]], 
                      grids_yfslarvae14[[3]], grids_yfslarvae14[[4]], grids_yfslarvae14[[5]],
                      grids_yfslarvae15[[1]], grids_yfslarvae15[[2]], grids_yfslarvae15[[3]], 
                      grids_yfslarvae15[[4]], grids_yfslarvae15[[5]], grids_yfslarvae16[[1]], 
                      grids_yfslarvae16[[2]], grids_yfslarvae16[[3]], grids_yfslarvae16[[4]],
                      grids_yfslarvae17[[5]], grids_yfslarvae17[[1]], grids_yfslarvae17[[2]],
                      grids_yfslarvae17[[3]], grids_yfslarvae17[[4]], grids_yfslarvae17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_yfslarvae6), fixed = T)
df_yfslarvae_avg6_miroc585 <- data.frame(lat = df_yfslarvae6$lat, 
                                         lon = df_yfslarvae6$lon, 
                                         dist = df_yfslarvae6$dist,
                                         avg_pred = rowSums(df_yfslarvae6[, x])/30)
saveRDS(df_yfslarvae_avg6_miroc585, file = here("data", "df_yfslarvae_avg6_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae_avg6_miroc585, "Forecasted Distribution 2070 - 2099 \n MIROC SSP585")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_miroc_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

### Average Predictions ------------------------------------------------------------------------------------------------------------


#### 2015-2039 ---------------------------------------------------------------------------------------------------------------------
##### Eggs
df_yfsegg_avg1_cesm126 <- readRDS(here('data', 'df_yfsegg_avg1_cesm126.rds'))
df_yfsegg_avg4_cesm585 <- readRDS(here('data', 'df_yfsegg_avg4_cesm585.rds'))
df_yfsegg_avg1_gfdl126 <- readRDS(here('data', 'df_yfsegg_avg1_gfdl126.rds'))
df_yfsegg_avg4_gfdl585 <- readRDS(here('data', 'df_yfsegg_avg4_gfdl585.rds'))
df_yfsegg_avg1_miroc126 <- readRDS(here('data', 'df_yfsegg_avg1_miroc126.rds'))
df_yfsegg_avg4_miroc585 <- readRDS(here('data', 'df_yfsegg_avg4_miroc585.rds'))

df_yfsegg_merged1 <- list(df_yfsegg_avg1_cesm126, df_yfsegg_avg4_cesm585,
                          df_yfsegg_avg1_gfdl126, df_yfsegg_avg4_gfdl585,
                          df_yfsegg_avg1_miroc126, df_yfsegg_avg4_miroc585) %>%
  reduce(inner_join, by = c("lon", "lat", "dist"))

x <- grepl("pred", names(df_yfsegg_merged1), fixed = T)
df_yfsegg_final1 <- data.frame(lat = df_yfsegg_merged1$lat,
                               lon = df_yfsegg_merged1$lon,
                               avg_pred = (rowSums(df_yfsegg_merged1[, x])/6))

windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfsegg_final1, "Forecasted Distribution 2015 - 2039")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_egg_avg1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

##### Larvae
df_yfslarvae_avg1_cesm126 <- readRDS(here('data', 'df_yfslarvae_avg1_cesm126.rds'))
df_yfslarvae_avg4_cesm585 <- readRDS(here('data', 'df_yfslarvae_avg4_cesm585.rds'))
df_yfslarvae_avg1_gfdl126 <- readRDS(here('data', 'df_yfslarvae_avg1_gfdl126.rds'))
df_yfslarvae_avg4_gfdl585 <- readRDS(here('data', 'df_yfslarvae_avg4_gfdl585.rds'))
df_yfslarvae_avg1_miroc126 <- readRDS(here('data', 'df_yfslarvae_avg1_miroc126.rds'))
df_yfslarvae_avg4_miroc585 <- readRDS(here('data', 'df_yfslarvae_avg4_miroc585.rds'))

df_yfslarvae_merged1 <- list(df_yfslarvae_avg1_cesm126, df_yfslarvae_avg4_cesm585,
                             df_yfslarvae_avg1_gfdl126, df_yfslarvae_avg4_gfdl585,
                             df_yfslarvae_avg1_miroc126, df_yfslarvae_avg4_miroc585) %>%
  reduce(inner_join, by = c("lon", "lat", "dist"))

x <- grepl("pred", names(df_yfslarvae_merged1), fixed = T)
df_yfslarvae_final1 <- data.frame(lat = df_yfslarvae_merged1$lat,
                                  lon = df_yfslarvae_merged1$lon,
                                  avg_pred = (rowSums(df_yfslarvae_merged1[, x])/6))

windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae_final1, "Forecasted Distribution 2015 - 2039")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_avg1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


#### 2040-2069 ---------------------------------------------------------------------------------------------------------------------
##### Eggs
df_yfsegg_avg2_cesm126 <- readRDS(here('data', 'df_yfsegg_avg2_cesm126.rds'))
df_yfsegg_avg5_cesm585 <- readRDS(here('data', 'df_yfsegg_avg5_cesm585.rds'))
df_yfsegg_avg2_gfdl126 <- readRDS(here('data', 'df_yfsegg_avg2_gfdl126.rds'))
df_yfsegg_avg5_gfdl585 <- readRDS(here('data', 'df_yfsegg_avg5_gfdl585.rds'))
df_yfsegg_avg2_miroc126 <- readRDS(here('data', 'df_yfsegg_avg2_miroc126.rds'))
df_yfsegg_avg5_miroc585 <- readRDS(here('data', 'df_yfsegg_avg5_miroc585.rds'))

df_yfsegg_merged2 <- list(df_yfsegg_avg2_cesm126, df_yfsegg_avg5_cesm585,
                          df_yfsegg_avg2_gfdl126, df_yfsegg_avg5_gfdl585,
                          df_yfsegg_avg2_miroc126, df_yfsegg_avg5_miroc585) %>%
  reduce(inner_join, by = c("lon", "lat", "dist"))

x <- grepl("pred", names(df_yfsegg_merged2), fixed = T)
df_yfsegg_final2 <- data.frame(lat = df_yfsegg_merged2$lat,
                               lon = df_yfsegg_merged2$lon,
                               avg_pred = (rowSums(df_yfsegg_merged2[, x])/6))

windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfsegg_final2, "Forecasted Distribution 2040 - 2069")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_egg_avg2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

##### Larvae
df_yfslarvae_avg2_cesm126 <- readRDS(here('data', 'df_yfslarvae_avg2_cesm126.rds'))
df_yfslarvae_avg5_cesm585 <- readRDS(here('data', 'df_yfslarvae_avg5_cesm585.rds'))
df_yfslarvae_avg2_gfdl126 <- readRDS(here('data', 'df_yfslarvae_avg2_gfdl126.rds'))
df_yfslarvae_avg5_gfdl585 <- readRDS(here('data', 'df_yfslarvae_avg5_gfdl585.rds'))
df_yfslarvae_avg2_miroc126 <- readRDS(here('data', 'df_yfslarvae_avg2_miroc126.rds'))
df_yfslarvae_avg5_miroc585 <- readRDS(here('data', 'df_yfslarvae_avg5_miroc585.rds'))

df_yfslarvae_merged2 <- list(df_yfslarvae_avg2_cesm126, df_yfslarvae_avg5_cesm585,
                             df_yfslarvae_avg2_gfdl126, df_yfslarvae_avg5_gfdl585,
                             df_yfslarvae_avg2_miroc126, df_yfslarvae_avg5_miroc585) %>%
  reduce(inner_join, by = c("lon", "lat", "dist"))

x <- grepl("pred", names(df_yfslarvae_merged2), fixed = T)
df_yfslarvae_final2 <- data.frame(lat = df_yfslarvae_merged2$lat,
                                  lon = df_yfslarvae_merged2$lon,
                                  avg_pred = (rowSums(df_yfslarvae_merged2[, x])/6))

windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae_final2, "Forecasted Distribution 2040 - 2069")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_avg2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


#### 2070-2099----------------------------------------------------------------------------------------------------------------------
df_yfsegg_avg3_cesm126 <- readRDS(here('data', 'df_yfsegg_avg3_cesm126.rds'))
df_yfsegg_avg6_cesm585 <- readRDS(here('data', 'df_yfsegg_avg6_cesm585.rds'))
df_yfsegg_avg3_gfdl126 <- readRDS(here('data', 'df_yfsegg_avg3_gfdl126.rds'))
df_yfsegg_avg6_gfdl585 <- readRDS(here('data', 'df_yfsegg_avg6_gfdl585.rds'))
df_yfsegg_avg3_miroc126 <- readRDS(here('data', 'df_yfsegg_avg3_miroc126.rds'))
df_yfsegg_avg6_miroc585 <- readRDS(here('data', 'df_yfsegg_avg6_miroc585.rds'))

df_yfsegg_merged3 <- list(df_yfsegg_avg3_cesm126, df_yfsegg_avg6_cesm585,
                          df_yfsegg_avg3_gfdl126, df_yfsegg_avg6_gfdl585,
                          df_yfsegg_avg3_miroc126, df_yfsegg_avg6_miroc585) %>%
  reduce(inner_join, by = c("lon", "lat", "dist"))

x <- grepl("pred", names(df_yfsegg_merged3), fixed = T)
df_yfsegg_final3 <- data.frame(lat = df_yfsegg_merged3$lat,
                               lon = df_yfsegg_merged3$lon,
                               avg_pred = (rowSums(df_yfsegg_merged3[, x])/6))

windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfsegg_final3, "Forecasted Distribution 2070 - 2099")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_egg_avg3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

df_yfslarvae_avg3_cesm126 <- readRDS(here('data', 'df_yfslarvae_avg3_cesm126.rds'))
df_yfslarvae_avg6_cesm585 <- readRDS(here('data', 'df_yfslarvae_avg6_cesm585.rds'))
df_yfslarvae_avg3_gfdl126 <- readRDS(here('data', 'df_yfslarvae_avg3_gfdl126.rds'))
df_yfslarvae_avg6_gfdl585 <- readRDS(here('data', 'df_yfslarvae_avg6_gfdl585.rds'))
df_yfslarvae_avg3_miroc126 <- readRDS(here('data', 'df_yfslarvae_avg3_miroc126.rds'))
df_yfslarvae_avg6_miroc585 <- readRDS(here('data', 'df_yfslarvae_avg6_miroc585.rds'))

df_yfslarvae_merged3 <- list(df_yfslarvae_avg3_cesm126, df_yfslarvae_avg6_cesm585,
                             df_yfslarvae_avg3_gfdl126, df_yfslarvae_avg6_gfdl585,
                             df_yfslarvae_avg3_miroc126, df_yfslarvae_avg6_miroc585) %>%
  reduce(inner_join, by = c("lon", "lat", "dist"))

x <- grepl("pred", names(df_yfslarvae_merged3), fixed = T)
df_yfslarvae_final3 <- data.frame(lat = df_yfslarvae_merged3$lat,
                                  lon = df_yfslarvae_merged3$lon,
                                  avg_pred = (rowSums(df_yfslarvae_merged3[, x])/6))

windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae_final3, "Forecasted Distribution 2070 - 2099")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_avg3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()