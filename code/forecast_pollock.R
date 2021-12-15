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
pk_egg <- as.data.frame(filter((readRDS(here('data', 'pk_egg.rds'))),
                               lat >= 52 & lat <= 62,
                               lon >= -176.5 & lon <= -156.5))
pk_egg$mean_temp <- roms_temps$mean[match(pk_egg$year, roms_temps$year)]

pk_larvae <- as.data.frame(filter(readRDS(here('data', 'pk_larvae.rds')),
                                  lat >= 52 & lat <= 62,
                                  lon >= -176.5 & lon <= -156.5))
pk_larvae$mean_temp <- roms_temps$mean[match(pk_larvae$year, roms_temps$year)]




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
        zlim = c(-13.6, 8.5),
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
             zlim = c(-13.6, 8.5),
             legend.args = list("Avg. Predicted \n Occurrence",
                                side = 2, cex = 1))
}


# Function to get predictions
get_preds <- function(data, year, date, doy, 
                      start_date, end_date,
                      temp_output, salt_output,
                      list){
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
  
  # Select out box on the shelf
  temp_data <- as.data.frame(cbind(lon = as.vector(temp_output[[list]][[1]]), 
                                   lat = as.vector(temp_output[[list]][[2]]), 
                                   temp = as.vector(temp_array)))
  temp_filtered <- temp_data %>% filter(lon >= -170 & lon <= -165, lat >= 56 & lat <= 58)
  mean <- mean(temp_filtered$temp, na.rm = T)
  
  grid_extent$mean_temp <- mean
  
  gam <- gam(larvalcatchper10m2 + 1 ~ s(year) +
               s(doy, k = 8) +
               s(lon, lat) +
               s(roms_temperature, k = 6) +
               s(roms_salinity, k = 6) +
               s(lat, lon, by = mean_temp, k = 6),
             data = data,
             family = tw(link = 'log'),
             method = 'REML')
  
  # Predict on forecasted output
  grid_extent$pred <- predict(gam, 
                              newdata = grid_extent,
                              type = "link")
  grid_extent$pred[grid_extent$dist > 30000] <- NA
  
  return(grid_extent)
  
}

# Function to loop through years
pred_loop <- function(range, data, doy, 
                      temp_output, salt_output,
                      list){
  grids_pkegg <- list()
  for(j in range) {
    date1 <- paste(j, "-05-10", sep = "")
    date2 <- paste(j, "-02-01", sep = "")
    date3 <- paste(j, "-04-30", sep = "")
    grid <- get_preds(data, j, date1, doy,
                      date2, date3,
                      temp_output, salt_output,
                      list)
    grids_pkegg[[paste("year", j, sep = "")]] <- grid
  }
  return(grids_pkegg)
}


### Pollock Eggs --------------------------------------------------------------------------------------------------------------------------
#### Forecast and average into 3 time periods ---------------------------------------------------------------------------------------------
##### CESM 126 ----------------------------------------------------------------------------------------------------------------------------
temps_cesm_ssp126 <- readRDS(here('data', 'temps_cesm_ssp126.rds'))
salts_cesm_ssp126 <- readRDS(here('data', 'salts_cesm_ssp126.rds'))

## 2015 - 2039
grids_pkegg1 <- pred_loop(2015:2019, pk_egg, 130, 
                          temps_cesm_ssp126,
                          salts_cesm_ssp126, 1)
grids_pkegg2 <- pred_loop(2020:2024, pk_egg, 130, 
                          temps_cesm_ssp126,
                          salts_cesm_ssp126, 2)
grids_pkegg3 <- pred_loop(2025:2029, pk_egg, 130,
                          temps_cesm_ssp126,
                          salts_cesm_ssp126, 3)
grids_pkegg4 <- pred_loop(2030:2034, pk_egg, 130, 
                          temps_cesm_ssp126,
                          salts_cesm_ssp126, 4)
grids_pkegg5 <- pred_loop(2035:2039, pk_egg, 130, 
                          temps_cesm_ssp126,
                          salts_cesm_ssp126, 5)

# Combine into one data frame
df_pkegg1 <- list(grids_pkegg1[[1]], grids_pkegg1[[2]], grids_pkegg1[[3]], 
                  grids_pkegg1[[4]], grids_pkegg1[[5]], grids_pkegg2[[1]], 
                  grids_pkegg2[[2]], grids_pkegg2[[3]], grids_pkegg2[[4]],
                  grids_pkegg2[[5]], grids_pkegg3[[1]], grids_pkegg3[[2]], 
                  grids_pkegg3[[3]], grids_pkegg3[[4]], grids_pkegg3[[5]],
                  grids_pkegg4[[1]], grids_pkegg4[[2]], grids_pkegg4[[3]], 
                  grids_pkegg4[[4]], grids_pkegg4[[5]], grids_pkegg5[[1]], 
                  grids_pkegg5[[2]], grids_pkegg5[[3]], grids_pkegg5[[4]],
                  grids_pkegg5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_pkegg1), fixed = T)
df_pkegg_avg1 <- data.frame(lat = df_pkegg1$lat, 
                           lon = df_pkegg1$lon, 
                           dist = df_pkegg1$dist,
                           avg_pred = rowSums(df_pkegg1[, x])/25)
saveRDS(df_pkegg_avg1, file = here("data", "df_pkegg_avg1.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pkegg_avg1, "Forecasted Distribution 2015 - 2039 \n CESM SSP126")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_egg_cesm_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_pkegg6 <- pred_loop(2040:2044, pk_egg, 130, 
                          temps_cesm_ssp126,
                          salts_cesm_ssp126, 6)
grids_pkegg7 <- pred_loop(2045:2049, pk_egg, 130, 
                          temps_cesm_ssp126,
                          salts_cesm_ssp126, 7)
grids_pkegg8 <- pred_loop(2050:2054, pk_egg, 130, 
                          temps_cesm_ssp126,
                          salts_cesm_ssp126, 8)
grids_pkegg9 <- pred_loop(2055:2059, pk_egg, 130, 
                          temps_cesm_ssp126,
                          salts_cesm_ssp126, 9)
grids_pkegg10 <- pred_loop(2060:2064, pk_egg, 130, 
                          temps_cesm_ssp126,
                          salts_cesm_ssp126, 10)
grids_pkegg11 <- pred_loop(2065:2069, pk_egg, 130, 
                           temps_cesm_ssp126,
                           salts_cesm_ssp126, 11)

# Combine into one data frame
df_pkegg2 <- list(grids_pkegg6[[1]], grids_pkegg6[[2]], grids_pkegg6[[3]], 
                  grids_pkegg6[[4]], grids_pkegg6[[5]], grids_pkegg7[[1]], 
                  grids_pkegg7[[2]], grids_pkegg7[[3]], grids_pkegg7[[4]],
                  grids_pkegg7[[5]], grids_pkegg8[[1]], grids_pkegg8[[2]], 
                  grids_pkegg8[[3]], grids_pkegg8[[4]], grids_pkegg8[[5]],
                  grids_pkegg9[[1]], grids_pkegg9[[2]], grids_pkegg9[[3]], 
                  grids_pkegg9[[4]], grids_pkegg9[[5]], grids_pkegg10[[1]], 
                  grids_pkegg10[[2]], grids_pkegg10[[3]], grids_pkegg10[[4]],
                  grids_pkegg11[[5]], grids_pkegg11[[1]], grids_pkegg11[[2]],
                  grids_pkegg11[[3]], grids_pkegg11[[4]], grids_pkegg11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_pkegg2), fixed = T)
df_pkegg_avg2 <- data.frame(lat = df_pkegg2$lat, 
                           lon = df_pkegg2$lon, 
                           dist = df_pkegg2$dist,
                           avg_pred = rowSums(df_pkegg2[, x])/30)
saveRDS(df_pkegg_avg2, file = here("data", "df_pkegg2_avg.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pkegg_avg2, "Forecasted Distribution 2040 - 2069 \n CESM SSP126")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_egg_cesm_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_pkegg12 <- pred_loop(2070:2074, pk_egg, 130,
                           temps_cesm_ssp126,
                           salts_cesm_ssp126, 12)
grids_pkegg13 <- pred_loop(2075:2079, pk_egg, 130,
                           temps_cesm_ssp126,
                           salts_cesm_ssp126, 13)
grids_pkegg14 <- pred_loop(2080:2084, pk_egg, 130,
                           temps_cesm_ssp126,
                           salts_cesm_ssp126, 14)
grids_pkegg15 <- pred_loop(2085:2089, pk_egg, 130,
                           temps_cesm_ssp126,
                           salts_cesm_ssp126, 15)
grids_pkegg16 <- pred_loop(2090:2094, pk_egg, 130,
                           temps_cesm_ssp126,
                           salts_cesm_ssp126, 16)
grids_pkegg17 <- pred_loop(2095:2099, pk_egg, 130, 
                           temps_cesm_ssp126,
                           salts_cesm_ssp126, 17)

# Combine into one data frame
df_pkegg3 <- list(grids_pkegg12[[1]], grids_pkegg12[[2]], grids_pkegg12[[3]], 
                  grids_pkegg12[[4]], grids_pkegg12[[5]], grids_pkegg13[[1]], 
                  grids_pkegg13[[2]], grids_pkegg13[[3]], grids_pkegg13[[4]],
                  grids_pkegg13[[5]], grids_pkegg14[[1]], grids_pkegg14[[2]], 
                  grids_pkegg14[[3]], grids_pkegg14[[4]], grids_pkegg14[[5]],
                  grids_pkegg15[[1]], grids_pkegg15[[2]], grids_pkegg15[[3]], 
                  grids_pkegg15[[4]], grids_pkegg15[[5]], grids_pkegg16[[1]], 
                  grids_pkegg16[[2]], grids_pkegg16[[3]], grids_pkegg16[[4]],
                  grids_pkegg17[[5]], grids_pkegg17[[1]], grids_pkegg17[[2]],
                  grids_pkegg17[[3]], grids_pkegg17[[4]], grids_pkegg17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_pkegg3), fixed = T)
df_pkegg_avg3 <- data.frame(lat = df_pkegg3$lat, 
                            lon = df_pkegg3$lon, 
                            dist = df_pkegg3$dist,
                            avg_pred = rowSums(df_pkegg3[, x])/30)
saveRDS(df_pkegg_avg3, file = here("data", "df_pkegg3_avg.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pkegg_avg3, "Forecasted Distribution 2070 - 2099 \n CESM SSP126")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_egg_cesm_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

##### CESM 585 -----------------------------------------------------------------------------------------------------------------
temps_cesm_ssp585 <- readRDS(here('data', 'temps_cesm_ssp585.rds'))
salts_cesm_ssp585 <- readRDS(here('data', 'salts_cesm_ssp585.rds'))

## 2015 - 2039
grids_pkegg1 <- pred_loop(2015:2019, pk_egg, 130, 
                          temps_cesm_ssp585,
                          salts_cesm_ssp585, 1)
grids_pkegg2 <- pred_loop(2020:2024, pk_egg, 130, 
                          temps_cesm_ssp585,
                          salts_cesm_ssp585, 2)
grids_pkegg3 <- pred_loop(2025:2029, pk_egg, 130,
                          temps_cesm_ssp585,
                          salts_cesm_ssp585, 3)
grids_pkegg4 <- pred_loop(2030:2034, pk_egg, 130, 
                          temps_cesm_ssp585,
                          salts_cesm_ssp585, 4)
grids_pkegg5 <- pred_loop(2035:2039, pk_egg, 130, 
                          temps_cesm_ssp585,
                          salts_cesm_ssp585, 5)

# Combine into one data frame
df_pkegg4 <- list(grids_pkegg1[[1]], grids_pkegg1[[2]], grids_pkegg1[[3]], 
                  grids_pkegg1[[4]], grids_pkegg1[[5]], grids_pkegg2[[1]], 
                  grids_pkegg2[[2]], grids_pkegg2[[3]], grids_pkegg2[[4]],
                  grids_pkegg2[[5]], grids_pkegg3[[1]], grids_pkegg3[[2]], 
                  grids_pkegg3[[3]], grids_pkegg3[[4]], grids_pkegg3[[5]],
                  grids_pkegg4[[1]], grids_pkegg4[[2]], grids_pkegg4[[3]], 
                  grids_pkegg4[[4]], grids_pkegg4[[5]], grids_pkegg5[[1]], 
                  grids_pkegg5[[2]], grids_pkegg5[[3]], grids_pkegg5[[4]],
                  grids_pkegg5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_pkegg4), fixed = T)
df_pkegg_avg4 <- data.frame(lat = df_pkegg4$lat, 
                            lon = df_pkegg4$lon, 
                            dist = df_pkegg4$dist,
                            avg_pred = rowSums(df_pkegg4[, x])/25)
saveRDS(df_pkegg_avg4, file = here("data", "df_pkegg_avg4.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pkegg_avg4, "Forecasted Distribution 2015 - 2039 \n CESM SSP585")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_egg_cesm_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_pkegg6 <- pred_loop(2040:2044, pk_egg, 130, 
                          temps_cesm_ssp585,
                          salts_cesm_ssp585, 6)
grids_pkegg7 <- pred_loop(2045:2049, pk_egg, 130, 
                          temps_cesm_ssp585,
                          salts_cesm_ssp585, 7)
grids_pkegg8 <- pred_loop(2050:2054, pk_egg, 130, 
                          temps_cesm_ssp585,
                          salts_cesm_ssp585, 8)
grids_pkegg9 <- pred_loop(2055:2059, pk_egg, 130, 
                          temps_cesm_ssp585,
                          salts_cesm_ssp585, 9)
grids_pkegg10 <- pred_loop(2060:2064, pk_egg, 130, 
                           temps_cesm_ssp585,
                           salts_cesm_ssp585, 10)
grids_pkegg11 <- pred_loop(2065:2069, pk_egg, 130, 
                           temps_cesm_ssp585,
                           salts_cesm_ssp585, 11)

# Combine into one data frame
df_pkegg5 <- list(grids_pkegg6[[1]], grids_pkegg6[[2]], grids_pkegg6[[3]], 
                  grids_pkegg6[[4]], grids_pkegg6[[5]], grids_pkegg7[[1]], 
                  grids_pkegg7[[2]], grids_pkegg7[[3]], grids_pkegg7[[4]],
                  grids_pkegg7[[5]], grids_pkegg8[[1]], grids_pkegg8[[2]], 
                  grids_pkegg8[[3]], grids_pkegg8[[4]], grids_pkegg8[[5]],
                  grids_pkegg9[[1]], grids_pkegg9[[2]], grids_pkegg9[[3]], 
                  grids_pkegg9[[4]], grids_pkegg9[[5]], grids_pkegg10[[1]], 
                  grids_pkegg10[[2]], grids_pkegg10[[3]], grids_pkegg10[[4]],
                  grids_pkegg11[[5]], grids_pkegg11[[1]], grids_pkegg11[[2]],
                  grids_pkegg11[[3]], grids_pkegg11[[4]], grids_pkegg11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_pkegg5), fixed = T)
df_pkegg_avg5 <- data.frame(lat = df_pkegg5$lat, 
                            lon = df_pkegg5$lon, 
                            dist = df_pkegg5$dist,
                            avg_pred = rowSums(df_pkegg5[, x])/30)
saveRDS(df_pkegg_avg5, file = here("data", "df_pkegg5_avg.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pkegg_avg5, "Forecasted Distribution 5040 - 5069 \n CESM SSP585")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_egg_cesm_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_pkegg12 <- pred_loop(2070:2074, pk_egg, 130,
                           temps_cesm_ssp585,
                           salts_cesm_ssp585, 12)
grids_pkegg13 <- pred_loop(2075:2079, pk_egg, 130,
                           temps_cesm_ssp585,
                           salts_cesm_ssp585, 13)
grids_pkegg14 <- pred_loop(2080:2084, pk_egg, 130,
                           temps_cesm_ssp585,
                           salts_cesm_ssp585, 14)
grids_pkegg15 <- pred_loop(2085:2089, pk_egg, 130,
                           temps_cesm_ssp585,
                           salts_cesm_ssp585, 15)
grids_pkegg16 <- pred_loop(2090:2094, pk_egg, 130,
                           temps_cesm_ssp585,
                           salts_cesm_ssp585, 16)
grids_pkegg17 <- pred_loop(2095:2099, pk_egg, 130, 
                           temps_cesm_ssp585,
                           salts_cesm_ssp585, 17)

# Combine into one data frame
df_pkegg6 <- list(grids_pkegg12[[1]], grids_pkegg12[[2]], grids_pkegg12[[3]], 
                  grids_pkegg12[[4]], grids_pkegg12[[5]], grids_pkegg13[[1]], 
                  grids_pkegg13[[2]], grids_pkegg13[[3]], grids_pkegg13[[4]],
                  grids_pkegg13[[5]], grids_pkegg14[[1]], grids_pkegg14[[2]], 
                  grids_pkegg14[[3]], grids_pkegg14[[4]], grids_pkegg14[[5]],
                  grids_pkegg15[[1]], grids_pkegg15[[2]], grids_pkegg15[[3]], 
                  grids_pkegg15[[4]], grids_pkegg15[[5]], grids_pkegg16[[1]], 
                  grids_pkegg16[[2]], grids_pkegg16[[3]], grids_pkegg16[[4]],
                  grids_pkegg17[[5]], grids_pkegg17[[1]], grids_pkegg17[[2]],
                  grids_pkegg17[[3]], grids_pkegg17[[4]], grids_pkegg17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_pkegg6), fixed = T)
df_pkegg_avg6 <- data.frame(lat = df_pkegg6$lat, 
                            lon = df_pkegg6$lon, 
                            dist = df_pkegg6$dist,
                            avg_pred = rowSums(df_pkegg6[, x])/60)
saveRDS(df_pkegg_avg6, file = here("data", "df_pkegg6_avg.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pkegg_avg6, "Forecasted Distribution 2070 - 2099 \n CESM SSP585")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_egg_cesm_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()