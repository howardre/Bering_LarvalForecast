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
fhs_egg <- as.data.frame(filter((readRDS(here('data', 'fhs_egg.rds'))),
                               lat >= 52 & lat <= 62,
                               lon >= -176.5 & lon <= -156.5))
fhs_egg$mean_temp <- roms_temps$mean[match(fhs_egg$year, roms_temps$year)]
fhs_egg$catch <- fhs_egg$larvalcatchper10m2 + 1

fhs_larvae <- as.data.frame(filter(readRDS(here('data', 'fhs_larvae.rds')),
                                  lat >= 52 & lat <= 62,
                                  lon >= -176.5 & lon <= -156.5))
fhs_larvae$mean_temp <- roms_temps$mean[match(fhs_larvae$year, roms_temps$year)]
fhs_larvae$catch <- fhs_larvae$larvalcatchper10m2 + 1

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
                     s(roms_salinity, k = 6) +
                     s(doy, by = mean_temp, k = 6),
                   data = fhs_egg,
                   family = tw(link = 'log'),
                   method = 'REML')

larval_formula <- gam(catch ~ s(year) +
                        s(doy, k = 8) +
                        s(lon, lat) +
                        s(roms_temperature, k = 6) +
                        s(roms_salinity, k = 6) +
                        s(lat, lon, by = mean_temp, k = 6),
                      data = fhs_larvae,
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
  grids_fhsegg <- list()
  for(j in range) {
    date1 <- paste(j, "-05-10", sep = "")
    date2 <- paste(j, "-02-01", sep = "")
    date3 <- paste(j, "-04-30", sep = "")
    grid <- get_preds(data, j, date1, doy,
                      date2, date3,
                      temp_output, salt_output,
                      list, formula)
    grids_fhsegg[[paste("year", j, sep = "")]] <- grid
  }
  return(grids_fhsegg)
}


### Flathead Eggs --------------------------------------------------------------------------------------------------------------------------
#### Forecast and average into 3 time periods ---------------------------------------------------------------------------------------------
##### CESM 126 ----------------------------------------------------------------------------------------------------------------------------
temps_cesm_ssp126 <- readRDS(here('data', 'temps_cesm_ssp126.rds'))
salts_cesm_ssp126 <- readRDS(here('data', 'salts_cesm_ssp126.rds'))

## 2015 - 2039
grids_fhsegg1 <- pred_loop(2015:2019, fhs_egg, 130, 
                          temps_cesm_ssp126,
                          salts_cesm_ssp126, 1,
                          egg_formula)
grids_fhsegg2 <- pred_loop(2020:2024, fhs_egg, 130, 
                          temps_cesm_ssp126,
                          salts_cesm_ssp126, 2,
                          egg_formula)
grids_fhsegg3 <- pred_loop(2025:2029, fhs_egg, 130,
                          temps_cesm_ssp126,
                          salts_cesm_ssp126, 3,
                          egg_formula)
grids_fhsegg4 <- pred_loop(2030:2034, fhs_egg, 130, 
                          temps_cesm_ssp126,
                          salts_cesm_ssp126, 4,
                          egg_formula)
grids_fhsegg5 <- pred_loop(2035:2039, fhs_egg, 130, 
                          temps_cesm_ssp126,
                          salts_cesm_ssp126, 5,
                          egg_formula)

# Combine into one data frame
df_fhsegg1 <- list(grids_fhsegg1[[1]], grids_fhsegg1[[2]], grids_fhsegg1[[3]], 
                  grids_fhsegg1[[4]], grids_fhsegg1[[5]], grids_fhsegg2[[1]], 
                  grids_fhsegg2[[2]], grids_fhsegg2[[3]], grids_fhsegg2[[4]],
                  grids_fhsegg2[[5]], grids_fhsegg3[[1]], grids_fhsegg3[[2]], 
                  grids_fhsegg3[[3]], grids_fhsegg3[[4]], grids_fhsegg3[[5]],
                  grids_fhsegg4[[1]], grids_fhsegg4[[2]], grids_fhsegg4[[3]], 
                  grids_fhsegg4[[4]], grids_fhsegg4[[5]], grids_fhsegg5[[1]], 
                  grids_fhsegg5[[2]], grids_fhsegg5[[3]], grids_fhsegg5[[4]],
                  grids_fhsegg5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhsegg1), fixed = T)
df_fhsegg_avg1_cesm126 <- data.frame(lat = df_fhsegg1$lat, 
                                    lon = df_fhsegg1$lon, 
                                    dist = df_fhsegg1$dist,
                                    avg_pred = rowSums(df_fhsegg1[, x])/25)
saveRDS(df_fhsegg_avg1_cesm126, file = here("data", "df_fhsegg_avg1_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg_avg1_cesm126, "Forecasted Distribution 2015 - 2039 \n CESM SSP126")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_cesm_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_fhsegg6 <- pred_loop(2040:2044, fhs_egg, 130, 
                          temps_cesm_ssp126,
                          salts_cesm_ssp126, 6,
                          egg_formula)
grids_fhsegg7 <- pred_loop(2045:2049, fhs_egg, 130, 
                          temps_cesm_ssp126,
                          salts_cesm_ssp126, 7,
                          egg_formula)
grids_fhsegg8 <- pred_loop(2050:2054, fhs_egg, 130, 
                          temps_cesm_ssp126,
                          salts_cesm_ssp126, 8,
                          egg_formula)
grids_fhsegg9 <- pred_loop(2055:2059, fhs_egg, 130, 
                          temps_cesm_ssp126,
                          salts_cesm_ssp126, 9,
                          egg_formula)
grids_fhsegg10 <- pred_loop(2060:2064, fhs_egg, 130, 
                           temps_cesm_ssp126,
                           salts_cesm_ssp126, 10,
                           egg_formula)
grids_fhsegg11 <- pred_loop(2065:2069, fhs_egg, 130, 
                           temps_cesm_ssp126,
                           salts_cesm_ssp126, 11,
                           egg_formula)

# Combine into one data frame
df_fhsegg2 <- list(grids_fhsegg6[[1]], grids_fhsegg6[[2]], grids_fhsegg6[[3]], 
                  grids_fhsegg6[[4]], grids_fhsegg6[[5]], grids_fhsegg7[[1]], 
                  grids_fhsegg7[[2]], grids_fhsegg7[[3]], grids_fhsegg7[[4]],
                  grids_fhsegg7[[5]], grids_fhsegg8[[1]], grids_fhsegg8[[2]], 
                  grids_fhsegg8[[3]], grids_fhsegg8[[4]], grids_fhsegg8[[5]],
                  grids_fhsegg9[[1]], grids_fhsegg9[[2]], grids_fhsegg9[[3]], 
                  grids_fhsegg9[[4]], grids_fhsegg9[[5]], grids_fhsegg10[[1]], 
                  grids_fhsegg10[[2]], grids_fhsegg10[[3]], grids_fhsegg10[[4]],
                  grids_fhsegg11[[5]], grids_fhsegg11[[1]], grids_fhsegg11[[2]],
                  grids_fhsegg11[[3]], grids_fhsegg11[[4]], grids_fhsegg11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhsegg2), fixed = T)
df_fhsegg_avg2_cesm126 <- data.frame(lat = df_fhsegg2$lat, 
                                    lon = df_fhsegg2$lon, 
                                    dist = df_fhsegg2$dist,
                                    avg_pred = rowSums(df_fhsegg2[, x])/30)
saveRDS(df_fhsegg_avg2_cesm126, file = here("data", "df_fhsegg_avg2_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg_avg2_cesm126, "Forecasted Distribution 2040 - 2069 \n CESM SSP126")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_cesm_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_fhsegg12 <- pred_loop(2070:2074, fhs_egg, 130,
                           temps_cesm_ssp126,
                           salts_cesm_ssp126, 12,
                           egg_formula)
grids_fhsegg13 <- pred_loop(2075:2079, fhs_egg, 130,
                           temps_cesm_ssp126,
                           salts_cesm_ssp126, 13,
                           egg_formula)
grids_fhsegg14 <- pred_loop(2080:2084, fhs_egg, 130,
                           temps_cesm_ssp126,
                           salts_cesm_ssp126, 14,
                           egg_formula)
grids_fhsegg15 <- pred_loop(2085:2089, fhs_egg, 130,
                           temps_cesm_ssp126,
                           salts_cesm_ssp126, 15,
                           egg_formula)
grids_fhsegg16 <- pred_loop(2090:2094, fhs_egg, 130,
                           temps_cesm_ssp126,
                           salts_cesm_ssp126, 16,
                           egg_formula)
grids_fhsegg17 <- pred_loop(2095:2099, fhs_egg, 130, 
                           temps_cesm_ssp126,
                           salts_cesm_ssp126, 17,
                           egg_formula)

# Combine into one data frame
df_fhsegg3 <- list(grids_fhsegg12[[1]], grids_fhsegg12[[2]], grids_fhsegg12[[3]], 
                  grids_fhsegg12[[4]], grids_fhsegg12[[5]], grids_fhsegg13[[1]], 
                  grids_fhsegg13[[2]], grids_fhsegg13[[3]], grids_fhsegg13[[4]],
                  grids_fhsegg13[[5]], grids_fhsegg14[[1]], grids_fhsegg14[[2]], 
                  grids_fhsegg14[[3]], grids_fhsegg14[[4]], grids_fhsegg14[[5]],
                  grids_fhsegg15[[1]], grids_fhsegg15[[2]], grids_fhsegg15[[3]], 
                  grids_fhsegg15[[4]], grids_fhsegg15[[5]], grids_fhsegg16[[1]], 
                  grids_fhsegg16[[2]], grids_fhsegg16[[3]], grids_fhsegg16[[4]],
                  grids_fhsegg17[[5]], grids_fhsegg17[[1]], grids_fhsegg17[[2]],
                  grids_fhsegg17[[3]], grids_fhsegg17[[4]], grids_fhsegg17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhsegg3), fixed = T)
df_fhsegg_avg3_cesm126 <- data.frame(lat = df_fhsegg3$lat, 
                                    lon = df_fhsegg3$lon, 
                                    dist = df_fhsegg3$dist,
                                    avg_pred = rowSums(df_fhsegg3[, x])/30)
saveRDS(df_fhsegg_avg3_cesm126, file = here("data", "df_fhsegg_avg3_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg_avg3_cesm126, "Forecasted Distribution 2070 - 2099 \n CESM SSP126")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_cesm_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

##### CESM 585 -----------------------------------------------------------------------------------------------------------------
temps_cesm_ssp585 <- readRDS(here('data', 'temps_cesm_ssp585.rds'))
salts_cesm_ssp585 <- readRDS(here('data', 'salts_cesm_ssp585.rds'))

## 2015 - 2039
grids_fhsegg1 <- pred_loop(2015:2019, fhs_egg, 130, 
                          temps_cesm_ssp585,
                          salts_cesm_ssp585, 1,
                          egg_formula)
grids_fhsegg2 <- pred_loop(2020:2024, fhs_egg, 130, 
                          temps_cesm_ssp585,
                          salts_cesm_ssp585, 2,
                          egg_formula)
grids_fhsegg3 <- pred_loop(2025:2029, fhs_egg, 130,
                          temps_cesm_ssp585,
                          salts_cesm_ssp585, 3,
                          egg_formula)
grids_fhsegg4 <- pred_loop(2030:2034, fhs_egg, 130, 
                          temps_cesm_ssp585,
                          salts_cesm_ssp585, 4,
                          egg_formula)
grids_fhsegg5 <- pred_loop(2035:2039, fhs_egg, 130, 
                          temps_cesm_ssp585,
                          salts_cesm_ssp585, 5,
                          egg_formula)

# Combine into one data frame
df_fhsegg4 <- list(grids_fhsegg1[[1]], grids_fhsegg1[[2]], grids_fhsegg1[[3]], 
                  grids_fhsegg1[[4]], grids_fhsegg1[[5]], grids_fhsegg2[[1]], 
                  grids_fhsegg2[[2]], grids_fhsegg2[[3]], grids_fhsegg2[[4]],
                  grids_fhsegg2[[5]], grids_fhsegg3[[1]], grids_fhsegg3[[2]], 
                  grids_fhsegg3[[3]], grids_fhsegg3[[4]], grids_fhsegg3[[5]],
                  grids_fhsegg4[[1]], grids_fhsegg4[[2]], grids_fhsegg4[[3]], 
                  grids_fhsegg4[[4]], grids_fhsegg4[[5]], grids_fhsegg5[[1]], 
                  grids_fhsegg5[[2]], grids_fhsegg5[[3]], grids_fhsegg5[[4]],
                  grids_fhsegg5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhsegg4), fixed = T)
df_fhsegg_avg4_cesm585 <- data.frame(lat = df_fhsegg4$lat, 
                                    lon = df_fhsegg4$lon, 
                                    dist = df_fhsegg4$dist,
                                    avg_pred = rowSums(df_fhsegg4[, x])/25)
saveRDS(df_fhsegg_avg4_cesm585, file = here("data", "df_fhsegg_avg4_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg_avg4_cesm585, "Forecasted Distribution 2015 - 2039 \n CESM SSP585")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_cesm_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_fhsegg6 <- pred_loop(2040:2044, fhs_egg, 130, 
                          temps_cesm_ssp585,
                          salts_cesm_ssp585, 6,
                          egg_formula)
grids_fhsegg7 <- pred_loop(2045:2049, fhs_egg, 130, 
                          temps_cesm_ssp585,
                          salts_cesm_ssp585, 7,
                          egg_formula)
grids_fhsegg8 <- pred_loop(2050:2054, fhs_egg, 130, 
                          temps_cesm_ssp585,
                          salts_cesm_ssp585, 8,
                          egg_formula)
grids_fhsegg9 <- pred_loop(2055:2059, fhs_egg, 130, 
                          temps_cesm_ssp585,
                          salts_cesm_ssp585, 9,
                          egg_formula)
grids_fhsegg10 <- pred_loop(2060:2064, fhs_egg, 130, 
                           temps_cesm_ssp585,
                           salts_cesm_ssp585, 10,
                           egg_formula)
grids_fhsegg11 <- pred_loop(2065:2069, fhs_egg, 130, 
                           temps_cesm_ssp585,
                           salts_cesm_ssp585, 11,
                           egg_formula)

# Combine into one data frame
df_fhsegg5 <- list(grids_fhsegg6[[1]], grids_fhsegg6[[2]], grids_fhsegg6[[3]], 
                  grids_fhsegg6[[4]], grids_fhsegg6[[5]], grids_fhsegg7[[1]], 
                  grids_fhsegg7[[2]], grids_fhsegg7[[3]], grids_fhsegg7[[4]],
                  grids_fhsegg7[[5]], grids_fhsegg8[[1]], grids_fhsegg8[[2]], 
                  grids_fhsegg8[[3]], grids_fhsegg8[[4]], grids_fhsegg8[[5]],
                  grids_fhsegg9[[1]], grids_fhsegg9[[2]], grids_fhsegg9[[3]], 
                  grids_fhsegg9[[4]], grids_fhsegg9[[5]], grids_fhsegg10[[1]], 
                  grids_fhsegg10[[2]], grids_fhsegg10[[3]], grids_fhsegg10[[4]],
                  grids_fhsegg11[[5]], grids_fhsegg11[[1]], grids_fhsegg11[[2]],
                  grids_fhsegg11[[3]], grids_fhsegg11[[4]], grids_fhsegg11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhsegg5), fixed = T)
df_fhsegg_avg5_cesm585 <- data.frame(lat = df_fhsegg5$lat, 
                                    lon = df_fhsegg5$lon, 
                                    dist = df_fhsegg5$dist,
                                    avg_pred = rowSums(df_fhsegg5[, x])/30)
saveRDS(df_fhsegg_avg5_cesm585, file = here("data", "df_fhsegg_avg5_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg_avg5_cesm585, "Forecasted Distribution 2040 - 2069 \n CESM SSP585")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_cesm_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_fhsegg12 <- pred_loop(2070:2074, fhs_egg, 130,
                           temps_cesm_ssp585,
                           salts_cesm_ssp585, 12,
                           egg_formula)
grids_fhsegg13 <- pred_loop(2075:2079, fhs_egg, 130,
                           temps_cesm_ssp585,
                           salts_cesm_ssp585, 13,
                           egg_formula)
grids_fhsegg14 <- pred_loop(2080:2084, fhs_egg, 130,
                           temps_cesm_ssp585,
                           salts_cesm_ssp585, 14,
                           egg_formula)
grids_fhsegg15 <- pred_loop(2085:2089, fhs_egg, 130,
                           temps_cesm_ssp585,
                           salts_cesm_ssp585, 15,
                           egg_formula)
grids_fhsegg16 <- pred_loop(2090:2094, fhs_egg, 130,
                           temps_cesm_ssp585,
                           salts_cesm_ssp585, 16,
                           egg_formula)
grids_fhsegg17 <- pred_loop(2095:2099, fhs_egg, 130, 
                           temps_cesm_ssp585,
                           salts_cesm_ssp585, 17,
                           egg_formula)

# Combine into one data frame
df_fhsegg6 <- list(grids_fhsegg12[[1]], grids_fhsegg12[[2]], grids_fhsegg12[[3]], 
                  grids_fhsegg12[[4]], grids_fhsegg12[[5]], grids_fhsegg13[[1]], 
                  grids_fhsegg13[[2]], grids_fhsegg13[[3]], grids_fhsegg13[[4]],
                  grids_fhsegg13[[5]], grids_fhsegg14[[1]], grids_fhsegg14[[2]], 
                  grids_fhsegg14[[3]], grids_fhsegg14[[4]], grids_fhsegg14[[5]],
                  grids_fhsegg15[[1]], grids_fhsegg15[[2]], grids_fhsegg15[[3]], 
                  grids_fhsegg15[[4]], grids_fhsegg15[[5]], grids_fhsegg16[[1]], 
                  grids_fhsegg16[[2]], grids_fhsegg16[[3]], grids_fhsegg16[[4]],
                  grids_fhsegg17[[5]], grids_fhsegg17[[1]], grids_fhsegg17[[2]],
                  grids_fhsegg17[[3]], grids_fhsegg17[[4]], grids_fhsegg17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhsegg6), fixed = T)
df_fhsegg_avg6_cesm585 <- data.frame(lat = df_fhsegg6$lat, 
                                    lon = df_fhsegg6$lon, 
                                    dist = df_fhsegg6$dist,
                                    avg_pred = rowSums(df_fhsegg6[, x])/30)
saveRDS(df_fhsegg_avg6_cesm585, file = here("data", "df_fhsegg_avg6_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg_avg6_cesm585, "Forecasted Distribution 2070 - 2099 \n CESM SSP585")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_cesm_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### GFDL 126 ----------------------------------------------------------------------------------------------------------------------------
temps_gfdl_ssp126 <- readRDS(here('data', 'temps_gfdl_ssp126.rds'))
salts_gfdl_ssp126 <- readRDS(here('data', 'salts_gfdl_ssp126.rds'))

## 2015 - 2039
grids_fhsegg1 <- pred_loop(2015:2019, fhs_egg, 130, 
                          temps_gfdl_ssp126,
                          salts_gfdl_ssp126, 1,
                          egg_formula)
grids_fhsegg2 <- pred_loop(2020:2024, fhs_egg, 130, 
                          temps_gfdl_ssp126,
                          salts_gfdl_ssp126, 2,
                          egg_formula)
grids_fhsegg3 <- pred_loop(2025:2029, fhs_egg, 130,
                          temps_gfdl_ssp126,
                          salts_gfdl_ssp126, 3,
                          egg_formula)
grids_fhsegg4 <- pred_loop(2030:2034, fhs_egg, 130, 
                          temps_gfdl_ssp126,
                          salts_gfdl_ssp126, 4,
                          egg_formula)
grids_fhsegg5 <- pred_loop(2035:2039, fhs_egg, 130, 
                          temps_gfdl_ssp126,
                          salts_gfdl_ssp126, 5,
                          egg_formula)

# Combine into one data frame
df_fhsegg1 <- list(grids_fhsegg1[[1]], grids_fhsegg1[[2]], grids_fhsegg1[[3]], 
                  grids_fhsegg1[[4]], grids_fhsegg1[[5]], grids_fhsegg2[[1]], 
                  grids_fhsegg2[[2]], grids_fhsegg2[[3]], grids_fhsegg2[[4]],
                  grids_fhsegg2[[5]], grids_fhsegg3[[1]], grids_fhsegg3[[2]], 
                  grids_fhsegg3[[3]], grids_fhsegg3[[4]], grids_fhsegg3[[5]],
                  grids_fhsegg4[[1]], grids_fhsegg4[[2]], grids_fhsegg4[[3]], 
                  grids_fhsegg4[[4]], grids_fhsegg4[[5]], grids_fhsegg5[[1]], 
                  grids_fhsegg5[[2]], grids_fhsegg5[[3]], grids_fhsegg5[[4]],
                  grids_fhsegg5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhsegg1), fixed = T)
df_fhsegg_avg1_gfdl126 <- data.frame(lat = df_fhsegg1$lat, 
                                    lon = df_fhsegg1$lon, 
                                    dist = df_fhsegg1$dist,
                                    avg_pred = rowSums(df_fhsegg1[, x])/25)
saveRDS(df_fhsegg_avg1_gfdl126, file = here("data", "df_fhsegg_avg1_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg_avg1_gfdl126, "Forecasted Distribution 2015 - 2039 \n GFDL SSP126")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_gfdl_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_fhsegg6 <- pred_loop(2040:2044, fhs_egg, 130, 
                          temps_gfdl_ssp126,
                          salts_gfdl_ssp126, 6,
                          egg_formula)
grids_fhsegg7 <- pred_loop(2045:2049, fhs_egg, 130, 
                          temps_gfdl_ssp126,
                          salts_gfdl_ssp126, 7,
                          egg_formula)
grids_fhsegg8 <- pred_loop(2050:2054, fhs_egg, 130, 
                          temps_gfdl_ssp126,
                          salts_gfdl_ssp126, 8,
                          egg_formula)
grids_fhsegg9 <- pred_loop(2055:2059, fhs_egg, 130, 
                          temps_gfdl_ssp126,
                          salts_gfdl_ssp126, 9,
                          egg_formula)
grids_fhsegg10 <- pred_loop(2060:2064, fhs_egg, 130, 
                           temps_gfdl_ssp126,
                           salts_gfdl_ssp126, 10,
                           egg_formula)
grids_fhsegg11 <- pred_loop(2065:2069, fhs_egg, 130, 
                           temps_gfdl_ssp126,
                           salts_gfdl_ssp126, 11,
                           egg_formula)

# Combine into one data frame
df_fhsegg2 <- list(grids_fhsegg6[[1]], grids_fhsegg6[[2]], grids_fhsegg6[[3]], 
                  grids_fhsegg6[[4]], grids_fhsegg6[[5]], grids_fhsegg7[[1]], 
                  grids_fhsegg7[[2]], grids_fhsegg7[[3]], grids_fhsegg7[[4]],
                  grids_fhsegg7[[5]], grids_fhsegg8[[1]], grids_fhsegg8[[2]], 
                  grids_fhsegg8[[3]], grids_fhsegg8[[4]], grids_fhsegg8[[5]],
                  grids_fhsegg9[[1]], grids_fhsegg9[[2]], grids_fhsegg9[[3]], 
                  grids_fhsegg9[[4]], grids_fhsegg9[[5]], grids_fhsegg10[[1]], 
                  grids_fhsegg10[[2]], grids_fhsegg10[[3]], grids_fhsegg10[[4]],
                  grids_fhsegg11[[5]], grids_fhsegg11[[1]], grids_fhsegg11[[2]],
                  grids_fhsegg11[[3]], grids_fhsegg11[[4]], grids_fhsegg11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhsegg2), fixed = T)
df_fhsegg_avg2_gfdl126 <- data.frame(lat = df_fhsegg2$lat, 
                                    lon = df_fhsegg2$lon, 
                                    dist = df_fhsegg2$dist,
                                    avg_pred = rowSums(df_fhsegg2[, x])/30)
saveRDS(df_fhsegg_avg2_gfdl126, file = here("data", "df_fhsegg_avg2_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg_avg2_gfdl126, "Forecasted Distribution 2040 - 2069 \n GFDL SSP126")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_gfdl_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_fhsegg12 <- pred_loop(2070:2074, fhs_egg, 130,
                           temps_gfdl_ssp126,
                           salts_gfdl_ssp126, 12,
                           egg_formula)
grids_fhsegg13 <- pred_loop(2075:2079, fhs_egg, 130,
                           temps_gfdl_ssp126,
                           salts_gfdl_ssp126, 13,
                           egg_formula)
grids_fhsegg14 <- pred_loop(2080:2084, fhs_egg, 130,
                           temps_gfdl_ssp126,
                           salts_gfdl_ssp126, 14,
                           egg_formula)
grids_fhsegg15 <- pred_loop(2085:2089, fhs_egg, 130,
                           temps_gfdl_ssp126,
                           salts_gfdl_ssp126, 15,
                           egg_formula)
grids_fhsegg16 <- pred_loop(2090:2094, fhs_egg, 130,
                           temps_gfdl_ssp126,
                           salts_gfdl_ssp126, 16,
                           egg_formula)
grids_fhsegg17 <- pred_loop(2095:2099, fhs_egg, 130, 
                           temps_gfdl_ssp126,
                           salts_gfdl_ssp126, 17,
                           egg_formula)

# Combine into one data frame
df_fhsegg3 <- list(grids_fhsegg12[[1]], grids_fhsegg12[[2]], grids_fhsegg12[[3]], 
                  grids_fhsegg12[[4]], grids_fhsegg12[[5]], grids_fhsegg13[[1]], 
                  grids_fhsegg13[[2]], grids_fhsegg13[[3]], grids_fhsegg13[[4]],
                  grids_fhsegg13[[5]], grids_fhsegg14[[1]], grids_fhsegg14[[2]], 
                  grids_fhsegg14[[3]], grids_fhsegg14[[4]], grids_fhsegg14[[5]],
                  grids_fhsegg15[[1]], grids_fhsegg15[[2]], grids_fhsegg15[[3]], 
                  grids_fhsegg15[[4]], grids_fhsegg15[[5]], grids_fhsegg16[[1]], 
                  grids_fhsegg16[[2]], grids_fhsegg16[[3]], grids_fhsegg16[[4]],
                  grids_fhsegg17[[5]], grids_fhsegg17[[1]], grids_fhsegg17[[2]],
                  grids_fhsegg17[[3]], grids_fhsegg17[[4]], grids_fhsegg17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhsegg3), fixed = T)
df_fhsegg_avg3_gfdl126 <- data.frame(lat = df_fhsegg3$lat, 
                                    lon = df_fhsegg3$lon, 
                                    dist = df_fhsegg3$dist,
                                    avg_pred = rowSums(df_fhsegg3[, x])/30)
saveRDS(df_fhsegg_avg3_gfdl126, file = here("data", "df_fhsegg_avg3_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg_avg3_gfdl126, "Forecasted Distribution 2070 - 2099 \n GFDL SSP126")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_gfdl_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

##### GFDL 585 -----------------------------------------------------------------------------------------------------------------
temps_gfdl_ssp585 <- readRDS(here('data', 'temps_gfdl_ssp585.rds'))
salts_gfdl_ssp585 <- readRDS(here('data', 'salts_gfdl_ssp585.rds'))

## 2015 - 2039
grids_fhsegg1 <- pred_loop(2015:2019, fhs_egg, 130, 
                          temps_gfdl_ssp585,
                          salts_gfdl_ssp585, 1,
                          egg_formula)
grids_fhsegg2 <- pred_loop(2020:2024, fhs_egg, 130, 
                          temps_gfdl_ssp585,
                          salts_gfdl_ssp585, 2,
                          egg_formula)
grids_fhsegg3 <- pred_loop(2025:2029, fhs_egg, 130,
                          temps_gfdl_ssp585,
                          salts_gfdl_ssp585, 3,
                          egg_formula)
grids_fhsegg4 <- pred_loop(2030:2034, fhs_egg, 130, 
                          temps_gfdl_ssp585,
                          salts_gfdl_ssp585, 4,
                          egg_formula)
grids_fhsegg5 <- pred_loop(2035:2039, fhs_egg, 130, 
                          temps_gfdl_ssp585,
                          salts_gfdl_ssp585, 5,
                          egg_formula)

# Combine into one data frame
df_fhsegg4 <- list(grids_fhsegg1[[1]], grids_fhsegg1[[2]], grids_fhsegg1[[3]], 
                  grids_fhsegg1[[4]], grids_fhsegg1[[5]], grids_fhsegg2[[1]], 
                  grids_fhsegg2[[2]], grids_fhsegg2[[3]], grids_fhsegg2[[4]],
                  grids_fhsegg2[[5]], grids_fhsegg3[[1]], grids_fhsegg3[[2]], 
                  grids_fhsegg3[[3]], grids_fhsegg3[[4]], grids_fhsegg3[[5]],
                  grids_fhsegg4[[1]], grids_fhsegg4[[2]], grids_fhsegg4[[3]], 
                  grids_fhsegg4[[4]], grids_fhsegg4[[5]], grids_fhsegg5[[1]], 
                  grids_fhsegg5[[2]], grids_fhsegg5[[3]], grids_fhsegg5[[4]],
                  grids_fhsegg5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhsegg4), fixed = T)
df_fhsegg_avg4_gfdl585 <- data.frame(lat = df_fhsegg4$lat, 
                                    lon = df_fhsegg4$lon, 
                                    dist = df_fhsegg4$dist,
                                    avg_pred = rowSums(df_fhsegg4[, x])/25)
saveRDS(df_fhsegg_avg4_gfdl585, file = here("data", "df_fhsegg_avg4_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg_avg4_gfdl585, "Forecasted Distribution 2015 - 2039 \n GFDL SSP585")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_gfdl_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_fhsegg6 <- pred_loop(2040:2044, fhs_egg, 130, 
                          temps_gfdl_ssp585,
                          salts_gfdl_ssp585, 6,
                          egg_formula)
grids_fhsegg7 <- pred_loop(2045:2049, fhs_egg, 130, 
                          temps_gfdl_ssp585,
                          salts_gfdl_ssp585, 7,
                          egg_formula)
grids_fhsegg8 <- pred_loop(2050:2054, fhs_egg, 130, 
                          temps_gfdl_ssp585,
                          salts_gfdl_ssp585, 8,
                          egg_formula)
grids_fhsegg9 <- pred_loop(2055:2059, fhs_egg, 130, 
                          temps_gfdl_ssp585,
                          salts_gfdl_ssp585, 9,
                          egg_formula)
grids_fhsegg10 <- pred_loop(2060:2064, fhs_egg, 130, 
                           temps_gfdl_ssp585,
                           salts_gfdl_ssp585, 10,
                           egg_formula)
grids_fhsegg11 <- pred_loop(2065:2069, fhs_egg, 130, 
                           temps_gfdl_ssp585,
                           salts_gfdl_ssp585, 11,
                           egg_formula)

# Combine into one data frame
df_fhsegg5 <- list(grids_fhsegg6[[1]], grids_fhsegg6[[2]], grids_fhsegg6[[3]], 
                  grids_fhsegg6[[4]], grids_fhsegg6[[5]], grids_fhsegg7[[1]], 
                  grids_fhsegg7[[2]], grids_fhsegg7[[3]], grids_fhsegg7[[4]],
                  grids_fhsegg7[[5]], grids_fhsegg8[[1]], grids_fhsegg8[[2]], 
                  grids_fhsegg8[[3]], grids_fhsegg8[[4]], grids_fhsegg8[[5]],
                  grids_fhsegg9[[1]], grids_fhsegg9[[2]], grids_fhsegg9[[3]], 
                  grids_fhsegg9[[4]], grids_fhsegg9[[5]], grids_fhsegg10[[1]], 
                  grids_fhsegg10[[2]], grids_fhsegg10[[3]], grids_fhsegg10[[4]],
                  grids_fhsegg11[[5]], grids_fhsegg11[[1]], grids_fhsegg11[[2]],
                  grids_fhsegg11[[3]], grids_fhsegg11[[4]], grids_fhsegg11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhsegg5), fixed = T)
df_fhsegg_avg5_gfdl585 <- data.frame(lat = df_fhsegg5$lat, 
                                    lon = df_fhsegg5$lon, 
                                    dist = df_fhsegg5$dist,
                                    avg_pred = rowSums(df_fhsegg5[, x])/30)
saveRDS(df_fhsegg_avg5_gfdl585, file = here("data", "df_fhsegg_avg5_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg_avg5_gfdl585, "Forecasted Distribution 2040 - 2069 \n GFDL SSP585")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_gfdl_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_fhsegg12 <- pred_loop(2070:2074, fhs_egg, 130,
                           temps_gfdl_ssp585,
                           salts_gfdl_ssp585, 12,
                           egg_formula)
grids_fhsegg13 <- pred_loop(2075:2079, fhs_egg, 130,
                           temps_gfdl_ssp585,
                           salts_gfdl_ssp585, 13,
                           egg_formula)
grids_fhsegg14 <- pred_loop(2080:2084, fhs_egg, 130,
                           temps_gfdl_ssp585,
                           salts_gfdl_ssp585, 14,
                           egg_formula)
grids_fhsegg15 <- pred_loop(2085:2089, fhs_egg, 130,
                           temps_gfdl_ssp585,
                           salts_gfdl_ssp585, 15,
                           egg_formula)
grids_fhsegg16 <- pred_loop(2090:2094, fhs_egg, 130,
                           temps_gfdl_ssp585,
                           salts_gfdl_ssp585, 16,
                           egg_formula)
grids_fhsegg17 <- pred_loop(2095:2099, fhs_egg, 130, 
                           temps_gfdl_ssp585,
                           salts_gfdl_ssp585, 17,
                           egg_formula)

# Combine into one data frame
df_fhsegg6 <- list(grids_fhsegg12[[1]], grids_fhsegg12[[2]], grids_fhsegg12[[3]], 
                  grids_fhsegg12[[4]], grids_fhsegg12[[5]], grids_fhsegg13[[1]], 
                  grids_fhsegg13[[2]], grids_fhsegg13[[3]], grids_fhsegg13[[4]],
                  grids_fhsegg13[[5]], grids_fhsegg14[[1]], grids_fhsegg14[[2]], 
                  grids_fhsegg14[[3]], grids_fhsegg14[[4]], grids_fhsegg14[[5]],
                  grids_fhsegg15[[1]], grids_fhsegg15[[2]], grids_fhsegg15[[3]], 
                  grids_fhsegg15[[4]], grids_fhsegg15[[5]], grids_fhsegg16[[1]], 
                  grids_fhsegg16[[2]], grids_fhsegg16[[3]], grids_fhsegg16[[4]],
                  grids_fhsegg17[[5]], grids_fhsegg17[[1]], grids_fhsegg17[[2]],
                  grids_fhsegg17[[3]], grids_fhsegg17[[4]], grids_fhsegg17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhsegg6), fixed = T)
df_fhsegg_avg6_gfdl585 <- data.frame(lat = df_fhsegg6$lat, 
                                    lon = df_fhsegg6$lon, 
                                    dist = df_fhsegg6$dist,
                                    avg_pred = rowSums(df_fhsegg6[, x])/30)
saveRDS(df_fhsegg_avg6_gfdl585, file = here("data", "df_fhsegg_avg6_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg_avg6_gfdl585, "Forecasted Distribution 2070 - 2099 \n GFDL SSP585")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_gfdl_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### MIROC 126 ----------------------------------------------------------------------------------------------------------------------------
temps_miroc_ssp126 <- readRDS(here('data', 'temps_miroc_ssp126.rds'))
salts_miroc_ssp126 <- readRDS(here('data', 'salts_miroc_ssp126.rds'))

## 2015 - 2039
grids_fhsegg1 <- pred_loop(2015:2019, fhs_egg, 130, 
                          temps_miroc_ssp126,
                          salts_miroc_ssp126, 1,
                          egg_formula)
grids_fhsegg2 <- pred_loop(2020:2024, fhs_egg, 130, 
                          temps_miroc_ssp126,
                          salts_miroc_ssp126, 2,
                          egg_formula)
grids_fhsegg3 <- pred_loop(2025:2029, fhs_egg, 130,
                          temps_miroc_ssp126,
                          salts_miroc_ssp126, 3,
                          egg_formula)
grids_fhsegg4 <- pred_loop(2030:2034, fhs_egg, 130, 
                          temps_miroc_ssp126,
                          salts_miroc_ssp126, 4,
                          egg_formula)
grids_fhsegg5 <- pred_loop(2035:2039, fhs_egg, 130, 
                          temps_miroc_ssp126,
                          salts_miroc_ssp126, 5,
                          egg_formula)

# Combine into one data frame
df_fhsegg1 <- list(grids_fhsegg1[[1]], grids_fhsegg1[[2]], grids_fhsegg1[[3]], 
                  grids_fhsegg1[[4]], grids_fhsegg1[[5]], grids_fhsegg2[[1]], 
                  grids_fhsegg2[[2]], grids_fhsegg2[[3]], grids_fhsegg2[[4]],
                  grids_fhsegg2[[5]], grids_fhsegg3[[1]], grids_fhsegg3[[2]], 
                  grids_fhsegg3[[3]], grids_fhsegg3[[4]], grids_fhsegg3[[5]],
                  grids_fhsegg4[[1]], grids_fhsegg4[[2]], grids_fhsegg4[[3]], 
                  grids_fhsegg4[[4]], grids_fhsegg4[[5]], grids_fhsegg5[[1]], 
                  grids_fhsegg5[[2]], grids_fhsegg5[[3]], grids_fhsegg5[[4]],
                  grids_fhsegg5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhsegg1), fixed = T)
df_fhsegg_avg1_miroc126 <- data.frame(lat = df_fhsegg1$lat, 
                                     lon = df_fhsegg1$lon, 
                                     dist = df_fhsegg1$dist,
                                     avg_pred = rowSums(df_fhsegg1[, x])/25)
saveRDS(df_fhsegg_avg1_miroc126, file = here("data", "df_fhsegg_avg1_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg_avg1_miroc126, "Forecasted Distribution 2015 - 2039 \n MIROC SSP126")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_miroc_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_fhsegg6 <- pred_loop(2040:2044, fhs_egg, 130, 
                          temps_miroc_ssp126,
                          salts_miroc_ssp126, 6,
                          egg_formula)
grids_fhsegg7 <- pred_loop(2045:2049, fhs_egg, 130, 
                          temps_miroc_ssp126,
                          salts_miroc_ssp126, 7,
                          egg_formula)
grids_fhsegg8 <- pred_loop(2050:2054, fhs_egg, 130, 
                          temps_miroc_ssp126,
                          salts_miroc_ssp126, 8,
                          egg_formula)
grids_fhsegg9 <- pred_loop(2055:2059, fhs_egg, 130, 
                          temps_miroc_ssp126,
                          salts_miroc_ssp126, 9,
                          egg_formula)
grids_fhsegg10 <- pred_loop(2060:2064, fhs_egg, 130, 
                           temps_miroc_ssp126,
                           salts_miroc_ssp126, 10,
                           egg_formula)
grids_fhsegg11 <- pred_loop(2065:2069, fhs_egg, 130, 
                           temps_miroc_ssp126,
                           salts_miroc_ssp126, 11,
                           egg_formula)

# Combine into one data frame
df_fhsegg2 <- list(grids_fhsegg6[[1]], grids_fhsegg6[[2]], grids_fhsegg6[[3]], 
                  grids_fhsegg6[[4]], grids_fhsegg6[[5]], grids_fhsegg7[[1]], 
                  grids_fhsegg7[[2]], grids_fhsegg7[[3]], grids_fhsegg7[[4]],
                  grids_fhsegg7[[5]], grids_fhsegg8[[1]], grids_fhsegg8[[2]], 
                  grids_fhsegg8[[3]], grids_fhsegg8[[4]], grids_fhsegg8[[5]],
                  grids_fhsegg9[[1]], grids_fhsegg9[[2]], grids_fhsegg9[[3]], 
                  grids_fhsegg9[[4]], grids_fhsegg9[[5]], grids_fhsegg10[[1]], 
                  grids_fhsegg10[[2]], grids_fhsegg10[[3]], grids_fhsegg10[[4]],
                  grids_fhsegg11[[5]], grids_fhsegg11[[1]], grids_fhsegg11[[2]],
                  grids_fhsegg11[[3]], grids_fhsegg11[[4]], grids_fhsegg11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhsegg2), fixed = T)
df_fhsegg_avg2_miroc126 <- data.frame(lat = df_fhsegg2$lat, 
                                     lon = df_fhsegg2$lon, 
                                     dist = df_fhsegg2$dist,
                                     avg_pred = rowSums(df_fhsegg2[, x])/30)
saveRDS(df_fhsegg_avg2_miroc126, file = here("data", "df_fhsegg_avg2_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg_avg2_miroc126, "Forecasted Distribution 2040 - 2069 \n MIROC SSP126")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_miroc_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_fhsegg12 <- pred_loop(2070:2074, fhs_egg, 130,
                           temps_miroc_ssp126,
                           salts_miroc_ssp126, 12,
                           egg_formula)
grids_fhsegg13 <- pred_loop(2075:2079, fhs_egg, 130,
                           temps_miroc_ssp126,
                           salts_miroc_ssp126, 13,
                           egg_formula)
grids_fhsegg14 <- pred_loop(2080:2084, fhs_egg, 130,
                           temps_miroc_ssp126,
                           salts_miroc_ssp126, 14,
                           egg_formula)
grids_fhsegg15 <- pred_loop(2085:2089, fhs_egg, 130,
                           temps_miroc_ssp126,
                           salts_miroc_ssp126, 15,
                           egg_formula)
grids_fhsegg16 <- pred_loop(2090:2094, fhs_egg, 130,
                           temps_miroc_ssp126,
                           salts_miroc_ssp126, 16,
                           egg_formula)
grids_fhsegg17 <- pred_loop(2095:2099, fhs_egg, 130, 
                           temps_miroc_ssp126,
                           salts_miroc_ssp126, 17,
                           egg_formula)

# Combine into one data frame
df_fhsegg3 <- list(grids_fhsegg12[[1]], grids_fhsegg12[[2]], grids_fhsegg12[[3]], 
                  grids_fhsegg12[[4]], grids_fhsegg12[[5]], grids_fhsegg13[[1]], 
                  grids_fhsegg13[[2]], grids_fhsegg13[[3]], grids_fhsegg13[[4]],
                  grids_fhsegg13[[5]], grids_fhsegg14[[1]], grids_fhsegg14[[2]], 
                  grids_fhsegg14[[3]], grids_fhsegg14[[4]], grids_fhsegg14[[5]],
                  grids_fhsegg15[[1]], grids_fhsegg15[[2]], grids_fhsegg15[[3]], 
                  grids_fhsegg15[[4]], grids_fhsegg15[[5]], grids_fhsegg16[[1]], 
                  grids_fhsegg16[[2]], grids_fhsegg16[[3]], grids_fhsegg16[[4]],
                  grids_fhsegg17[[5]], grids_fhsegg17[[1]], grids_fhsegg17[[2]],
                  grids_fhsegg17[[3]], grids_fhsegg17[[4]], grids_fhsegg17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhsegg3), fixed = T)
df_fhsegg_avg3_miroc126 <- data.frame(lat = df_fhsegg3$lat, 
                                     lon = df_fhsegg3$lon, 
                                     dist = df_fhsegg3$dist,
                                     avg_pred = rowSums(df_fhsegg3[, x])/30)
saveRDS(df_fhsegg_avg3_miroc126, file = here("data", "df_fhsegg_avg3_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg_avg3_miroc126, "Forecasted Distribution 2070 - 2099 \n MIROC SSP126")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_miroc_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

##### MIROC 585 -----------------------------------------------------------------------------------------------------------------
temps_miroc_ssp585 <- readRDS(here('data', 'temps_miroc_ssp585.rds'))
salts_miroc_ssp585 <- readRDS(here('data', 'salts_miroc_ssp585.rds'))

## 2015 - 2039
grids_fhsegg1 <- pred_loop(2015:2019, fhs_egg, 130, 
                          temps_miroc_ssp585,
                          salts_miroc_ssp585, 1,
                          egg_formula)
grids_fhsegg2 <- pred_loop(2020:2024, fhs_egg, 130, 
                          temps_miroc_ssp585,
                          salts_miroc_ssp585, 2,
                          egg_formula)
grids_fhsegg3 <- pred_loop(2025:2029, fhs_egg, 130,
                          temps_miroc_ssp585,
                          salts_miroc_ssp585, 3,
                          egg_formula)
grids_fhsegg4 <- pred_loop(2030:2034, fhs_egg, 130, 
                          temps_miroc_ssp585,
                          salts_miroc_ssp585, 4,
                          egg_formula)
grids_fhsegg5 <- pred_loop(2035:2039, fhs_egg, 130, 
                          temps_miroc_ssp585,
                          salts_miroc_ssp585, 5,
                          egg_formula)

# Combine into one data frame
df_fhsegg4 <- list(grids_fhsegg1[[1]], grids_fhsegg1[[2]], grids_fhsegg1[[3]], 
                  grids_fhsegg1[[4]], grids_fhsegg1[[5]], grids_fhsegg2[[1]], 
                  grids_fhsegg2[[2]], grids_fhsegg2[[3]], grids_fhsegg2[[4]],
                  grids_fhsegg2[[5]], grids_fhsegg3[[1]], grids_fhsegg3[[2]], 
                  grids_fhsegg3[[3]], grids_fhsegg3[[4]], grids_fhsegg3[[5]],
                  grids_fhsegg4[[1]], grids_fhsegg4[[2]], grids_fhsegg4[[3]], 
                  grids_fhsegg4[[4]], grids_fhsegg4[[5]], grids_fhsegg5[[1]], 
                  grids_fhsegg5[[2]], grids_fhsegg5[[3]], grids_fhsegg5[[4]],
                  grids_fhsegg5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhsegg4), fixed = T)
df_fhsegg_avg4_miroc585 <- data.frame(lat = df_fhsegg4$lat, 
                                     lon = df_fhsegg4$lon, 
                                     dist = df_fhsegg4$dist,
                                     avg_pred = rowSums(df_fhsegg4[, x])/25)
saveRDS(df_fhsegg_avg4_miroc585, file = here("data", "df_fhsegg_avg4_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg_avg4_miroc585, "Forecasted Distribution 2015 - 2039 \n MIROC SSP585")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_miroc_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_fhsegg6 <- pred_loop(2040:2044, fhs_egg, 130, 
                          temps_miroc_ssp585,
                          salts_miroc_ssp585, 6,
                          egg_formula)
grids_fhsegg7 <- pred_loop(2045:2049, fhs_egg, 130, 
                          temps_miroc_ssp585,
                          salts_miroc_ssp585, 7,
                          egg_formula)
grids_fhsegg8 <- pred_loop(2050:2054, fhs_egg, 130, 
                          temps_miroc_ssp585,
                          salts_miroc_ssp585, 8,
                          egg_formula)
grids_fhsegg9 <- pred_loop(2055:2059, fhs_egg, 130, 
                          temps_miroc_ssp585,
                          salts_miroc_ssp585, 9,
                          egg_formula)
grids_fhsegg10 <- pred_loop(2060:2064, fhs_egg, 130, 
                           temps_miroc_ssp585,
                           salts_miroc_ssp585, 10,
                           egg_formula)
grids_fhsegg11 <- pred_loop(2065:2069, fhs_egg, 130, 
                           temps_miroc_ssp585,
                           salts_miroc_ssp585, 11,
                           egg_formula)

# Combine into one data frame
df_fhsegg5 <- list(grids_fhsegg6[[1]], grids_fhsegg6[[2]], grids_fhsegg6[[3]], 
                  grids_fhsegg6[[4]], grids_fhsegg6[[5]], grids_fhsegg7[[1]], 
                  grids_fhsegg7[[2]], grids_fhsegg7[[3]], grids_fhsegg7[[4]],
                  grids_fhsegg7[[5]], grids_fhsegg8[[1]], grids_fhsegg8[[2]], 
                  grids_fhsegg8[[3]], grids_fhsegg8[[4]], grids_fhsegg8[[5]],
                  grids_fhsegg9[[1]], grids_fhsegg9[[2]], grids_fhsegg9[[3]], 
                  grids_fhsegg9[[4]], grids_fhsegg9[[5]], grids_fhsegg10[[1]], 
                  grids_fhsegg10[[2]], grids_fhsegg10[[3]], grids_fhsegg10[[4]],
                  grids_fhsegg11[[5]], grids_fhsegg11[[1]], grids_fhsegg11[[2]],
                  grids_fhsegg11[[3]], grids_fhsegg11[[4]], grids_fhsegg11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhsegg5), fixed = T)
df_fhsegg_avg5_miroc585 <- data.frame(lat = df_fhsegg5$lat, 
                                     lon = df_fhsegg5$lon, 
                                     dist = df_fhsegg5$dist,
                                     avg_pred = rowSums(df_fhsegg5[, x])/30)
saveRDS(df_fhsegg_avg5_miroc585, file = here("data", "df_fhsegg_avg5_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg_avg5_miroc585, "Forecasted Distribution 2040 - 2069 \n MIROC SSP585")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_miroc_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_fhsegg12 <- pred_loop(2070:2074, fhs_egg, 130,
                           temps_miroc_ssp585,
                           salts_miroc_ssp585, 12,
                           egg_formula)
grids_fhsegg13 <- pred_loop(2075:2079, fhs_egg, 130,
                           temps_miroc_ssp585,
                           salts_miroc_ssp585, 13,
                           egg_formula)
grids_fhsegg14 <- pred_loop(2080:2084, fhs_egg, 130,
                           temps_miroc_ssp585,
                           salts_miroc_ssp585, 14,
                           egg_formula)
grids_fhsegg15 <- pred_loop(2085:2089, fhs_egg, 130,
                           temps_miroc_ssp585,
                           salts_miroc_ssp585, 15,
                           egg_formula)
grids_fhsegg16 <- pred_loop(2090:2094, fhs_egg, 130,
                           temps_miroc_ssp585,
                           salts_miroc_ssp585, 16,
                           egg_formula)
grids_fhsegg17 <- pred_loop(2095:2099, fhs_egg, 130, 
                           temps_miroc_ssp585,
                           salts_miroc_ssp585, 17,
                           egg_formula)

# Combine into one data frame
df_fhsegg6 <- list(grids_fhsegg12[[1]], grids_fhsegg12[[2]], grids_fhsegg12[[3]], 
                  grids_fhsegg12[[4]], grids_fhsegg12[[5]], grids_fhsegg13[[1]], 
                  grids_fhsegg13[[2]], grids_fhsegg13[[3]], grids_fhsegg13[[4]],
                  grids_fhsegg13[[5]], grids_fhsegg14[[1]], grids_fhsegg14[[2]], 
                  grids_fhsegg14[[3]], grids_fhsegg14[[4]], grids_fhsegg14[[5]],
                  grids_fhsegg15[[1]], grids_fhsegg15[[2]], grids_fhsegg15[[3]], 
                  grids_fhsegg15[[4]], grids_fhsegg15[[5]], grids_fhsegg16[[1]], 
                  grids_fhsegg16[[2]], grids_fhsegg16[[3]], grids_fhsegg16[[4]],
                  grids_fhsegg17[[5]], grids_fhsegg17[[1]], grids_fhsegg17[[2]],
                  grids_fhsegg17[[3]], grids_fhsegg17[[4]], grids_fhsegg17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhsegg6), fixed = T)
df_fhsegg_avg6_miroc585 <- data.frame(lat = df_fhsegg6$lat, 
                                     lon = df_fhsegg6$lon, 
                                     dist = df_fhsegg6$dist,
                                     avg_pred = rowSums(df_fhsegg6[, x])/30)
saveRDS(df_fhsegg_avg6_miroc585, file = here("data", "df_fhsegg_avg6_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg_avg6_miroc585, "Forecasted Distribution 2070 - 2099 \n MIROC SSP585")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_miroc_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()



### Flathead Larvae --------------------------------------------------------------------------------------------------------------------------
#### Forecast and average into 3 time periods ---------------------------------------------------------------------------------------------
##### CESM 126 ----------------------------------------------------------------------------------------------------------------------------
temps_cesm_ssp126 <- readRDS(here('data', 'temps_cesm_ssp126.rds'))
salts_cesm_ssp126 <- readRDS(here('data', 'salts_cesm_ssp126.rds'))

## 2015 - 2039
grids_fhslarvae1 <- pred_loop(2015:2019, fhs_larvae, 130, 
                             temps_cesm_ssp126,
                             salts_cesm_ssp126, 1,
                             larval_formula)
grids_fhslarvae2 <- pred_loop(2020:2024, fhs_larvae, 130, 
                             temps_cesm_ssp126,
                             salts_cesm_ssp126, 2,
                             larval_formula)
grids_fhslarvae3 <- pred_loop(2025:2029, fhs_larvae, 130,
                             temps_cesm_ssp126,
                             salts_cesm_ssp126, 3,
                             larval_formula)
grids_fhslarvae4 <- pred_loop(2030:2034, fhs_larvae, 130, 
                             temps_cesm_ssp126,
                             salts_cesm_ssp126, 4,
                             larval_formula)
grids_fhslarvae5 <- pred_loop(2035:2039, fhs_larvae, 130, 
                             temps_cesm_ssp126,
                             salts_cesm_ssp126, 5,
                             larval_formula)

# Combine into one data frame
df_fhslarvae1 <- list(grids_fhslarvae1[[1]], grids_fhslarvae1[[2]], grids_fhslarvae1[[3]], 
                     grids_fhslarvae1[[4]], grids_fhslarvae1[[5]], grids_fhslarvae2[[1]], 
                     grids_fhslarvae2[[2]], grids_fhslarvae2[[3]], grids_fhslarvae2[[4]],
                     grids_fhslarvae2[[5]], grids_fhslarvae3[[1]], grids_fhslarvae3[[2]], 
                     grids_fhslarvae3[[3]], grids_fhslarvae3[[4]], grids_fhslarvae3[[5]],
                     grids_fhslarvae4[[1]], grids_fhslarvae4[[2]], grids_fhslarvae4[[3]], 
                     grids_fhslarvae4[[4]], grids_fhslarvae4[[5]], grids_fhslarvae5[[1]], 
                     grids_fhslarvae5[[2]], grids_fhslarvae5[[3]], grids_fhslarvae5[[4]],
                     grids_fhslarvae5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhslarvae1), fixed = T)
df_fhslarvae_avg1_cesm126 <- data.frame(lat = df_fhslarvae1$lat, 
                                       lon = df_fhslarvae1$lon, 
                                       dist = df_fhslarvae1$dist,
                                       avg_pred = rowSums(df_fhslarvae1[, x])/25)
saveRDS(df_fhslarvae_avg1_cesm126, file = here("data", "df_fhslarvae_avg1_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae_avg1_cesm126, "Forecasted Distribution 2015 - 2039 \n CESM SSP126")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_cesm_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_fhslarvae6 <- pred_loop(2040:2044, fhs_larvae, 130, 
                             temps_cesm_ssp126,
                             salts_cesm_ssp126, 6,
                             larval_formula)
grids_fhslarvae7 <- pred_loop(2045:2049, fhs_larvae, 130, 
                             temps_cesm_ssp126,
                             salts_cesm_ssp126, 7,
                             larval_formula)
grids_fhslarvae8 <- pred_loop(2050:2054, fhs_larvae, 130, 
                             temps_cesm_ssp126,
                             salts_cesm_ssp126, 8,
                             larval_formula)
grids_fhslarvae9 <- pred_loop(2055:2059, fhs_larvae, 130, 
                             temps_cesm_ssp126,
                             salts_cesm_ssp126, 9,
                             larval_formula)
grids_fhslarvae10 <- pred_loop(2060:2064, fhs_larvae, 130, 
                              temps_cesm_ssp126,
                              salts_cesm_ssp126, 10,
                              larval_formula)
grids_fhslarvae11 <- pred_loop(2065:2069, fhs_larvae, 130, 
                              temps_cesm_ssp126,
                              salts_cesm_ssp126, 11,
                              larval_formula)

# Combine into one data frame
df_fhslarvae2 <- list(grids_fhslarvae6[[1]], grids_fhslarvae6[[2]], grids_fhslarvae6[[3]], 
                     grids_fhslarvae6[[4]], grids_fhslarvae6[[5]], grids_fhslarvae7[[1]], 
                     grids_fhslarvae7[[2]], grids_fhslarvae7[[3]], grids_fhslarvae7[[4]],
                     grids_fhslarvae7[[5]], grids_fhslarvae8[[1]], grids_fhslarvae8[[2]], 
                     grids_fhslarvae8[[3]], grids_fhslarvae8[[4]], grids_fhslarvae8[[5]],
                     grids_fhslarvae9[[1]], grids_fhslarvae9[[2]], grids_fhslarvae9[[3]], 
                     grids_fhslarvae9[[4]], grids_fhslarvae9[[5]], grids_fhslarvae10[[1]], 
                     grids_fhslarvae10[[2]], grids_fhslarvae10[[3]], grids_fhslarvae10[[4]],
                     grids_fhslarvae11[[5]], grids_fhslarvae11[[1]], grids_fhslarvae11[[2]],
                     grids_fhslarvae11[[3]], grids_fhslarvae11[[4]], grids_fhslarvae11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhslarvae2), fixed = T)
df_fhslarvae_avg2_cesm126 <- data.frame(lat = df_fhslarvae2$lat, 
                                       lon = df_fhslarvae2$lon, 
                                       dist = df_fhslarvae2$dist,
                                       avg_pred = rowSums(df_fhslarvae2[, x])/30)
saveRDS(df_fhslarvae_avg2_cesm126, file = here("data", "df_fhslarvae_avg2_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae_avg2_cesm126, "Forecasted Distribution 2040 - 2069 \n CESM SSP126")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_cesm_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_fhslarvae12 <- pred_loop(2070:2074, fhs_larvae, 130,
                              temps_cesm_ssp126,
                              salts_cesm_ssp126, 12,
                              larval_formula)
grids_fhslarvae13 <- pred_loop(2075:2079, fhs_larvae, 130,
                              temps_cesm_ssp126,
                              salts_cesm_ssp126, 13,
                              larval_formula)
grids_fhslarvae14 <- pred_loop(2080:2084, fhs_larvae, 130,
                              temps_cesm_ssp126,
                              salts_cesm_ssp126, 14,
                              larval_formula)
grids_fhslarvae15 <- pred_loop(2085:2089, fhs_larvae, 130,
                              temps_cesm_ssp126,
                              salts_cesm_ssp126, 15,
                              larval_formula)
grids_fhslarvae16 <- pred_loop(2090:2094, fhs_larvae, 130,
                              temps_cesm_ssp126,
                              salts_cesm_ssp126, 16,
                              larval_formula)
grids_fhslarvae17 <- pred_loop(2095:2099, fhs_larvae, 130, 
                              temps_cesm_ssp126,
                              salts_cesm_ssp126, 17,
                              larval_formula)

# Combine into one data frame
df_fhslarvae3 <- list(grids_fhslarvae12[[1]], grids_fhslarvae12[[2]], grids_fhslarvae12[[3]], 
                     grids_fhslarvae12[[4]], grids_fhslarvae12[[5]], grids_fhslarvae13[[1]], 
                     grids_fhslarvae13[[2]], grids_fhslarvae13[[3]], grids_fhslarvae13[[4]],
                     grids_fhslarvae13[[5]], grids_fhslarvae14[[1]], grids_fhslarvae14[[2]], 
                     grids_fhslarvae14[[3]], grids_fhslarvae14[[4]], grids_fhslarvae14[[5]],
                     grids_fhslarvae15[[1]], grids_fhslarvae15[[2]], grids_fhslarvae15[[3]], 
                     grids_fhslarvae15[[4]], grids_fhslarvae15[[5]], grids_fhslarvae16[[1]], 
                     grids_fhslarvae16[[2]], grids_fhslarvae16[[3]], grids_fhslarvae16[[4]],
                     grids_fhslarvae17[[5]], grids_fhslarvae17[[1]], grids_fhslarvae17[[2]],
                     grids_fhslarvae17[[3]], grids_fhslarvae17[[4]], grids_fhslarvae17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhslarvae3), fixed = T)
df_fhslarvae_avg3_cesm126 <- data.frame(lat = df_fhslarvae3$lat, 
                                       lon = df_fhslarvae3$lon, 
                                       dist = df_fhslarvae3$dist,
                                       avg_pred = rowSums(df_fhslarvae3[, x])/30)
saveRDS(df_fhslarvae_avg3_cesm126, file = here("data", "df_fhslarvae_avg3_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae_avg3_cesm126, "Forecasted Distribution 2070 - 2099 \n CESM SSP126")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_cesm_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

##### CESM 585 -----------------------------------------------------------------------------------------------------------------
temps_cesm_ssp585 <- readRDS(here('data', 'temps_cesm_ssp585.rds'))
salts_cesm_ssp585 <- readRDS(here('data', 'salts_cesm_ssp585.rds'))

## 2015 - 2039
grids_fhslarvae1 <- pred_loop(2015:2019, fhs_larvae, 130, 
                             temps_cesm_ssp585,
                             salts_cesm_ssp585, 1,
                             larval_formula)
grids_fhslarvae2 <- pred_loop(2020:2024, fhs_larvae, 130, 
                             temps_cesm_ssp585,
                             salts_cesm_ssp585, 2,
                             larval_formula)
grids_fhslarvae3 <- pred_loop(2025:2029, fhs_larvae, 130,
                             temps_cesm_ssp585,
                             salts_cesm_ssp585, 3,
                             larval_formula)
grids_fhslarvae4 <- pred_loop(2030:2034, fhs_larvae, 130, 
                             temps_cesm_ssp585,
                             salts_cesm_ssp585, 4,
                             larval_formula)
grids_fhslarvae5 <- pred_loop(2035:2039, fhs_larvae, 130, 
                             temps_cesm_ssp585,
                             salts_cesm_ssp585, 5,
                             larval_formula)

# Combine into one data frame
df_fhslarvae4 <- list(grids_fhslarvae1[[1]], grids_fhslarvae1[[2]], grids_fhslarvae1[[3]], 
                     grids_fhslarvae1[[4]], grids_fhslarvae1[[5]], grids_fhslarvae2[[1]], 
                     grids_fhslarvae2[[2]], grids_fhslarvae2[[3]], grids_fhslarvae2[[4]],
                     grids_fhslarvae2[[5]], grids_fhslarvae3[[1]], grids_fhslarvae3[[2]], 
                     grids_fhslarvae3[[3]], grids_fhslarvae3[[4]], grids_fhslarvae3[[5]],
                     grids_fhslarvae4[[1]], grids_fhslarvae4[[2]], grids_fhslarvae4[[3]], 
                     grids_fhslarvae4[[4]], grids_fhslarvae4[[5]], grids_fhslarvae5[[1]], 
                     grids_fhslarvae5[[2]], grids_fhslarvae5[[3]], grids_fhslarvae5[[4]],
                     grids_fhslarvae5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhslarvae4), fixed = T)
df_fhslarvae_avg4_cesm585 <- data.frame(lat = df_fhslarvae4$lat, 
                                       lon = df_fhslarvae4$lon, 
                                       dist = df_fhslarvae4$dist,
                                       avg_pred = rowSums(df_fhslarvae4[, x])/25)
saveRDS(df_fhslarvae_avg4_cesm585, file = here("data", "df_fhslarvae_avg4_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae_avg4_cesm585, "Forecasted Distribution 2015 - 2039 \n CESM SSP585")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_cesm_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_fhslarvae6 <- pred_loop(2040:2044, fhs_larvae, 130, 
                             temps_cesm_ssp585,
                             salts_cesm_ssp585, 6,
                             larval_formula)
grids_fhslarvae7 <- pred_loop(2045:2049, fhs_larvae, 130, 
                             temps_cesm_ssp585,
                             salts_cesm_ssp585, 7,
                             larval_formula)
grids_fhslarvae8 <- pred_loop(2050:2054, fhs_larvae, 130, 
                             temps_cesm_ssp585,
                             salts_cesm_ssp585, 8,
                             larval_formula)
grids_fhslarvae9 <- pred_loop(2055:2059, fhs_larvae, 130, 
                             temps_cesm_ssp585,
                             salts_cesm_ssp585, 9,
                             larval_formula)
grids_fhslarvae10 <- pred_loop(2060:2064, fhs_larvae, 130, 
                              temps_cesm_ssp585,
                              salts_cesm_ssp585, 10,
                              larval_formula)
grids_fhslarvae11 <- pred_loop(2065:2069, fhs_larvae, 130, 
                              temps_cesm_ssp585,
                              salts_cesm_ssp585, 11,
                              larval_formula)

# Combine into one data frame
df_fhslarvae5 <- list(grids_fhslarvae6[[1]], grids_fhslarvae6[[2]], grids_fhslarvae6[[3]], 
                     grids_fhslarvae6[[4]], grids_fhslarvae6[[5]], grids_fhslarvae7[[1]], 
                     grids_fhslarvae7[[2]], grids_fhslarvae7[[3]], grids_fhslarvae7[[4]],
                     grids_fhslarvae7[[5]], grids_fhslarvae8[[1]], grids_fhslarvae8[[2]], 
                     grids_fhslarvae8[[3]], grids_fhslarvae8[[4]], grids_fhslarvae8[[5]],
                     grids_fhslarvae9[[1]], grids_fhslarvae9[[2]], grids_fhslarvae9[[3]], 
                     grids_fhslarvae9[[4]], grids_fhslarvae9[[5]], grids_fhslarvae10[[1]], 
                     grids_fhslarvae10[[2]], grids_fhslarvae10[[3]], grids_fhslarvae10[[4]],
                     grids_fhslarvae11[[5]], grids_fhslarvae11[[1]], grids_fhslarvae11[[2]],
                     grids_fhslarvae11[[3]], grids_fhslarvae11[[4]], grids_fhslarvae11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhslarvae5), fixed = T)
df_fhslarvae_avg5_cesm585 <- data.frame(lat = df_fhslarvae5$lat, 
                                       lon = df_fhslarvae5$lon, 
                                       dist = df_fhslarvae5$dist,
                                       avg_pred = rowSums(df_fhslarvae5[, x])/30)
saveRDS(df_fhslarvae_avg5_cesm585, file = here("data", "df_fhslarvae_avg5_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae_avg5_cesm585, "Forecasted Distribution 2040 - 2069 \n CESM SSP585")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_cesm_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_fhslarvae12 <- pred_loop(2070:2074, fhs_larvae, 130,
                              temps_cesm_ssp585,
                              salts_cesm_ssp585, 12,
                              larval_formula)
grids_fhslarvae13 <- pred_loop(2075:2079, fhs_larvae, 130,
                              temps_cesm_ssp585,
                              salts_cesm_ssp585, 13,
                              larval_formula)
grids_fhslarvae14 <- pred_loop(2080:2084, fhs_larvae, 130,
                              temps_cesm_ssp585,
                              salts_cesm_ssp585, 14,
                              larval_formula)
grids_fhslarvae15 <- pred_loop(2085:2089, fhs_larvae, 130,
                              temps_cesm_ssp585,
                              salts_cesm_ssp585, 15,
                              larval_formula)
grids_fhslarvae16 <- pred_loop(2090:2094, fhs_larvae, 130,
                              temps_cesm_ssp585,
                              salts_cesm_ssp585, 16,
                              larval_formula)
grids_fhslarvae17 <- pred_loop(2095:2099, fhs_larvae, 130, 
                              temps_cesm_ssp585,
                              salts_cesm_ssp585, 17,
                              larval_formula)

# Combine into one data frame
df_fhslarvae6 <- list(grids_fhslarvae12[[1]], grids_fhslarvae12[[2]], grids_fhslarvae12[[3]], 
                     grids_fhslarvae12[[4]], grids_fhslarvae12[[5]], grids_fhslarvae13[[1]], 
                     grids_fhslarvae13[[2]], grids_fhslarvae13[[3]], grids_fhslarvae13[[4]],
                     grids_fhslarvae13[[5]], grids_fhslarvae14[[1]], grids_fhslarvae14[[2]], 
                     grids_fhslarvae14[[3]], grids_fhslarvae14[[4]], grids_fhslarvae14[[5]],
                     grids_fhslarvae15[[1]], grids_fhslarvae15[[2]], grids_fhslarvae15[[3]], 
                     grids_fhslarvae15[[4]], grids_fhslarvae15[[5]], grids_fhslarvae16[[1]], 
                     grids_fhslarvae16[[2]], grids_fhslarvae16[[3]], grids_fhslarvae16[[4]],
                     grids_fhslarvae17[[5]], grids_fhslarvae17[[1]], grids_fhslarvae17[[2]],
                     grids_fhslarvae17[[3]], grids_fhslarvae17[[4]], grids_fhslarvae17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhslarvae6), fixed = T)
df_fhslarvae_avg6_cesm585 <- data.frame(lat = df_fhslarvae6$lat, 
                                       lon = df_fhslarvae6$lon, 
                                       dist = df_fhslarvae6$dist,
                                       avg_pred = rowSums(df_fhslarvae6[, x])/30)
saveRDS(df_fhslarvae_avg6_cesm585, file = here("data", "df_fhslarvae_avg6_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae_avg6_cesm585, "Forecasted Distribution 2070 - 2099 \n CESM SSP585")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_cesm_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### GFDL 126 ----------------------------------------------------------------------------------------------------------------------------
temps_gfdl_ssp126 <- readRDS(here('data', 'temps_gfdl_ssp126.rds'))
salts_gfdl_ssp126 <- readRDS(here('data', 'salts_gfdl_ssp126.rds'))

## 2015 - 2039
grids_fhslarvae1 <- pred_loop(2015:2019, fhs_larvae, 130, 
                             temps_gfdl_ssp126,
                             salts_gfdl_ssp126, 1,
                             larval_formula)
grids_fhslarvae2 <- pred_loop(2020:2024, fhs_larvae, 130, 
                             temps_gfdl_ssp126,
                             salts_gfdl_ssp126, 2,
                             larval_formula)
grids_fhslarvae3 <- pred_loop(2025:2029, fhs_larvae, 130,
                             temps_gfdl_ssp126,
                             salts_gfdl_ssp126, 3,
                             larval_formula)
grids_fhslarvae4 <- pred_loop(2030:2034, fhs_larvae, 130, 
                             temps_gfdl_ssp126,
                             salts_gfdl_ssp126, 4,
                             larval_formula)
grids_fhslarvae5 <- pred_loop(2035:2039, fhs_larvae, 130, 
                             temps_gfdl_ssp126,
                             salts_gfdl_ssp126, 5,
                             larval_formula)

# Combine into one data frame
df_fhslarvae1 <- list(grids_fhslarvae1[[1]], grids_fhslarvae1[[2]], grids_fhslarvae1[[3]], 
                     grids_fhslarvae1[[4]], grids_fhslarvae1[[5]], grids_fhslarvae2[[1]], 
                     grids_fhslarvae2[[2]], grids_fhslarvae2[[3]], grids_fhslarvae2[[4]],
                     grids_fhslarvae2[[5]], grids_fhslarvae3[[1]], grids_fhslarvae3[[2]], 
                     grids_fhslarvae3[[3]], grids_fhslarvae3[[4]], grids_fhslarvae3[[5]],
                     grids_fhslarvae4[[1]], grids_fhslarvae4[[2]], grids_fhslarvae4[[3]], 
                     grids_fhslarvae4[[4]], grids_fhslarvae4[[5]], grids_fhslarvae5[[1]], 
                     grids_fhslarvae5[[2]], grids_fhslarvae5[[3]], grids_fhslarvae5[[4]],
                     grids_fhslarvae5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhslarvae1), fixed = T)
df_fhslarvae_avg1_gfdl126 <- data.frame(lat = df_fhslarvae1$lat, 
                                       lon = df_fhslarvae1$lon, 
                                       dist = df_fhslarvae1$dist,
                                       avg_pred = rowSums(df_fhslarvae1[, x])/25)
saveRDS(df_fhslarvae_avg1_gfdl126, file = here("data", "df_fhslarvae_avg1_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae_avg1_gfdl126, "Forecasted Distribution 2015 - 2039 \n GFDL SSP126")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_gfdl_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_fhslarvae6 <- pred_loop(2040:2044, fhs_larvae, 130, 
                             temps_gfdl_ssp126,
                             salts_gfdl_ssp126, 6,
                             larval_formula)
grids_fhslarvae7 <- pred_loop(2045:2049, fhs_larvae, 130, 
                             temps_gfdl_ssp126,
                             salts_gfdl_ssp126, 7,
                             larval_formula)
grids_fhslarvae8 <- pred_loop(2050:2054, fhs_larvae, 130, 
                             temps_gfdl_ssp126,
                             salts_gfdl_ssp126, 8,
                             larval_formula)
grids_fhslarvae9 <- pred_loop(2055:2059, fhs_larvae, 130, 
                             temps_gfdl_ssp126,
                             salts_gfdl_ssp126, 9,
                             larval_formula)
grids_fhslarvae10 <- pred_loop(2060:2064, fhs_larvae, 130, 
                              temps_gfdl_ssp126,
                              salts_gfdl_ssp126, 10,
                              larval_formula)
grids_fhslarvae11 <- pred_loop(2065:2069, fhs_larvae, 130, 
                              temps_gfdl_ssp126,
                              salts_gfdl_ssp126, 11,
                              larval_formula)

# Combine into one data frame
df_fhslarvae2 <- list(grids_fhslarvae6[[1]], grids_fhslarvae6[[2]], grids_fhslarvae6[[3]], 
                     grids_fhslarvae6[[4]], grids_fhslarvae6[[5]], grids_fhslarvae7[[1]], 
                     grids_fhslarvae7[[2]], grids_fhslarvae7[[3]], grids_fhslarvae7[[4]],
                     grids_fhslarvae7[[5]], grids_fhslarvae8[[1]], grids_fhslarvae8[[2]], 
                     grids_fhslarvae8[[3]], grids_fhslarvae8[[4]], grids_fhslarvae8[[5]],
                     grids_fhslarvae9[[1]], grids_fhslarvae9[[2]], grids_fhslarvae9[[3]], 
                     grids_fhslarvae9[[4]], grids_fhslarvae9[[5]], grids_fhslarvae10[[1]], 
                     grids_fhslarvae10[[2]], grids_fhslarvae10[[3]], grids_fhslarvae10[[4]],
                     grids_fhslarvae11[[5]], grids_fhslarvae11[[1]], grids_fhslarvae11[[2]],
                     grids_fhslarvae11[[3]], grids_fhslarvae11[[4]], grids_fhslarvae11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhslarvae2), fixed = T)
df_fhslarvae_avg2_gfdl126 <- data.frame(lat = df_fhslarvae2$lat, 
                                       lon = df_fhslarvae2$lon, 
                                       dist = df_fhslarvae2$dist,
                                       avg_pred = rowSums(df_fhslarvae2[, x])/30)
saveRDS(df_fhslarvae_avg2_gfdl126, file = here("data", "df_fhslarvae_avg2_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae_avg2_gfdl126, "Forecasted Distribution 2040 - 2069 \n GFDL SSP126")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_gfdl_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_fhslarvae12 <- pred_loop(2070:2074, fhs_larvae, 130,
                              temps_gfdl_ssp126,
                              salts_gfdl_ssp126, 12,
                              larval_formula)
grids_fhslarvae13 <- pred_loop(2075:2079, fhs_larvae, 130,
                              temps_gfdl_ssp126,
                              salts_gfdl_ssp126, 13,
                              larval_formula)
grids_fhslarvae14 <- pred_loop(2080:2084, fhs_larvae, 130,
                              temps_gfdl_ssp126,
                              salts_gfdl_ssp126, 14,
                              larval_formula)
grids_fhslarvae15 <- pred_loop(2085:2089, fhs_larvae, 130,
                              temps_gfdl_ssp126,
                              salts_gfdl_ssp126, 15,
                              larval_formula)
grids_fhslarvae16 <- pred_loop(2090:2094, fhs_larvae, 130,
                              temps_gfdl_ssp126,
                              salts_gfdl_ssp126, 16,
                              larval_formula)
grids_fhslarvae17 <- pred_loop(2095:2099, fhs_larvae, 130, 
                              temps_gfdl_ssp126,
                              salts_gfdl_ssp126, 17,
                              larval_formula)

# Combine into one data frame
df_fhslarvae3 <- list(grids_fhslarvae12[[1]], grids_fhslarvae12[[2]], grids_fhslarvae12[[3]], 
                     grids_fhslarvae12[[4]], grids_fhslarvae12[[5]], grids_fhslarvae13[[1]], 
                     grids_fhslarvae13[[2]], grids_fhslarvae13[[3]], grids_fhslarvae13[[4]],
                     grids_fhslarvae13[[5]], grids_fhslarvae14[[1]], grids_fhslarvae14[[2]], 
                     grids_fhslarvae14[[3]], grids_fhslarvae14[[4]], grids_fhslarvae14[[5]],
                     grids_fhslarvae15[[1]], grids_fhslarvae15[[2]], grids_fhslarvae15[[3]], 
                     grids_fhslarvae15[[4]], grids_fhslarvae15[[5]], grids_fhslarvae16[[1]], 
                     grids_fhslarvae16[[2]], grids_fhslarvae16[[3]], grids_fhslarvae16[[4]],
                     grids_fhslarvae17[[5]], grids_fhslarvae17[[1]], grids_fhslarvae17[[2]],
                     grids_fhslarvae17[[3]], grids_fhslarvae17[[4]], grids_fhslarvae17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhslarvae3), fixed = T)
df_fhslarvae_avg3_gfdl126 <- data.frame(lat = df_fhslarvae3$lat, 
                                       lon = df_fhslarvae3$lon, 
                                       dist = df_fhslarvae3$dist,
                                       avg_pred = rowSums(df_fhslarvae3[, x])/30)
saveRDS(df_fhslarvae_avg3_gfdl126, file = here("data", "df_fhslarvae_avg3_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae_avg3_gfdl126, "Forecasted Distribution 2070 - 2099 \n GFDL SSP126")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_gfdl_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

##### GFDL 585 -----------------------------------------------------------------------------------------------------------------
temps_gfdl_ssp585 <- readRDS(here('data', 'temps_gfdl_ssp585.rds'))
salts_gfdl_ssp585 <- readRDS(here('data', 'salts_gfdl_ssp585.rds'))

## 2015 - 2039
grids_fhslarvae1 <- pred_loop(2015:2019, fhs_larvae, 130, 
                             temps_gfdl_ssp585,
                             salts_gfdl_ssp585, 1,
                             larval_formula)
grids_fhslarvae2 <- pred_loop(2020:2024, fhs_larvae, 130, 
                             temps_gfdl_ssp585,
                             salts_gfdl_ssp585, 2,
                             larval_formula)
grids_fhslarvae3 <- pred_loop(2025:2029, fhs_larvae, 130,
                             temps_gfdl_ssp585,
                             salts_gfdl_ssp585, 3,
                             larval_formula)
grids_fhslarvae4 <- pred_loop(2030:2034, fhs_larvae, 130, 
                             temps_gfdl_ssp585,
                             salts_gfdl_ssp585, 4,
                             larval_formula)
grids_fhslarvae5 <- pred_loop(2035:2039, fhs_larvae, 130, 
                             temps_gfdl_ssp585,
                             salts_gfdl_ssp585, 5,
                             larval_formula)

# Combine into one data frame
df_fhslarvae4 <- list(grids_fhslarvae1[[1]], grids_fhslarvae1[[2]], grids_fhslarvae1[[3]], 
                     grids_fhslarvae1[[4]], grids_fhslarvae1[[5]], grids_fhslarvae2[[1]], 
                     grids_fhslarvae2[[2]], grids_fhslarvae2[[3]], grids_fhslarvae2[[4]],
                     grids_fhslarvae2[[5]], grids_fhslarvae3[[1]], grids_fhslarvae3[[2]], 
                     grids_fhslarvae3[[3]], grids_fhslarvae3[[4]], grids_fhslarvae3[[5]],
                     grids_fhslarvae4[[1]], grids_fhslarvae4[[2]], grids_fhslarvae4[[3]], 
                     grids_fhslarvae4[[4]], grids_fhslarvae4[[5]], grids_fhslarvae5[[1]], 
                     grids_fhslarvae5[[2]], grids_fhslarvae5[[3]], grids_fhslarvae5[[4]],
                     grids_fhslarvae5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhslarvae4), fixed = T)
df_fhslarvae_avg4_gfdl585 <- data.frame(lat = df_fhslarvae4$lat, 
                                       lon = df_fhslarvae4$lon, 
                                       dist = df_fhslarvae4$dist,
                                       avg_pred = rowSums(df_fhslarvae4[, x])/25)
saveRDS(df_fhslarvae_avg4_gfdl585, file = here("data", "df_fhslarvae_avg4_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae_avg4_gfdl585, "Forecasted Distribution 2015 - 2039 \n GFDL SSP585")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_gfdl_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_fhslarvae6 <- pred_loop(2040:2044, fhs_larvae, 130, 
                             temps_gfdl_ssp585,
                             salts_gfdl_ssp585, 6,
                             larval_formula)
grids_fhslarvae7 <- pred_loop(2045:2049, fhs_larvae, 130, 
                             temps_gfdl_ssp585,
                             salts_gfdl_ssp585, 7,
                             larval_formula)
grids_fhslarvae8 <- pred_loop(2050:2054, fhs_larvae, 130, 
                             temps_gfdl_ssp585,
                             salts_gfdl_ssp585, 8,
                             larval_formula)
grids_fhslarvae9 <- pred_loop(2055:2059, fhs_larvae, 130, 
                             temps_gfdl_ssp585,
                             salts_gfdl_ssp585, 9,
                             larval_formula)
grids_fhslarvae10 <- pred_loop(2060:2064, fhs_larvae, 130, 
                              temps_gfdl_ssp585,
                              salts_gfdl_ssp585, 10,
                              larval_formula)
grids_fhslarvae11 <- pred_loop(2065:2069, fhs_larvae, 130, 
                              temps_gfdl_ssp585,
                              salts_gfdl_ssp585, 11,
                              larval_formula)

# Combine into one data frame
df_fhslarvae5 <- list(grids_fhslarvae6[[1]], grids_fhslarvae6[[2]], grids_fhslarvae6[[3]], 
                     grids_fhslarvae6[[4]], grids_fhslarvae6[[5]], grids_fhslarvae7[[1]], 
                     grids_fhslarvae7[[2]], grids_fhslarvae7[[3]], grids_fhslarvae7[[4]],
                     grids_fhslarvae7[[5]], grids_fhslarvae8[[1]], grids_fhslarvae8[[2]], 
                     grids_fhslarvae8[[3]], grids_fhslarvae8[[4]], grids_fhslarvae8[[5]],
                     grids_fhslarvae9[[1]], grids_fhslarvae9[[2]], grids_fhslarvae9[[3]], 
                     grids_fhslarvae9[[4]], grids_fhslarvae9[[5]], grids_fhslarvae10[[1]], 
                     grids_fhslarvae10[[2]], grids_fhslarvae10[[3]], grids_fhslarvae10[[4]],
                     grids_fhslarvae11[[5]], grids_fhslarvae11[[1]], grids_fhslarvae11[[2]],
                     grids_fhslarvae11[[3]], grids_fhslarvae11[[4]], grids_fhslarvae11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhslarvae5), fixed = T)
df_fhslarvae_avg5_gfdl585 <- data.frame(lat = df_fhslarvae5$lat, 
                                       lon = df_fhslarvae5$lon, 
                                       dist = df_fhslarvae5$dist,
                                       avg_pred = rowSums(df_fhslarvae5[, x])/30)
saveRDS(df_fhslarvae_avg5_gfdl585, file = here("data", "df_fhslarvae_avg5_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae_avg5_gfdl585, "Forecasted Distribution 2040 - 2069 \n GFDL SSP585")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_gfdl_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_fhslarvae12 <- pred_loop(2070:2074, fhs_larvae, 130,
                              temps_gfdl_ssp585,
                              salts_gfdl_ssp585, 12,
                              larval_formula)
grids_fhslarvae13 <- pred_loop(2075:2079, fhs_larvae, 130,
                              temps_gfdl_ssp585,
                              salts_gfdl_ssp585, 13,
                              larval_formula)
grids_fhslarvae14 <- pred_loop(2080:2084, fhs_larvae, 130,
                              temps_gfdl_ssp585,
                              salts_gfdl_ssp585, 14,
                              larval_formula)
grids_fhslarvae15 <- pred_loop(2085:2089, fhs_larvae, 130,
                              temps_gfdl_ssp585,
                              salts_gfdl_ssp585, 15,
                              larval_formula)
grids_fhslarvae16 <- pred_loop(2090:2094, fhs_larvae, 130,
                              temps_gfdl_ssp585,
                              salts_gfdl_ssp585, 16,
                              larval_formula)
grids_fhslarvae17 <- pred_loop(2095:2099, fhs_larvae, 130, 
                              temps_gfdl_ssp585,
                              salts_gfdl_ssp585, 17,
                              larval_formula)

# Combine into one data frame
df_fhslarvae6 <- list(grids_fhslarvae12[[1]], grids_fhslarvae12[[2]], grids_fhslarvae12[[3]], 
                     grids_fhslarvae12[[4]], grids_fhslarvae12[[5]], grids_fhslarvae13[[1]], 
                     grids_fhslarvae13[[2]], grids_fhslarvae13[[3]], grids_fhslarvae13[[4]],
                     grids_fhslarvae13[[5]], grids_fhslarvae14[[1]], grids_fhslarvae14[[2]], 
                     grids_fhslarvae14[[3]], grids_fhslarvae14[[4]], grids_fhslarvae14[[5]],
                     grids_fhslarvae15[[1]], grids_fhslarvae15[[2]], grids_fhslarvae15[[3]], 
                     grids_fhslarvae15[[4]], grids_fhslarvae15[[5]], grids_fhslarvae16[[1]], 
                     grids_fhslarvae16[[2]], grids_fhslarvae16[[3]], grids_fhslarvae16[[4]],
                     grids_fhslarvae17[[5]], grids_fhslarvae17[[1]], grids_fhslarvae17[[2]],
                     grids_fhslarvae17[[3]], grids_fhslarvae17[[4]], grids_fhslarvae17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhslarvae6), fixed = T)
df_fhslarvae_avg6_gfdl585 <- data.frame(lat = df_fhslarvae6$lat, 
                                       lon = df_fhslarvae6$lon, 
                                       dist = df_fhslarvae6$dist,
                                       avg_pred = rowSums(df_fhslarvae6[, x])/30)
saveRDS(df_fhslarvae_avg6_gfdl585, file = here("data", "df_fhslarvae_avg6_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae_avg6_gfdl585, "Forecasted Distribution 2070 - 2099 \n GFDL SSP585")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_gfdl_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### MIROC 126 ----------------------------------------------------------------------------------------------------------------------------
temps_miroc_ssp126 <- readRDS(here('data', 'temps_miroc_ssp126.rds'))
salts_miroc_ssp126 <- readRDS(here('data', 'salts_miroc_ssp126.rds'))

## 2015 - 2039
grids_fhslarvae1 <- pred_loop(2015:2019, fhs_larvae, 130, 
                             temps_miroc_ssp126,
                             salts_miroc_ssp126, 1,
                             larval_formula)
grids_fhslarvae2 <- pred_loop(2020:2024, fhs_larvae, 130, 
                             temps_miroc_ssp126,
                             salts_miroc_ssp126, 2,
                             larval_formula)
grids_fhslarvae3 <- pred_loop(2025:2029, fhs_larvae, 130,
                             temps_miroc_ssp126,
                             salts_miroc_ssp126, 3,
                             larval_formula)
grids_fhslarvae4 <- pred_loop(2030:2034, fhs_larvae, 130, 
                             temps_miroc_ssp126,
                             salts_miroc_ssp126, 4,
                             larval_formula)
grids_fhslarvae5 <- pred_loop(2035:2039, fhs_larvae, 130, 
                             temps_miroc_ssp126,
                             salts_miroc_ssp126, 5,
                             larval_formula)

# Combine into one data frame
df_fhslarvae1 <- list(grids_fhslarvae1[[1]], grids_fhslarvae1[[2]], grids_fhslarvae1[[3]], 
                     grids_fhslarvae1[[4]], grids_fhslarvae1[[5]], grids_fhslarvae2[[1]], 
                     grids_fhslarvae2[[2]], grids_fhslarvae2[[3]], grids_fhslarvae2[[4]],
                     grids_fhslarvae2[[5]], grids_fhslarvae3[[1]], grids_fhslarvae3[[2]], 
                     grids_fhslarvae3[[3]], grids_fhslarvae3[[4]], grids_fhslarvae3[[5]],
                     grids_fhslarvae4[[1]], grids_fhslarvae4[[2]], grids_fhslarvae4[[3]], 
                     grids_fhslarvae4[[4]], grids_fhslarvae4[[5]], grids_fhslarvae5[[1]], 
                     grids_fhslarvae5[[2]], grids_fhslarvae5[[3]], grids_fhslarvae5[[4]],
                     grids_fhslarvae5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhslarvae1), fixed = T)
df_fhslarvae_avg1_miroc126 <- data.frame(lat = df_fhslarvae1$lat, 
                                        lon = df_fhslarvae1$lon, 
                                        dist = df_fhslarvae1$dist,
                                        avg_pred = rowSums(df_fhslarvae1[, x])/25)
saveRDS(df_fhslarvae_avg1_miroc126, file = here("data", "df_fhslarvae_avg1_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae_avg1_miroc126, "Forecasted Distribution 2015 - 2039 \n MIROC SSP126")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_miroc_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_fhslarvae6 <- pred_loop(2040:2044, fhs_larvae, 130, 
                             temps_miroc_ssp126,
                             salts_miroc_ssp126, 6,
                             larval_formula)
grids_fhslarvae7 <- pred_loop(2045:2049, fhs_larvae, 130, 
                             temps_miroc_ssp126,
                             salts_miroc_ssp126, 7,
                             larval_formula)
grids_fhslarvae8 <- pred_loop(2050:2054, fhs_larvae, 130, 
                             temps_miroc_ssp126,
                             salts_miroc_ssp126, 8,
                             larval_formula)
grids_fhslarvae9 <- pred_loop(2055:2059, fhs_larvae, 130, 
                             temps_miroc_ssp126,
                             salts_miroc_ssp126, 9,
                             larval_formula)
grids_fhslarvae10 <- pred_loop(2060:2064, fhs_larvae, 130, 
                              temps_miroc_ssp126,
                              salts_miroc_ssp126, 10,
                              larval_formula)
grids_fhslarvae11 <- pred_loop(2065:2069, fhs_larvae, 130, 
                              temps_miroc_ssp126,
                              salts_miroc_ssp126, 11,
                              larval_formula)

# Combine into one data frame
df_fhslarvae2 <- list(grids_fhslarvae6[[1]], grids_fhslarvae6[[2]], grids_fhslarvae6[[3]], 
                     grids_fhslarvae6[[4]], grids_fhslarvae6[[5]], grids_fhslarvae7[[1]], 
                     grids_fhslarvae7[[2]], grids_fhslarvae7[[3]], grids_fhslarvae7[[4]],
                     grids_fhslarvae7[[5]], grids_fhslarvae8[[1]], grids_fhslarvae8[[2]], 
                     grids_fhslarvae8[[3]], grids_fhslarvae8[[4]], grids_fhslarvae8[[5]],
                     grids_fhslarvae9[[1]], grids_fhslarvae9[[2]], grids_fhslarvae9[[3]], 
                     grids_fhslarvae9[[4]], grids_fhslarvae9[[5]], grids_fhslarvae10[[1]], 
                     grids_fhslarvae10[[2]], grids_fhslarvae10[[3]], grids_fhslarvae10[[4]],
                     grids_fhslarvae11[[5]], grids_fhslarvae11[[1]], grids_fhslarvae11[[2]],
                     grids_fhslarvae11[[3]], grids_fhslarvae11[[4]], grids_fhslarvae11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhslarvae2), fixed = T)
df_fhslarvae_avg2_miroc126 <- data.frame(lat = df_fhslarvae2$lat, 
                                        lon = df_fhslarvae2$lon, 
                                        dist = df_fhslarvae2$dist,
                                        avg_pred = rowSums(df_fhslarvae2[, x])/30)
saveRDS(df_fhslarvae_avg2_miroc126, file = here("data", "df_fhslarvae_avg2_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae_avg2_miroc126, "Forecasted Distribution 2040 - 2069 \n MIROC SSP126")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_miroc_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_fhslarvae12 <- pred_loop(2070:2074, fhs_larvae, 130,
                              temps_miroc_ssp126,
                              salts_miroc_ssp126, 12,
                              larval_formula)
grids_fhslarvae13 <- pred_loop(2075:2079, fhs_larvae, 130,
                              temps_miroc_ssp126,
                              salts_miroc_ssp126, 13,
                              larval_formula)
grids_fhslarvae14 <- pred_loop(2080:2084, fhs_larvae, 130,
                              temps_miroc_ssp126,
                              salts_miroc_ssp126, 14,
                              larval_formula)
grids_fhslarvae15 <- pred_loop(2085:2089, fhs_larvae, 130,
                              temps_miroc_ssp126,
                              salts_miroc_ssp126, 15,
                              larval_formula)
grids_fhslarvae16 <- pred_loop(2090:2094, fhs_larvae, 130,
                              temps_miroc_ssp126,
                              salts_miroc_ssp126, 16,
                              larval_formula)
grids_fhslarvae17 <- pred_loop(2095:2099, fhs_larvae, 130, 
                              temps_miroc_ssp126,
                              salts_miroc_ssp126, 17,
                              larval_formula)

# Combine into one data frame
df_fhslarvae3 <- list(grids_fhslarvae12[[1]], grids_fhslarvae12[[2]], grids_fhslarvae12[[3]], 
                     grids_fhslarvae12[[4]], grids_fhslarvae12[[5]], grids_fhslarvae13[[1]], 
                     grids_fhslarvae13[[2]], grids_fhslarvae13[[3]], grids_fhslarvae13[[4]],
                     grids_fhslarvae13[[5]], grids_fhslarvae14[[1]], grids_fhslarvae14[[2]], 
                     grids_fhslarvae14[[3]], grids_fhslarvae14[[4]], grids_fhslarvae14[[5]],
                     grids_fhslarvae15[[1]], grids_fhslarvae15[[2]], grids_fhslarvae15[[3]], 
                     grids_fhslarvae15[[4]], grids_fhslarvae15[[5]], grids_fhslarvae16[[1]], 
                     grids_fhslarvae16[[2]], grids_fhslarvae16[[3]], grids_fhslarvae16[[4]],
                     grids_fhslarvae17[[5]], grids_fhslarvae17[[1]], grids_fhslarvae17[[2]],
                     grids_fhslarvae17[[3]], grids_fhslarvae17[[4]], grids_fhslarvae17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhslarvae3), fixed = T)
df_fhslarvae_avg3_miroc126 <- data.frame(lat = df_fhslarvae3$lat, 
                                        lon = df_fhslarvae3$lon, 
                                        dist = df_fhslarvae3$dist,
                                        avg_pred = rowSums(df_fhslarvae3[, x])/30)
saveRDS(df_fhslarvae_avg3_miroc126, file = here("data", "df_fhslarvae_avg3_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae_avg3_miroc126, "Forecasted Distribution 2070 - 2099 \n MIROC SSP126")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_miroc_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

##### MIROC 585 -----------------------------------------------------------------------------------------------------------------
temps_miroc_ssp585 <- readRDS(here('data', 'temps_miroc_ssp585.rds'))
salts_miroc_ssp585 <- readRDS(here('data', 'salts_miroc_ssp585.rds'))

## 2015 - 2039
grids_fhslarvae1 <- pred_loop(2015:2019, fhs_larvae, 130, 
                             temps_miroc_ssp585,
                             salts_miroc_ssp585, 1,
                             larval_formula)
grids_fhslarvae2 <- pred_loop(2020:2024, fhs_larvae, 130, 
                             temps_miroc_ssp585,
                             salts_miroc_ssp585, 2,
                             larval_formula)
grids_fhslarvae3 <- pred_loop(2025:2029, fhs_larvae, 130,
                             temps_miroc_ssp585,
                             salts_miroc_ssp585, 3,
                             larval_formula)
grids_fhslarvae4 <- pred_loop(2030:2034, fhs_larvae, 130, 
                             temps_miroc_ssp585,
                             salts_miroc_ssp585, 4,
                             larval_formula)
grids_fhslarvae5 <- pred_loop(2035:2039, fhs_larvae, 130, 
                             temps_miroc_ssp585,
                             salts_miroc_ssp585, 5,
                             larval_formula)

# Combine into one data frame
df_fhslarvae4 <- list(grids_fhslarvae1[[1]], grids_fhslarvae1[[2]], grids_fhslarvae1[[3]], 
                     grids_fhslarvae1[[4]], grids_fhslarvae1[[5]], grids_fhslarvae2[[1]], 
                     grids_fhslarvae2[[2]], grids_fhslarvae2[[3]], grids_fhslarvae2[[4]],
                     grids_fhslarvae2[[5]], grids_fhslarvae3[[1]], grids_fhslarvae3[[2]], 
                     grids_fhslarvae3[[3]], grids_fhslarvae3[[4]], grids_fhslarvae3[[5]],
                     grids_fhslarvae4[[1]], grids_fhslarvae4[[2]], grids_fhslarvae4[[3]], 
                     grids_fhslarvae4[[4]], grids_fhslarvae4[[5]], grids_fhslarvae5[[1]], 
                     grids_fhslarvae5[[2]], grids_fhslarvae5[[3]], grids_fhslarvae5[[4]],
                     grids_fhslarvae5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhslarvae4), fixed = T)
df_fhslarvae_avg4_miroc585 <- data.frame(lat = df_fhslarvae4$lat, 
                                        lon = df_fhslarvae4$lon, 
                                        dist = df_fhslarvae4$dist,
                                        avg_pred = rowSums(df_fhslarvae4[, x])/25)
saveRDS(df_fhslarvae_avg4_miroc585, file = here("data", "df_fhslarvae_avg4_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae_avg4_miroc585, "Forecasted Distribution 2015 - 2039 \n MIROC SSP585")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_miroc_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_fhslarvae6 <- pred_loop(2040:2044, fhs_larvae, 130, 
                             temps_miroc_ssp585,
                             salts_miroc_ssp585, 6,
                             larval_formula)
grids_fhslarvae7 <- pred_loop(2045:2049, fhs_larvae, 130, 
                             temps_miroc_ssp585,
                             salts_miroc_ssp585, 7,
                             larval_formula)
grids_fhslarvae8 <- pred_loop(2050:2054, fhs_larvae, 130, 
                             temps_miroc_ssp585,
                             salts_miroc_ssp585, 8,
                             larval_formula)
grids_fhslarvae9 <- pred_loop(2055:2059, fhs_larvae, 130, 
                             temps_miroc_ssp585,
                             salts_miroc_ssp585, 9,
                             larval_formula)
grids_fhslarvae10 <- pred_loop(2060:2064, fhs_larvae, 130, 
                              temps_miroc_ssp585,
                              salts_miroc_ssp585, 10,
                              larval_formula)
grids_fhslarvae11 <- pred_loop(2065:2069, fhs_larvae, 130, 
                              temps_miroc_ssp585,
                              salts_miroc_ssp585, 11,
                              larval_formula)

# Combine into one data frame
df_fhslarvae5 <- list(grids_fhslarvae6[[1]], grids_fhslarvae6[[2]], grids_fhslarvae6[[3]], 
                     grids_fhslarvae6[[4]], grids_fhslarvae6[[5]], grids_fhslarvae7[[1]], 
                     grids_fhslarvae7[[2]], grids_fhslarvae7[[3]], grids_fhslarvae7[[4]],
                     grids_fhslarvae7[[5]], grids_fhslarvae8[[1]], grids_fhslarvae8[[2]], 
                     grids_fhslarvae8[[3]], grids_fhslarvae8[[4]], grids_fhslarvae8[[5]],
                     grids_fhslarvae9[[1]], grids_fhslarvae9[[2]], grids_fhslarvae9[[3]], 
                     grids_fhslarvae9[[4]], grids_fhslarvae9[[5]], grids_fhslarvae10[[1]], 
                     grids_fhslarvae10[[2]], grids_fhslarvae10[[3]], grids_fhslarvae10[[4]],
                     grids_fhslarvae11[[5]], grids_fhslarvae11[[1]], grids_fhslarvae11[[2]],
                     grids_fhslarvae11[[3]], grids_fhslarvae11[[4]], grids_fhslarvae11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhslarvae5), fixed = T)
df_fhslarvae_avg5_miroc585 <- data.frame(lat = df_fhslarvae5$lat, 
                                        lon = df_fhslarvae5$lon, 
                                        dist = df_fhslarvae5$dist,
                                        avg_pred = rowSums(df_fhslarvae5[, x])/30)
saveRDS(df_fhslarvae_avg5_miroc585, file = here("data", "df_fhslarvae_avg5_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae_avg5_miroc585, "Forecasted Distribution 2040 - 2069 \n MIROC SSP585")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_miroc_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_fhslarvae12 <- pred_loop(2070:2074, fhs_larvae, 130,
                              temps_miroc_ssp585,
                              salts_miroc_ssp585, 12,
                              larval_formula)
grids_fhslarvae13 <- pred_loop(2075:2079, fhs_larvae, 130,
                              temps_miroc_ssp585,
                              salts_miroc_ssp585, 13,
                              larval_formula)
grids_fhslarvae14 <- pred_loop(2080:2084, fhs_larvae, 130,
                              temps_miroc_ssp585,
                              salts_miroc_ssp585, 14,
                              larval_formula)
grids_fhslarvae15 <- pred_loop(2085:2089, fhs_larvae, 130,
                              temps_miroc_ssp585,
                              salts_miroc_ssp585, 15,
                              larval_formula)
grids_fhslarvae16 <- pred_loop(2090:2094, fhs_larvae, 130,
                              temps_miroc_ssp585,
                              salts_miroc_ssp585, 16,
                              larval_formula)
grids_fhslarvae17 <- pred_loop(2095:2099, fhs_larvae, 130, 
                              temps_miroc_ssp585,
                              salts_miroc_ssp585, 17,
                              larval_formula)

# Combine into one data frame
df_fhslarvae6 <- list(grids_fhslarvae12[[1]], grids_fhslarvae12[[2]], grids_fhslarvae12[[3]], 
                     grids_fhslarvae12[[4]], grids_fhslarvae12[[5]], grids_fhslarvae13[[1]], 
                     grids_fhslarvae13[[2]], grids_fhslarvae13[[3]], grids_fhslarvae13[[4]],
                     grids_fhslarvae13[[5]], grids_fhslarvae14[[1]], grids_fhslarvae14[[2]], 
                     grids_fhslarvae14[[3]], grids_fhslarvae14[[4]], grids_fhslarvae14[[5]],
                     grids_fhslarvae15[[1]], grids_fhslarvae15[[2]], grids_fhslarvae15[[3]], 
                     grids_fhslarvae15[[4]], grids_fhslarvae15[[5]], grids_fhslarvae16[[1]], 
                     grids_fhslarvae16[[2]], grids_fhslarvae16[[3]], grids_fhslarvae16[[4]],
                     grids_fhslarvae17[[5]], grids_fhslarvae17[[1]], grids_fhslarvae17[[2]],
                     grids_fhslarvae17[[3]], grids_fhslarvae17[[4]], grids_fhslarvae17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhslarvae6), fixed = T)
df_fhslarvae_avg6_miroc585 <- data.frame(lat = df_fhslarvae6$lat, 
                                        lon = df_fhslarvae6$lon, 
                                        dist = df_fhslarvae6$dist,
                                        avg_pred = rowSums(df_fhslarvae6[, x])/30)
saveRDS(df_fhslarvae_avg6_miroc585, file = here("data", "df_fhslarvae_avg6_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae_avg6_miroc585, "Forecasted Distribution 2070 - 2099 \n MIROC SSP585")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_miroc_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

### Average Predictions ------------------------------------------------------------------------------------------------------------


#### 2015-2039 ---------------------------------------------------------------------------------------------------------------------
##### Eggs
df_fhsegg_avg1_cesm126 <- readRDS(here('data', 'df_fhsegg_avg1_cesm126.rds'))
df_fhsegg_avg4_cesm585 <- readRDS(here('data', 'df_fhsegg_avg4_cesm585.rds'))
df_fhsegg_avg1_gfdl126 <- readRDS(here('data', 'df_fhsegg_avg1_gfdl126.rds'))
df_fhsegg_avg4_gfdl585 <- readRDS(here('data', 'df_fhsegg_avg4_gfdl585.rds'))
df_fhsegg_avg1_miroc126 <- readRDS(here('data', 'df_fhsegg_avg1_miroc126.rds'))
df_fhsegg_avg4_miroc585 <- readRDS(here('data', 'df_fhsegg_avg4_miroc585.rds'))

df_fhsegg_merged1 <- list(df_fhsegg_avg1_cesm126, df_fhsegg_avg4_cesm585,
                         df_fhsegg_avg1_gfdl126, df_fhsegg_avg4_gfdl585,
                         df_fhsegg_avg1_miroc126, df_fhsegg_avg4_miroc585) %>%
  reduce(inner_join, by = c("lon", "lat", "dist"))

x <- grepl("pred", names(df_fhsegg_merged1), fixed = T)
df_fhsegg_final1 <- data.frame(lat = df_fhsegg_merged1$lat,
                              lon = df_fhsegg_merged1$lon,
                              avg_pred = (rowSums(df_fhsegg_merged1[, x])/6))

windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg_final1, "Forecasted Distribution 2015 - 2039")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_avg1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

##### Larvae
df_fhslarvae_avg1_cesm126 <- readRDS(here('data', 'df_fhslarvae_avg1_cesm126.rds'))
df_fhslarvae_avg4_cesm585 <- readRDS(here('data', 'df_fhslarvae_avg4_cesm585.rds'))
df_fhslarvae_avg1_gfdl126 <- readRDS(here('data', 'df_fhslarvae_avg1_gfdl126.rds'))
df_fhslarvae_avg4_gfdl585 <- readRDS(here('data', 'df_fhslarvae_avg4_gfdl585.rds'))
df_fhslarvae_avg1_miroc126 <- readRDS(here('data', 'df_fhslarvae_avg1_miroc126.rds'))
df_fhslarvae_avg4_miroc585 <- readRDS(here('data', 'df_fhslarvae_avg4_miroc585.rds'))

df_fhslarvae_merged1 <- list(df_fhslarvae_avg1_cesm126, df_fhslarvae_avg4_cesm585,
                            df_fhslarvae_avg1_gfdl126, df_fhslarvae_avg4_gfdl585,
                            df_fhslarvae_avg1_miroc126, df_fhslarvae_avg4_miroc585) %>%
  reduce(inner_join, by = c("lon", "lat", "dist"))

x <- grepl("pred", names(df_fhslarvae_merged1), fixed = T)
df_fhslarvae_final1 <- data.frame(lat = df_fhslarvae_merged1$lat,
                                 lon = df_fhslarvae_merged1$lon,
                                 avg_pred = (rowSums(df_fhslarvae_merged1[, x])/6))

windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae_final1, "Forecasted Distribution 2015 - 2039")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_avg1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


#### 2040-2069 ---------------------------------------------------------------------------------------------------------------------
##### Eggs
df_fhsegg_avg2_cesm126 <- readRDS(here('data', 'df_fhsegg_avg2_cesm126.rds'))
df_fhsegg_avg5_cesm585 <- readRDS(here('data', 'df_fhsegg_avg5_cesm585.rds'))
df_fhsegg_avg2_gfdl126 <- readRDS(here('data', 'df_fhsegg_avg2_gfdl126.rds'))
df_fhsegg_avg5_gfdl585 <- readRDS(here('data', 'df_fhsegg_avg5_gfdl585.rds'))
df_fhsegg_avg2_miroc126 <- readRDS(here('data', 'df_fhsegg_avg2_miroc126.rds'))
df_fhsegg_avg5_miroc585 <- readRDS(here('data', 'df_fhsegg_avg5_miroc585.rds'))

df_fhsegg_merged2 <- list(df_fhsegg_avg2_cesm126, df_fhsegg_avg5_cesm585,
                         df_fhsegg_avg2_gfdl126, df_fhsegg_avg5_gfdl585,
                         df_fhsegg_avg2_miroc126, df_fhsegg_avg5_miroc585) %>%
  reduce(inner_join, by = c("lon", "lat", "dist"))

x <- grepl("pred", names(df_fhsegg_merged2), fixed = T)
df_fhsegg_final2 <- data.frame(lat = df_fhsegg_merged2$lat,
                              lon = df_fhsegg_merged2$lon,
                              avg_pred = (rowSums(df_fhsegg_merged2[, x])/6))

windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg_final2, "Forecasted Distribution 2040 - 2069")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_avg2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

##### Larvae
df_fhslarvae_avg2_cesm126 <- readRDS(here('data', 'df_fhslarvae_avg2_cesm126.rds'))
df_fhslarvae_avg5_cesm585 <- readRDS(here('data', 'df_fhslarvae_avg5_cesm585.rds'))
df_fhslarvae_avg2_gfdl126 <- readRDS(here('data', 'df_fhslarvae_avg2_gfdl126.rds'))
df_fhslarvae_avg5_gfdl585 <- readRDS(here('data', 'df_fhslarvae_avg5_gfdl585.rds'))
df_fhslarvae_avg2_miroc126 <- readRDS(here('data', 'df_fhslarvae_avg2_miroc126.rds'))
df_fhslarvae_avg5_miroc585 <- readRDS(here('data', 'df_fhslarvae_avg5_miroc585.rds'))

df_fhslarvae_merged2 <- list(df_fhslarvae_avg2_cesm126, df_fhslarvae_avg5_cesm585,
                            df_fhslarvae_avg2_gfdl126, df_fhslarvae_avg5_gfdl585,
                            df_fhslarvae_avg2_miroc126, df_fhslarvae_avg5_miroc585) %>%
  reduce(inner_join, by = c("lon", "lat", "dist"))

x <- grepl("pred", names(df_fhslarvae_merged2), fixed = T)
df_fhslarvae_final2 <- data.frame(lat = df_fhslarvae_merged2$lat,
                                 lon = df_fhslarvae_merged2$lon,
                                 avg_pred = (rowSums(df_fhslarvae_merged2[, x])/6))

windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae_final2, "Forecasted Distribution 2040 - 2069")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_avg2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


#### 2070-2099----------------------------------------------------------------------------------------------------------------------
df_fhsegg_avg3_cesm126 <- readRDS(here('data', 'df_fhsegg_avg3_cesm126.rds'))
df_fhsegg_avg6_cesm585 <- readRDS(here('data', 'df_fhsegg_avg6_cesm585.rds'))
df_fhsegg_avg3_gfdl126 <- readRDS(here('data', 'df_fhsegg_avg3_gfdl126.rds'))
df_fhsegg_avg6_gfdl585 <- readRDS(here('data', 'df_fhsegg_avg6_gfdl585.rds'))
df_fhsegg_avg3_miroc126 <- readRDS(here('data', 'df_fhsegg_avg3_miroc126.rds'))
df_fhsegg_avg6_miroc585 <- readRDS(here('data', 'df_fhsegg_avg6_miroc585.rds'))

df_fhsegg_merged3 <- list(df_fhsegg_avg3_cesm126, df_fhsegg_avg6_cesm585,
                         df_fhsegg_avg3_gfdl126, df_fhsegg_avg6_gfdl585,
                         df_fhsegg_avg3_miroc126, df_fhsegg_avg6_miroc585) %>%
  reduce(inner_join, by = c("lon", "lat", "dist"))

x <- grepl("pred", names(df_fhsegg_merged3), fixed = T)
df_fhsegg_final3 <- data.frame(lat = df_fhsegg_merged3$lat,
                              lon = df_fhsegg_merged3$lon,
                              avg_pred = (rowSums(df_fhsegg_merged3[, x])/6))

windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg_final3, "Forecasted Distribution 2070 - 2099")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_avg3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

df_fhslarvae_avg3_cesm126 <- readRDS(here('data', 'df_fhslarvae_avg3_cesm126.rds'))
df_fhslarvae_avg6_cesm585 <- readRDS(here('data', 'df_fhslarvae_avg6_cesm585.rds'))
df_fhslarvae_avg3_gfdl126 <- readRDS(here('data', 'df_fhslarvae_avg3_gfdl126.rds'))
df_fhslarvae_avg6_gfdl585 <- readRDS(here('data', 'df_fhslarvae_avg6_gfdl585.rds'))
df_fhslarvae_avg3_miroc126 <- readRDS(here('data', 'df_fhslarvae_avg3_miroc126.rds'))
df_fhslarvae_avg6_miroc585 <- readRDS(here('data', 'df_fhslarvae_avg6_miroc585.rds'))

df_fhslarvae_merged3 <- list(df_fhslarvae_avg3_cesm126, df_fhslarvae_avg6_cesm585,
                            df_fhslarvae_avg3_gfdl126, df_fhslarvae_avg6_gfdl585,
                            df_fhslarvae_avg3_miroc126, df_fhslarvae_avg6_miroc585) %>%
  reduce(inner_join, by = c("lon", "lat", "dist"))

x <- grepl("pred", names(df_fhslarvae_merged3), fixed = T)
df_fhslarvae_final3 <- data.frame(lat = df_fhslarvae_merged3$lat,
                                 lon = df_fhslarvae_merged3$lon,
                                 avg_pred = (rowSums(df_fhslarvae_merged3[, x])/6))

windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae_final3, "Forecasted Distribution 2070 - 2099")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_avg3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()