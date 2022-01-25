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
nrs_larvae <- as.data.frame(filter(readRDS(here('data', 'nrs_larvae.rds')),
                                    lat >= 52 & lat <= 62,
                                    lon >= -176.5 & lon <= -156.5))
nrs_larvae$mean_temp <- roms_temps$mean[match(nrs_larvae$year, roms_temps$year)]
nrs_larvae$catch <- nrs_larvae$larvalcatchper10m2 + 1

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

larval_formula <- gam(catch ~ factor(year) + 
                        s(doy, k = 8) +
                        s(lon, lat) +
                        s(roms_temperature, k = 6) +
                        s(roms_salinity, k = 6) +
                        s(lat, lon, by = mean_temp, k = 6),
                      data = nrs_larvae,
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
                      list, formula, year){
  grids <- list()
  for(j in range) {
    date1 <- paste(j, "-05-10", sep = "")
    date2 <- paste(j, "-02-01", sep = "")
    date3 <- paste(j, "-04-30", sep = "")
    grid <- get_preds(data, year, date1, doy,
                      date2, date3,
                      temp_output, salt_output,
                      list, formula)
    grids[[paste("year", j, sep = "")]] <- grid
  }
  return(grids)
}


### Northern Rock Sole Larvae --------------------------------------------------------------------------------------------------------------------------
#### Forecast and average into 3 time periods ---------------------------------------------------------------------------------------------
##### CESM 126 ----------------------------------------------------------------------------------------------------------------------------
temps_cesm_ssp126 <- readRDS(here('data', 'temps_cesm_ssp126.rds'))
salts_cesm_ssp126 <- readRDS(here('data', 'salts_cesm_ssp126.rds'))

## 2015 - 2039
grids_nrslarvae1 <- pred_loop(2015:2019, nrs_larvae, 130, 
                               temps_cesm_ssp126,
                               salts_cesm_ssp126, 1,
                               larval_formula, 1994)
grids_nrslarvae2 <- pred_loop(2020:2024, nrs_larvae, 130, 
                               temps_cesm_ssp126,
                               salts_cesm_ssp126, 2,
                               larval_formula, 1994)
grids_nrslarvae3 <- pred_loop(2025:2029, nrs_larvae, 130,
                               temps_cesm_ssp126,
                               salts_cesm_ssp126, 3,
                               larval_formula, 1994)
grids_nrslarvae4 <- pred_loop(2030:2034, nrs_larvae, 130, 
                               temps_cesm_ssp126,
                               salts_cesm_ssp126, 4,
                               larval_formula, 1994)
grids_nrslarvae5 <- pred_loop(2035:2039, nrs_larvae, 130, 
                               temps_cesm_ssp126,
                               salts_cesm_ssp126, 5,
                               larval_formula, 1994)

# Combine into one data frame
df_nrslarvae1 <- list(grids_nrslarvae1[[1]], grids_nrslarvae1[[2]], grids_nrslarvae1[[3]], 
                       grids_nrslarvae1[[4]], grids_nrslarvae1[[5]], grids_nrslarvae2[[1]], 
                       grids_nrslarvae2[[2]], grids_nrslarvae2[[3]], grids_nrslarvae2[[4]],
                       grids_nrslarvae2[[5]], grids_nrslarvae3[[1]], grids_nrslarvae3[[2]], 
                       grids_nrslarvae3[[3]], grids_nrslarvae3[[4]], grids_nrslarvae3[[5]],
                       grids_nrslarvae4[[1]], grids_nrslarvae4[[2]], grids_nrslarvae4[[3]], 
                       grids_nrslarvae4[[4]], grids_nrslarvae4[[5]], grids_nrslarvae5[[1]], 
                       grids_nrslarvae5[[2]], grids_nrslarvae5[[3]], grids_nrslarvae5[[4]],
                       grids_nrslarvae5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_nrslarvae1), fixed = T)
df_nrslarvae_avg1_cesm126 <- data.frame(lat = df_nrslarvae1$lat, 
                                         lon = df_nrslarvae1$lon, 
                                         dist = df_nrslarvae1$dist,
                                         avg_pred = rowSums(df_nrslarvae1[, x])/25)
saveRDS(df_nrslarvae_avg1_cesm126, file = here("data", "df_nrslarvae_avg1_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae_avg1_cesm126, "Forecasted Distribution 2015 - 2039 \n CESM SSP126")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_cesm_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_nrslarvae6 <- pred_loop(2040:2044, nrs_larvae, 130, 
                               temps_cesm_ssp126,
                               salts_cesm_ssp126, 6,
                               larval_formula, 1994)
grids_nrslarvae7 <- pred_loop(2045:2049, nrs_larvae, 130, 
                               temps_cesm_ssp126,
                               salts_cesm_ssp126, 7,
                               larval_formula, 1994)
grids_nrslarvae8 <- pred_loop(2050:2054, nrs_larvae, 130, 
                               temps_cesm_ssp126,
                               salts_cesm_ssp126, 8,
                               larval_formula, 1994)
grids_nrslarvae9 <- pred_loop(2055:2059, nrs_larvae, 130, 
                               temps_cesm_ssp126,
                               salts_cesm_ssp126, 9,
                               larval_formula, 1994)
grids_nrslarvae10 <- pred_loop(2060:2064, nrs_larvae, 130, 
                                temps_cesm_ssp126,
                                salts_cesm_ssp126, 10,
                                larval_formula, 1994)
grids_nrslarvae11 <- pred_loop(2065:2069, nrs_larvae, 130, 
                                temps_cesm_ssp126,
                                salts_cesm_ssp126, 11,
                                larval_formula, 1994)

# Combine into one data frame
df_nrslarvae2 <- list(grids_nrslarvae6[[1]], grids_nrslarvae6[[2]], grids_nrslarvae6[[3]], 
                       grids_nrslarvae6[[4]], grids_nrslarvae6[[5]], grids_nrslarvae7[[1]], 
                       grids_nrslarvae7[[2]], grids_nrslarvae7[[3]], grids_nrslarvae7[[4]],
                       grids_nrslarvae7[[5]], grids_nrslarvae8[[1]], grids_nrslarvae8[[2]], 
                       grids_nrslarvae8[[3]], grids_nrslarvae8[[4]], grids_nrslarvae8[[5]],
                       grids_nrslarvae9[[1]], grids_nrslarvae9[[2]], grids_nrslarvae9[[3]], 
                       grids_nrslarvae9[[4]], grids_nrslarvae9[[5]], grids_nrslarvae10[[1]], 
                       grids_nrslarvae10[[2]], grids_nrslarvae10[[3]], grids_nrslarvae10[[4]],
                       grids_nrslarvae11[[5]], grids_nrslarvae11[[1]], grids_nrslarvae11[[2]],
                       grids_nrslarvae11[[3]], grids_nrslarvae11[[4]], grids_nrslarvae11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_nrslarvae2), fixed = T)
df_nrslarvae_avg2_cesm126 <- data.frame(lat = df_nrslarvae2$lat, 
                                         lon = df_nrslarvae2$lon, 
                                         dist = df_nrslarvae2$dist,
                                         avg_pred = rowSums(df_nrslarvae2[, x])/30)
saveRDS(df_nrslarvae_avg2_cesm126, file = here("data", "df_nrslarvae_avg2_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae_avg2_cesm126, "Forecasted Distribution 2040 - 2069 \n CESM SSP126")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_cesm_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_nrslarvae12 <- pred_loop(2070:2074, nrs_larvae, 130,
                                temps_cesm_ssp126,
                                salts_cesm_ssp126, 12,
                                larval_formula, 1994)
grids_nrslarvae13 <- pred_loop(2075:2079, nrs_larvae, 130,
                                temps_cesm_ssp126,
                                salts_cesm_ssp126, 13,
                                larval_formula, 1994)
grids_nrslarvae14 <- pred_loop(2080:2084, nrs_larvae, 130,
                                temps_cesm_ssp126,
                                salts_cesm_ssp126, 14,
                                larval_formula, 1994)
grids_nrslarvae15 <- pred_loop(2085:2089, nrs_larvae, 130,
                                temps_cesm_ssp126,
                                salts_cesm_ssp126, 15,
                                larval_formula, 1994)
grids_nrslarvae16 <- pred_loop(2090:2094, nrs_larvae, 130,
                                temps_cesm_ssp126,
                                salts_cesm_ssp126, 16,
                                larval_formula, 1994)
grids_nrslarvae17 <- pred_loop(2095:2099, nrs_larvae, 130, 
                                temps_cesm_ssp126,
                                salts_cesm_ssp126, 17,
                                larval_formula, 1994)

# Combine into one data frame
df_nrslarvae3 <- list(grids_nrslarvae12[[1]], grids_nrslarvae12[[2]], grids_nrslarvae12[[3]], 
                       grids_nrslarvae12[[4]], grids_nrslarvae12[[5]], grids_nrslarvae13[[1]], 
                       grids_nrslarvae13[[2]], grids_nrslarvae13[[3]], grids_nrslarvae13[[4]],
                       grids_nrslarvae13[[5]], grids_nrslarvae14[[1]], grids_nrslarvae14[[2]], 
                       grids_nrslarvae14[[3]], grids_nrslarvae14[[4]], grids_nrslarvae14[[5]],
                       grids_nrslarvae15[[1]], grids_nrslarvae15[[2]], grids_nrslarvae15[[3]], 
                       grids_nrslarvae15[[4]], grids_nrslarvae15[[5]], grids_nrslarvae16[[1]], 
                       grids_nrslarvae16[[2]], grids_nrslarvae16[[3]], grids_nrslarvae16[[4]],
                       grids_nrslarvae17[[5]], grids_nrslarvae17[[1]], grids_nrslarvae17[[2]],
                       grids_nrslarvae17[[3]], grids_nrslarvae17[[4]], grids_nrslarvae17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_nrslarvae3), fixed = T)
df_nrslarvae_avg3_cesm126 <- data.frame(lat = df_nrslarvae3$lat, 
                                         lon = df_nrslarvae3$lon, 
                                         dist = df_nrslarvae3$dist,
                                         avg_pred = rowSums(df_nrslarvae3[, x])/30)
saveRDS(df_nrslarvae_avg3_cesm126, file = here("data", "df_nrslarvae_avg3_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae_avg3_cesm126, "Forecasted Distribution 2070 - 2099 \n CESM SSP126")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_cesm_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

rm(temps_cesm_ssp126)
rm(salts_cesm_ssp126)

##### CESM 585 -----------------------------------------------------------------------------------------------------------------
temps_cesm_ssp585 <- readRDS(here('data', 'temps_cesm_ssp585.rds'))
salts_cesm_ssp585 <- readRDS(here('data', 'salts_cesm_ssp585.rds'))

## 2015 - 2039
grids_nrslarvae1 <- pred_loop(2015:2019, nrs_larvae, 130, 
                               temps_cesm_ssp585,
                               salts_cesm_ssp585, 1,
                               larval_formula, 1994)
grids_nrslarvae2 <- pred_loop(2020:2024, nrs_larvae, 130, 
                               temps_cesm_ssp585,
                               salts_cesm_ssp585, 2,
                               larval_formula, 1994)
grids_nrslarvae3 <- pred_loop(2025:2029, nrs_larvae, 130,
                               temps_cesm_ssp585,
                               salts_cesm_ssp585, 3,
                               larval_formula, 1994)
grids_nrslarvae4 <- pred_loop(2030:2034, nrs_larvae, 130, 
                               temps_cesm_ssp585,
                               salts_cesm_ssp585, 4,
                               larval_formula, 1994)
grids_nrslarvae5 <- pred_loop(2035:2039, nrs_larvae, 130, 
                               temps_cesm_ssp585,
                               salts_cesm_ssp585, 5,
                               larval_formula, 1994)

# Combine into one data frame
df_nrslarvae4 <- list(grids_nrslarvae1[[1]], grids_nrslarvae1[[2]], grids_nrslarvae1[[3]], 
                       grids_nrslarvae1[[4]], grids_nrslarvae1[[5]], grids_nrslarvae2[[1]], 
                       grids_nrslarvae2[[2]], grids_nrslarvae2[[3]], grids_nrslarvae2[[4]],
                       grids_nrslarvae2[[5]], grids_nrslarvae3[[1]], grids_nrslarvae3[[2]], 
                       grids_nrslarvae3[[3]], grids_nrslarvae3[[4]], grids_nrslarvae3[[5]],
                       grids_nrslarvae4[[1]], grids_nrslarvae4[[2]], grids_nrslarvae4[[3]], 
                       grids_nrslarvae4[[4]], grids_nrslarvae4[[5]], grids_nrslarvae5[[1]], 
                       grids_nrslarvae5[[2]], grids_nrslarvae5[[3]], grids_nrslarvae5[[4]],
                       grids_nrslarvae5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_nrslarvae4), fixed = T)
df_nrslarvae_avg4_cesm585 <- data.frame(lat = df_nrslarvae4$lat, 
                                         lon = df_nrslarvae4$lon, 
                                         dist = df_nrslarvae4$dist,
                                         avg_pred = rowSums(df_nrslarvae4[, x])/25)
saveRDS(df_nrslarvae_avg4_cesm585, file = here("data", "df_nrslarvae_avg4_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae_avg4_cesm585, "Forecasted Distribution 2015 - 2039 \n CESM SSP585")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_cesm_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_nrslarvae6 <- pred_loop(2040:2044, nrs_larvae, 130, 
                               temps_cesm_ssp585,
                               salts_cesm_ssp585, 6,
                               larval_formula, 1994)
grids_nrslarvae7 <- pred_loop(2045:2049, nrs_larvae, 130, 
                               temps_cesm_ssp585,
                               salts_cesm_ssp585, 7,
                               larval_formula, 1994)
grids_nrslarvae8 <- pred_loop(2050:2054, nrs_larvae, 130, 
                               temps_cesm_ssp585,
                               salts_cesm_ssp585, 8,
                               larval_formula, 1994)
grids_nrslarvae9 <- pred_loop(2055:2059, nrs_larvae, 130, 
                               temps_cesm_ssp585,
                               salts_cesm_ssp585, 9,
                               larval_formula, 1994)
grids_nrslarvae10 <- pred_loop(2060:2064, nrs_larvae, 130, 
                                temps_cesm_ssp585,
                                salts_cesm_ssp585, 10,
                                larval_formula, 1994)
grids_nrslarvae11 <- pred_loop(2065:2069, nrs_larvae, 130, 
                                temps_cesm_ssp585,
                                salts_cesm_ssp585, 11,
                                larval_formula, 1994)

# Combine into one data frame
df_nrslarvae5 <- list(grids_nrslarvae6[[1]], grids_nrslarvae6[[2]], grids_nrslarvae6[[3]], 
                       grids_nrslarvae6[[4]], grids_nrslarvae6[[5]], grids_nrslarvae7[[1]], 
                       grids_nrslarvae7[[2]], grids_nrslarvae7[[3]], grids_nrslarvae7[[4]],
                       grids_nrslarvae7[[5]], grids_nrslarvae8[[1]], grids_nrslarvae8[[2]], 
                       grids_nrslarvae8[[3]], grids_nrslarvae8[[4]], grids_nrslarvae8[[5]],
                       grids_nrslarvae9[[1]], grids_nrslarvae9[[2]], grids_nrslarvae9[[3]], 
                       grids_nrslarvae9[[4]], grids_nrslarvae9[[5]], grids_nrslarvae10[[1]], 
                       grids_nrslarvae10[[2]], grids_nrslarvae10[[3]], grids_nrslarvae10[[4]],
                       grids_nrslarvae11[[5]], grids_nrslarvae11[[1]], grids_nrslarvae11[[2]],
                       grids_nrslarvae11[[3]], grids_nrslarvae11[[4]], grids_nrslarvae11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_nrslarvae5), fixed = T)
df_nrslarvae_avg5_cesm585 <- data.frame(lat = df_nrslarvae5$lat, 
                                         lon = df_nrslarvae5$lon, 
                                         dist = df_nrslarvae5$dist,
                                         avg_pred = rowSums(df_nrslarvae5[, x])/30)
saveRDS(df_nrslarvae_avg5_cesm585, file = here("data", "df_nrslarvae_avg5_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae_avg5_cesm585, "Forecasted Distribution 2040 - 2069 \n CESM SSP585")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_cesm_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_nrslarvae12 <- pred_loop(2070:2074, nrs_larvae, 130,
                                temps_cesm_ssp585,
                                salts_cesm_ssp585, 12,
                                larval_formula, 1994)
grids_nrslarvae13 <- pred_loop(2075:2079, nrs_larvae, 130,
                                temps_cesm_ssp585,
                                salts_cesm_ssp585, 13,
                                larval_formula, 1994)
grids_nrslarvae14 <- pred_loop(2080:2084, nrs_larvae, 130,
                                temps_cesm_ssp585,
                                salts_cesm_ssp585, 14,
                                larval_formula, 1994)
grids_nrslarvae15 <- pred_loop(2085:2089, nrs_larvae, 130,
                                temps_cesm_ssp585,
                                salts_cesm_ssp585, 15,
                                larval_formula, 1994)
grids_nrslarvae16 <- pred_loop(2090:2094, nrs_larvae, 130,
                                temps_cesm_ssp585,
                                salts_cesm_ssp585, 16,
                                larval_formula, 1994)
grids_nrslarvae17 <- pred_loop(2095:2099, nrs_larvae, 130, 
                                temps_cesm_ssp585,
                                salts_cesm_ssp585, 17,
                                larval_formula, 1994)

# Combine into one data frame
df_nrslarvae6 <- list(grids_nrslarvae12[[1]], grids_nrslarvae12[[2]], grids_nrslarvae12[[3]], 
                       grids_nrslarvae12[[4]], grids_nrslarvae12[[5]], grids_nrslarvae13[[1]], 
                       grids_nrslarvae13[[2]], grids_nrslarvae13[[3]], grids_nrslarvae13[[4]],
                       grids_nrslarvae13[[5]], grids_nrslarvae14[[1]], grids_nrslarvae14[[2]], 
                       grids_nrslarvae14[[3]], grids_nrslarvae14[[4]], grids_nrslarvae14[[5]],
                       grids_nrslarvae15[[1]], grids_nrslarvae15[[2]], grids_nrslarvae15[[3]], 
                       grids_nrslarvae15[[4]], grids_nrslarvae15[[5]], grids_nrslarvae16[[1]], 
                       grids_nrslarvae16[[2]], grids_nrslarvae16[[3]], grids_nrslarvae16[[4]],
                       grids_nrslarvae17[[5]], grids_nrslarvae17[[1]], grids_nrslarvae17[[2]],
                       grids_nrslarvae17[[3]], grids_nrslarvae17[[4]], grids_nrslarvae17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_nrslarvae6), fixed = T)
df_nrslarvae_avg6_cesm585 <- data.frame(lat = df_nrslarvae6$lat, 
                                         lon = df_nrslarvae6$lon, 
                                         dist = df_nrslarvae6$dist,
                                         avg_pred = rowSums(df_nrslarvae6[, x])/30)
saveRDS(df_nrslarvae_avg6_cesm585, file = here("data", "df_nrslarvae_avg6_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae_avg6_cesm585, "Forecasted Distribution 2070 - 2099 \n CESM SSP585")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_cesm_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

rm(temps_cesm_ssp585)
rm(salts_cesm_ssp585)


##### GFDL 126 ----------------------------------------------------------------------------------------------------------------------------
temps_gfdl_ssp126 <- readRDS(here('data', 'temps_gfdl_ssp126.rds'))
salts_gfdl_ssp126 <- readRDS(here('data', 'salts_gfdl_ssp126.rds'))

## 2015 - 2039
grids_nrslarvae1 <- pred_loop(2015:2019, nrs_larvae, 130, 
                               temps_gfdl_ssp126,
                               salts_gfdl_ssp126, 1,
                               larval_formula, 1994)
grids_nrslarvae2 <- pred_loop(2020:2024, nrs_larvae, 130, 
                               temps_gfdl_ssp126,
                               salts_gfdl_ssp126, 2,
                               larval_formula, 1994)
grids_nrslarvae3 <- pred_loop(2025:2029, nrs_larvae, 130,
                               temps_gfdl_ssp126,
                               salts_gfdl_ssp126, 3,
                               larval_formula, 1994)
grids_nrslarvae4 <- pred_loop(2030:2034, nrs_larvae, 130, 
                               temps_gfdl_ssp126,
                               salts_gfdl_ssp126, 4,
                               larval_formula, 1994)
grids_nrslarvae5 <- pred_loop(2035:2039, nrs_larvae, 130, 
                               temps_gfdl_ssp126,
                               salts_gfdl_ssp126, 5,
                               larval_formula, 1994)

# Combine into one data frame
df_nrslarvae1 <- list(grids_nrslarvae1[[1]], grids_nrslarvae1[[2]], grids_nrslarvae1[[3]], 
                       grids_nrslarvae1[[4]], grids_nrslarvae1[[5]], grids_nrslarvae2[[1]], 
                       grids_nrslarvae2[[2]], grids_nrslarvae2[[3]], grids_nrslarvae2[[4]],
                       grids_nrslarvae2[[5]], grids_nrslarvae3[[1]], grids_nrslarvae3[[2]], 
                       grids_nrslarvae3[[3]], grids_nrslarvae3[[4]], grids_nrslarvae3[[5]],
                       grids_nrslarvae4[[1]], grids_nrslarvae4[[2]], grids_nrslarvae4[[3]], 
                       grids_nrslarvae4[[4]], grids_nrslarvae4[[5]], grids_nrslarvae5[[1]], 
                       grids_nrslarvae5[[2]], grids_nrslarvae5[[3]], grids_nrslarvae5[[4]],
                       grids_nrslarvae5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_nrslarvae1), fixed = T)
df_nrslarvae_avg1_gfdl126 <- data.frame(lat = df_nrslarvae1$lat, 
                                         lon = df_nrslarvae1$lon, 
                                         dist = df_nrslarvae1$dist,
                                         avg_pred = rowSums(df_nrslarvae1[, x])/25)
saveRDS(df_nrslarvae_avg1_gfdl126, file = here("data", "df_nrslarvae_avg1_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae_avg1_gfdl126, "Forecasted Distribution 2015 - 2039 \n GFDL SSP126")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_gfdl_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_nrslarvae6 <- pred_loop(2040:2044, nrs_larvae, 130, 
                               temps_gfdl_ssp126,
                               salts_gfdl_ssp126, 6,
                               larval_formula, 1994)
grids_nrslarvae7 <- pred_loop(2045:2049, nrs_larvae, 130, 
                               temps_gfdl_ssp126,
                               salts_gfdl_ssp126, 7,
                               larval_formula, 1994)
grids_nrslarvae8 <- pred_loop(2050:2054, nrs_larvae, 130, 
                               temps_gfdl_ssp126,
                               salts_gfdl_ssp126, 8,
                               larval_formula, 1994)
grids_nrslarvae9 <- pred_loop(2055:2059, nrs_larvae, 130, 
                               temps_gfdl_ssp126,
                               salts_gfdl_ssp126, 9,
                               larval_formula, 1994)
grids_nrslarvae10 <- pred_loop(2060:2064, nrs_larvae, 130, 
                                temps_gfdl_ssp126,
                                salts_gfdl_ssp126, 10,
                                larval_formula, 1994)
grids_nrslarvae11 <- pred_loop(2065:2069, nrs_larvae, 130, 
                                temps_gfdl_ssp126,
                                salts_gfdl_ssp126, 11,
                                larval_formula, 1994)

# Combine into one data frame
df_nrslarvae2 <- list(grids_nrslarvae6[[1]], grids_nrslarvae6[[2]], grids_nrslarvae6[[3]], 
                       grids_nrslarvae6[[4]], grids_nrslarvae6[[5]], grids_nrslarvae7[[1]], 
                       grids_nrslarvae7[[2]], grids_nrslarvae7[[3]], grids_nrslarvae7[[4]],
                       grids_nrslarvae7[[5]], grids_nrslarvae8[[1]], grids_nrslarvae8[[2]], 
                       grids_nrslarvae8[[3]], grids_nrslarvae8[[4]], grids_nrslarvae8[[5]],
                       grids_nrslarvae9[[1]], grids_nrslarvae9[[2]], grids_nrslarvae9[[3]], 
                       grids_nrslarvae9[[4]], grids_nrslarvae9[[5]], grids_nrslarvae10[[1]], 
                       grids_nrslarvae10[[2]], grids_nrslarvae10[[3]], grids_nrslarvae10[[4]],
                       grids_nrslarvae11[[5]], grids_nrslarvae11[[1]], grids_nrslarvae11[[2]],
                       grids_nrslarvae11[[3]], grids_nrslarvae11[[4]], grids_nrslarvae11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_nrslarvae2), fixed = T)
df_nrslarvae_avg2_gfdl126 <- data.frame(lat = df_nrslarvae2$lat, 
                                         lon = df_nrslarvae2$lon, 
                                         dist = df_nrslarvae2$dist,
                                         avg_pred = rowSums(df_nrslarvae2[, x])/30)
saveRDS(df_nrslarvae_avg2_gfdl126, file = here("data", "df_nrslarvae_avg2_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae_avg2_gfdl126, "Forecasted Distribution 2040 - 2069 \n GFDL SSP126")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_gfdl_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_nrslarvae12 <- pred_loop(2070:2074, nrs_larvae, 130,
                                temps_gfdl_ssp126,
                                salts_gfdl_ssp126, 12,
                                larval_formula, 1994)
grids_nrslarvae13 <- pred_loop(2075:2079, nrs_larvae, 130,
                                temps_gfdl_ssp126,
                                salts_gfdl_ssp126, 13,
                                larval_formula, 1994)
grids_nrslarvae14 <- pred_loop(2080:2084, nrs_larvae, 130,
                                temps_gfdl_ssp126,
                                salts_gfdl_ssp126, 14,
                                larval_formula, 1994)
grids_nrslarvae15 <- pred_loop(2085:2089, nrs_larvae, 130,
                                temps_gfdl_ssp126,
                                salts_gfdl_ssp126, 15,
                                larval_formula, 1994)
grids_nrslarvae16 <- pred_loop(2090:2094, nrs_larvae, 130,
                                temps_gfdl_ssp126,
                                salts_gfdl_ssp126, 16,
                                larval_formula, 1994)
grids_nrslarvae17 <- pred_loop(2095:2099, nrs_larvae, 130, 
                                temps_gfdl_ssp126,
                                salts_gfdl_ssp126, 17,
                                larval_formula, 1994)

# Combine into one data frame
df_nrslarvae3 <- list(grids_nrslarvae12[[1]], grids_nrslarvae12[[2]], grids_nrslarvae12[[3]], 
                       grids_nrslarvae12[[4]], grids_nrslarvae12[[5]], grids_nrslarvae13[[1]], 
                       grids_nrslarvae13[[2]], grids_nrslarvae13[[3]], grids_nrslarvae13[[4]],
                       grids_nrslarvae13[[5]], grids_nrslarvae14[[1]], grids_nrslarvae14[[2]], 
                       grids_nrslarvae14[[3]], grids_nrslarvae14[[4]], grids_nrslarvae14[[5]],
                       grids_nrslarvae15[[1]], grids_nrslarvae15[[2]], grids_nrslarvae15[[3]], 
                       grids_nrslarvae15[[4]], grids_nrslarvae15[[5]], grids_nrslarvae16[[1]], 
                       grids_nrslarvae16[[2]], grids_nrslarvae16[[3]], grids_nrslarvae16[[4]],
                       grids_nrslarvae17[[5]], grids_nrslarvae17[[1]], grids_nrslarvae17[[2]],
                       grids_nrslarvae17[[3]], grids_nrslarvae17[[4]], grids_nrslarvae17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_nrslarvae3), fixed = T)
df_nrslarvae_avg3_gfdl126 <- data.frame(lat = df_nrslarvae3$lat, 
                                         lon = df_nrslarvae3$lon, 
                                         dist = df_nrslarvae3$dist,
                                         avg_pred = rowSums(df_nrslarvae3[, x])/30)
saveRDS(df_nrslarvae_avg3_gfdl126, file = here("data", "df_nrslarvae_avg3_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae_avg3_gfdl126, "Forecasted Distribution 2070 - 2099 \n GFDL SSP126")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_gfdl_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

rm(temps_gfdl_ssp126)
rm(salts_gfdl_ssp126)

##### GFDL 585 -----------------------------------------------------------------------------------------------------------------
temps_gfdl_ssp585 <- readRDS(here('data', 'temps_gfdl_ssp585.rds'))
salts_gfdl_ssp585 <- readRDS(here('data', 'salts_gfdl_ssp585.rds'))

## 2015 - 2039
grids_nrslarvae1 <- pred_loop(2015:2019, nrs_larvae, 130, 
                               temps_gfdl_ssp585,
                               salts_gfdl_ssp585, 1,
                               larval_formula, 1994)
grids_nrslarvae2 <- pred_loop(2020:2024, nrs_larvae, 130, 
                               temps_gfdl_ssp585,
                               salts_gfdl_ssp585, 2,
                               larval_formula, 1994)
grids_nrslarvae3 <- pred_loop(2025:2029, nrs_larvae, 130,
                               temps_gfdl_ssp585,
                               salts_gfdl_ssp585, 3,
                               larval_formula, 1994)
grids_nrslarvae4 <- pred_loop(2030:2034, nrs_larvae, 130, 
                               temps_gfdl_ssp585,
                               salts_gfdl_ssp585, 4,
                               larval_formula, 1994)
grids_nrslarvae5 <- pred_loop(2035:2039, nrs_larvae, 130, 
                               temps_gfdl_ssp585,
                               salts_gfdl_ssp585, 5,
                               larval_formula, 1994)

# Combine into one data frame
df_nrslarvae4 <- list(grids_nrslarvae1[[1]], grids_nrslarvae1[[2]], grids_nrslarvae1[[3]], 
                       grids_nrslarvae1[[4]], grids_nrslarvae1[[5]], grids_nrslarvae2[[1]], 
                       grids_nrslarvae2[[2]], grids_nrslarvae2[[3]], grids_nrslarvae2[[4]],
                       grids_nrslarvae2[[5]], grids_nrslarvae3[[1]], grids_nrslarvae3[[2]], 
                       grids_nrslarvae3[[3]], grids_nrslarvae3[[4]], grids_nrslarvae3[[5]],
                       grids_nrslarvae4[[1]], grids_nrslarvae4[[2]], grids_nrslarvae4[[3]], 
                       grids_nrslarvae4[[4]], grids_nrslarvae4[[5]], grids_nrslarvae5[[1]], 
                       grids_nrslarvae5[[2]], grids_nrslarvae5[[3]], grids_nrslarvae5[[4]],
                       grids_nrslarvae5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_nrslarvae4), fixed = T)
df_nrslarvae_avg4_gfdl585 <- data.frame(lat = df_nrslarvae4$lat, 
                                         lon = df_nrslarvae4$lon, 
                                         dist = df_nrslarvae4$dist,
                                         avg_pred = rowSums(df_nrslarvae4[, x])/25)
saveRDS(df_nrslarvae_avg4_gfdl585, file = here("data", "df_nrslarvae_avg4_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae_avg4_gfdl585, "Forecasted Distribution 2015 - 2039 \n GFDL SSP585")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_gfdl_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_nrslarvae6 <- pred_loop(2040:2044, nrs_larvae, 130, 
                               temps_gfdl_ssp585,
                               salts_gfdl_ssp585, 6,
                               larval_formula, 1994)
grids_nrslarvae7 <- pred_loop(2045:2049, nrs_larvae, 130, 
                               temps_gfdl_ssp585,
                               salts_gfdl_ssp585, 7,
                               larval_formula, 1994)
grids_nrslarvae8 <- pred_loop(2050:2054, nrs_larvae, 130, 
                               temps_gfdl_ssp585,
                               salts_gfdl_ssp585, 8,
                               larval_formula, 1994)
grids_nrslarvae9 <- pred_loop(2055:2059, nrs_larvae, 130, 
                               temps_gfdl_ssp585,
                               salts_gfdl_ssp585, 9,
                               larval_formula, 1994)
grids_nrslarvae10 <- pred_loop(2060:2064, nrs_larvae, 130, 
                                temps_gfdl_ssp585,
                                salts_gfdl_ssp585, 10,
                                larval_formula, 1994)
grids_nrslarvae11 <- pred_loop(2065:2069, nrs_larvae, 130, 
                                temps_gfdl_ssp585,
                                salts_gfdl_ssp585, 11,
                                larval_formula, 1994)

# Combine into one data frame
df_nrslarvae5 <- list(grids_nrslarvae6[[1]], grids_nrslarvae6[[2]], grids_nrslarvae6[[3]], 
                       grids_nrslarvae6[[4]], grids_nrslarvae6[[5]], grids_nrslarvae7[[1]], 
                       grids_nrslarvae7[[2]], grids_nrslarvae7[[3]], grids_nrslarvae7[[4]],
                       grids_nrslarvae7[[5]], grids_nrslarvae8[[1]], grids_nrslarvae8[[2]], 
                       grids_nrslarvae8[[3]], grids_nrslarvae8[[4]], grids_nrslarvae8[[5]],
                       grids_nrslarvae9[[1]], grids_nrslarvae9[[2]], grids_nrslarvae9[[3]], 
                       grids_nrslarvae9[[4]], grids_nrslarvae9[[5]], grids_nrslarvae10[[1]], 
                       grids_nrslarvae10[[2]], grids_nrslarvae10[[3]], grids_nrslarvae10[[4]],
                       grids_nrslarvae11[[5]], grids_nrslarvae11[[1]], grids_nrslarvae11[[2]],
                       grids_nrslarvae11[[3]], grids_nrslarvae11[[4]], grids_nrslarvae11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_nrslarvae5), fixed = T)
df_nrslarvae_avg5_gfdl585 <- data.frame(lat = df_nrslarvae5$lat, 
                                         lon = df_nrslarvae5$lon, 
                                         dist = df_nrslarvae5$dist,
                                         avg_pred = rowSums(df_nrslarvae5[, x])/30)
saveRDS(df_nrslarvae_avg5_gfdl585, file = here("data", "df_nrslarvae_avg5_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae_avg5_gfdl585, "Forecasted Distribution 2040 - 2069 \n GFDL SSP585")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_gfdl_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_nrslarvae12 <- pred_loop(2070:2074, nrs_larvae, 130,
                                temps_gfdl_ssp585,
                                salts_gfdl_ssp585, 12,
                                larval_formula, 1994)
grids_nrslarvae13 <- pred_loop(2075:2079, nrs_larvae, 130,
                                temps_gfdl_ssp585,
                                salts_gfdl_ssp585, 13,
                                larval_formula, 1994)
grids_nrslarvae14 <- pred_loop(2080:2084, nrs_larvae, 130,
                                temps_gfdl_ssp585,
                                salts_gfdl_ssp585, 14,
                                larval_formula, 1994)
grids_nrslarvae15 <- pred_loop(2085:2089, nrs_larvae, 130,
                                temps_gfdl_ssp585,
                                salts_gfdl_ssp585, 15,
                                larval_formula, 1994)
grids_nrslarvae16 <- pred_loop(2090:2094, nrs_larvae, 130,
                                temps_gfdl_ssp585,
                                salts_gfdl_ssp585, 16,
                                larval_formula, 1994)
grids_nrslarvae17 <- pred_loop(2095:2099, nrs_larvae, 130, 
                                temps_gfdl_ssp585,
                                salts_gfdl_ssp585, 17,
                                larval_formula, 1994)

# Combine into one data frame
df_nrslarvae6 <- list(grids_nrslarvae12[[1]], grids_nrslarvae12[[2]], grids_nrslarvae12[[3]], 
                       grids_nrslarvae12[[4]], grids_nrslarvae12[[5]], grids_nrslarvae13[[1]], 
                       grids_nrslarvae13[[2]], grids_nrslarvae13[[3]], grids_nrslarvae13[[4]],
                       grids_nrslarvae13[[5]], grids_nrslarvae14[[1]], grids_nrslarvae14[[2]], 
                       grids_nrslarvae14[[3]], grids_nrslarvae14[[4]], grids_nrslarvae14[[5]],
                       grids_nrslarvae15[[1]], grids_nrslarvae15[[2]], grids_nrslarvae15[[3]], 
                       grids_nrslarvae15[[4]], grids_nrslarvae15[[5]], grids_nrslarvae16[[1]], 
                       grids_nrslarvae16[[2]], grids_nrslarvae16[[3]], grids_nrslarvae16[[4]],
                       grids_nrslarvae17[[5]], grids_nrslarvae17[[1]], grids_nrslarvae17[[2]],
                       grids_nrslarvae17[[3]], grids_nrslarvae17[[4]], grids_nrslarvae17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_nrslarvae6), fixed = T)
df_nrslarvae_avg6_gfdl585 <- data.frame(lat = df_nrslarvae6$lat, 
                                         lon = df_nrslarvae6$lon, 
                                         dist = df_nrslarvae6$dist,
                                         avg_pred = rowSums(df_nrslarvae6[, x])/30)
saveRDS(df_nrslarvae_avg6_gfdl585, file = here("data", "df_nrslarvae_avg6_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae_avg6_gfdl585, "Forecasted Distribution 2070 - 2099 \n GFDL SSP585")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_gfdl_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

rm(temps_gfdl_ssp585)
rm(salts_gfdl_ssp585)


##### MIROC 126 ----------------------------------------------------------------------------------------------------------------------------
temps_miroc_ssp126 <- readRDS(here('data', 'temps_miroc_ssp126.rds'))
salts_miroc_ssp126 <- readRDS(here('data', 'salts_miroc_ssp126.rds'))

## 2015 - 2039
grids_nrslarvae1 <- pred_loop(2015:2019, nrs_larvae, 130, 
                               temps_miroc_ssp126,
                               salts_miroc_ssp126, 1,
                               larval_formula, 1994)
grids_nrslarvae2 <- pred_loop(2020:2024, nrs_larvae, 130, 
                               temps_miroc_ssp126,
                               salts_miroc_ssp126, 2,
                               larval_formula, 1994)
grids_nrslarvae3 <- pred_loop(2025:2029, nrs_larvae, 130,
                               temps_miroc_ssp126,
                               salts_miroc_ssp126, 3,
                               larval_formula, 1994)
grids_nrslarvae4 <- pred_loop(2030:2034, nrs_larvae, 130, 
                               temps_miroc_ssp126,
                               salts_miroc_ssp126, 4,
                               larval_formula, 1994)
grids_nrslarvae5 <- pred_loop(2035:2039, nrs_larvae, 130, 
                               temps_miroc_ssp126,
                               salts_miroc_ssp126, 5,
                               larval_formula, 1994)

# Combine into one data frame
df_nrslarvae1 <- list(grids_nrslarvae1[[1]], grids_nrslarvae1[[2]], grids_nrslarvae1[[3]], 
                       grids_nrslarvae1[[4]], grids_nrslarvae1[[5]], grids_nrslarvae2[[1]], 
                       grids_nrslarvae2[[2]], grids_nrslarvae2[[3]], grids_nrslarvae2[[4]],
                       grids_nrslarvae2[[5]], grids_nrslarvae3[[1]], grids_nrslarvae3[[2]], 
                       grids_nrslarvae3[[3]], grids_nrslarvae3[[4]], grids_nrslarvae3[[5]],
                       grids_nrslarvae4[[1]], grids_nrslarvae4[[2]], grids_nrslarvae4[[3]], 
                       grids_nrslarvae4[[4]], grids_nrslarvae4[[5]], grids_nrslarvae5[[1]], 
                       grids_nrslarvae5[[2]], grids_nrslarvae5[[3]], grids_nrslarvae5[[4]],
                       grids_nrslarvae5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_nrslarvae1), fixed = T)
df_nrslarvae_avg1_miroc126 <- data.frame(lat = df_nrslarvae1$lat, 
                                          lon = df_nrslarvae1$lon, 
                                          dist = df_nrslarvae1$dist,
                                          avg_pred = rowSums(df_nrslarvae1[, x])/25)
saveRDS(df_nrslarvae_avg1_miroc126, file = here("data", "df_nrslarvae_avg1_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae_avg1_miroc126, "Forecasted Distribution 2015 - 2039 \n MIROC SSP126")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_miroc_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_nrslarvae6 <- pred_loop(2040:2044, nrs_larvae, 130, 
                               temps_miroc_ssp126,
                               salts_miroc_ssp126, 6,
                               larval_formula, 1994)
grids_nrslarvae7 <- pred_loop(2045:2049, nrs_larvae, 130, 
                               temps_miroc_ssp126,
                               salts_miroc_ssp126, 7,
                               larval_formula, 1994)
grids_nrslarvae8 <- pred_loop(2050:2054, nrs_larvae, 130, 
                               temps_miroc_ssp126,
                               salts_miroc_ssp126, 8,
                               larval_formula, 1994)
grids_nrslarvae9 <- pred_loop(2055:2059, nrs_larvae, 130, 
                               temps_miroc_ssp126,
                               salts_miroc_ssp126, 9,
                               larval_formula, 1994)
grids_nrslarvae10 <- pred_loop(2060:2064, nrs_larvae, 130, 
                                temps_miroc_ssp126,
                                salts_miroc_ssp126, 10,
                                larval_formula, 1994)
grids_nrslarvae11 <- pred_loop(2065:2069, nrs_larvae, 130, 
                                temps_miroc_ssp126,
                                salts_miroc_ssp126, 11,
                                larval_formula, 1994)

# Combine into one data frame
df_nrslarvae2 <- list(grids_nrslarvae6[[1]], grids_nrslarvae6[[2]], grids_nrslarvae6[[3]], 
                       grids_nrslarvae6[[4]], grids_nrslarvae6[[5]], grids_nrslarvae7[[1]], 
                       grids_nrslarvae7[[2]], grids_nrslarvae7[[3]], grids_nrslarvae7[[4]],
                       grids_nrslarvae7[[5]], grids_nrslarvae8[[1]], grids_nrslarvae8[[2]], 
                       grids_nrslarvae8[[3]], grids_nrslarvae8[[4]], grids_nrslarvae8[[5]],
                       grids_nrslarvae9[[1]], grids_nrslarvae9[[2]], grids_nrslarvae9[[3]], 
                       grids_nrslarvae9[[4]], grids_nrslarvae9[[5]], grids_nrslarvae10[[1]], 
                       grids_nrslarvae10[[2]], grids_nrslarvae10[[3]], grids_nrslarvae10[[4]],
                       grids_nrslarvae11[[5]], grids_nrslarvae11[[1]], grids_nrslarvae11[[2]],
                       grids_nrslarvae11[[3]], grids_nrslarvae11[[4]], grids_nrslarvae11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_nrslarvae2), fixed = T)
df_nrslarvae_avg2_miroc126 <- data.frame(lat = df_nrslarvae2$lat, 
                                          lon = df_nrslarvae2$lon, 
                                          dist = df_nrslarvae2$dist,
                                          avg_pred = rowSums(df_nrslarvae2[, x])/30)
saveRDS(df_nrslarvae_avg2_miroc126, file = here("data", "df_nrslarvae_avg2_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae_avg2_miroc126, "Forecasted Distribution 2040 - 2069 \n MIROC SSP126")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_miroc_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_nrslarvae12 <- pred_loop(2070:2074, nrs_larvae, 130,
                                temps_miroc_ssp126,
                                salts_miroc_ssp126, 12,
                                larval_formula, 1994)
grids_nrslarvae13 <- pred_loop(2075:2079, nrs_larvae, 130,
                                temps_miroc_ssp126,
                                salts_miroc_ssp126, 13,
                                larval_formula, 1994)
grids_nrslarvae14 <- pred_loop(2080:2084, nrs_larvae, 130,
                                temps_miroc_ssp126,
                                salts_miroc_ssp126, 14,
                                larval_formula, 1994)
grids_nrslarvae15 <- pred_loop(2085:2089, nrs_larvae, 130,
                                temps_miroc_ssp126,
                                salts_miroc_ssp126, 15,
                                larval_formula, 1994)
grids_nrslarvae16 <- pred_loop(2090:2094, nrs_larvae, 130,
                                temps_miroc_ssp126,
                                salts_miroc_ssp126, 16,
                                larval_formula, 1994)
grids_nrslarvae17 <- pred_loop(2095:2099, nrs_larvae, 130, 
                                temps_miroc_ssp126,
                                salts_miroc_ssp126, 17,
                                larval_formula, 1994)

# Combine into one data frame
df_nrslarvae3 <- list(grids_nrslarvae12[[1]], grids_nrslarvae12[[2]], grids_nrslarvae12[[3]], 
                       grids_nrslarvae12[[4]], grids_nrslarvae12[[5]], grids_nrslarvae13[[1]], 
                       grids_nrslarvae13[[2]], grids_nrslarvae13[[3]], grids_nrslarvae13[[4]],
                       grids_nrslarvae13[[5]], grids_nrslarvae14[[1]], grids_nrslarvae14[[2]], 
                       grids_nrslarvae14[[3]], grids_nrslarvae14[[4]], grids_nrslarvae14[[5]],
                       grids_nrslarvae15[[1]], grids_nrslarvae15[[2]], grids_nrslarvae15[[3]], 
                       grids_nrslarvae15[[4]], grids_nrslarvae15[[5]], grids_nrslarvae16[[1]], 
                       grids_nrslarvae16[[2]], grids_nrslarvae16[[3]], grids_nrslarvae16[[4]],
                       grids_nrslarvae17[[5]], grids_nrslarvae17[[1]], grids_nrslarvae17[[2]],
                       grids_nrslarvae17[[3]], grids_nrslarvae17[[4]], grids_nrslarvae17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_nrslarvae3), fixed = T)
df_nrslarvae_avg3_miroc126 <- data.frame(lat = df_nrslarvae3$lat, 
                                          lon = df_nrslarvae3$lon, 
                                          dist = df_nrslarvae3$dist,
                                          avg_pred = rowSums(df_nrslarvae3[, x])/30)
saveRDS(df_nrslarvae_avg3_miroc126, file = here("data", "df_nrslarvae_avg3_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae_avg3_miroc126, "Forecasted Distribution 2070 - 2099 \n MIROC SSP126")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_miroc_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

rm(temps_miroc_ssp126)
rm(salts_miroc_ssp126)

##### MIROC 585 -----------------------------------------------------------------------------------------------------------------
temps_miroc_ssp585 <- readRDS(here('data', 'temps_miroc_ssp585.rds'))
salts_miroc_ssp585 <- readRDS(here('data', 'salts_miroc_ssp585.rds'))

## 2015 - 2039
grids_nrslarvae1 <- pred_loop(2015:2019, nrs_larvae, 130, 
                               temps_miroc_ssp585,
                               salts_miroc_ssp585, 1,
                               larval_formula, 1994)
grids_nrslarvae2 <- pred_loop(2020:2024, nrs_larvae, 130, 
                               temps_miroc_ssp585,
                               salts_miroc_ssp585, 2,
                               larval_formula, 1994)
grids_nrslarvae3 <- pred_loop(2025:2029, nrs_larvae, 130,
                               temps_miroc_ssp585,
                               salts_miroc_ssp585, 3,
                               larval_formula, 1994)
grids_nrslarvae4 <- pred_loop(2030:2034, nrs_larvae, 130, 
                               temps_miroc_ssp585,
                               salts_miroc_ssp585, 4,
                               larval_formula, 1994)
grids_nrslarvae5 <- pred_loop(2035:2039, nrs_larvae, 130, 
                               temps_miroc_ssp585,
                               salts_miroc_ssp585, 5,
                               larval_formula, 1994)

# Combine into one data frame
df_nrslarvae4 <- list(grids_nrslarvae1[[1]], grids_nrslarvae1[[2]], grids_nrslarvae1[[3]], 
                       grids_nrslarvae1[[4]], grids_nrslarvae1[[5]], grids_nrslarvae2[[1]], 
                       grids_nrslarvae2[[2]], grids_nrslarvae2[[3]], grids_nrslarvae2[[4]],
                       grids_nrslarvae2[[5]], grids_nrslarvae3[[1]], grids_nrslarvae3[[2]], 
                       grids_nrslarvae3[[3]], grids_nrslarvae3[[4]], grids_nrslarvae3[[5]],
                       grids_nrslarvae4[[1]], grids_nrslarvae4[[2]], grids_nrslarvae4[[3]], 
                       grids_nrslarvae4[[4]], grids_nrslarvae4[[5]], grids_nrslarvae5[[1]], 
                       grids_nrslarvae5[[2]], grids_nrslarvae5[[3]], grids_nrslarvae5[[4]],
                       grids_nrslarvae5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_nrslarvae4), fixed = T)
df_nrslarvae_avg4_miroc585 <- data.frame(lat = df_nrslarvae4$lat, 
                                          lon = df_nrslarvae4$lon, 
                                          dist = df_nrslarvae4$dist,
                                          avg_pred = rowSums(df_nrslarvae4[, x])/25)
saveRDS(df_nrslarvae_avg4_miroc585, file = here("data", "df_nrslarvae_avg4_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae_avg4_miroc585, "Forecasted Distribution 2015 - 2039 \n MIROC SSP585")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_miroc_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_nrslarvae6 <- pred_loop(2040:2044, nrs_larvae, 130, 
                               temps_miroc_ssp585,
                               salts_miroc_ssp585, 6,
                               larval_formula, 1994)
grids_nrslarvae7 <- pred_loop(2045:2049, nrs_larvae, 130, 
                               temps_miroc_ssp585,
                               salts_miroc_ssp585, 7,
                               larval_formula, 1994)
grids_nrslarvae8 <- pred_loop(2050:2054, nrs_larvae, 130, 
                               temps_miroc_ssp585,
                               salts_miroc_ssp585, 8,
                               larval_formula, 1994)
grids_nrslarvae9 <- pred_loop(2055:2059, nrs_larvae, 130, 
                               temps_miroc_ssp585,
                               salts_miroc_ssp585, 9,
                               larval_formula, 1994)
grids_nrslarvae10 <- pred_loop(2060:2064, nrs_larvae, 130, 
                                temps_miroc_ssp585,
                                salts_miroc_ssp585, 10,
                                larval_formula, 1994)
grids_nrslarvae11 <- pred_loop(2065:2069, nrs_larvae, 130, 
                                temps_miroc_ssp585,
                                salts_miroc_ssp585, 11,
                                larval_formula, 1994)

# Combine into one data frame
df_nrslarvae5 <- list(grids_nrslarvae6[[1]], grids_nrslarvae6[[2]], grids_nrslarvae6[[3]], 
                       grids_nrslarvae6[[4]], grids_nrslarvae6[[5]], grids_nrslarvae7[[1]], 
                       grids_nrslarvae7[[2]], grids_nrslarvae7[[3]], grids_nrslarvae7[[4]],
                       grids_nrslarvae7[[5]], grids_nrslarvae8[[1]], grids_nrslarvae8[[2]], 
                       grids_nrslarvae8[[3]], grids_nrslarvae8[[4]], grids_nrslarvae8[[5]],
                       grids_nrslarvae9[[1]], grids_nrslarvae9[[2]], grids_nrslarvae9[[3]], 
                       grids_nrslarvae9[[4]], grids_nrslarvae9[[5]], grids_nrslarvae10[[1]], 
                       grids_nrslarvae10[[2]], grids_nrslarvae10[[3]], grids_nrslarvae10[[4]],
                       grids_nrslarvae11[[5]], grids_nrslarvae11[[1]], grids_nrslarvae11[[2]],
                       grids_nrslarvae11[[3]], grids_nrslarvae11[[4]], grids_nrslarvae11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_nrslarvae5), fixed = T)
df_nrslarvae_avg5_miroc585 <- data.frame(lat = df_nrslarvae5$lat, 
                                          lon = df_nrslarvae5$lon, 
                                          dist = df_nrslarvae5$dist,
                                          avg_pred = rowSums(df_nrslarvae5[, x])/30)
saveRDS(df_nrslarvae_avg5_miroc585, file = here("data", "df_nrslarvae_avg5_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae_avg5_miroc585, "Forecasted Distribution 2040 - 2069 \n MIROC SSP585")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_miroc_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_nrslarvae12 <- pred_loop(2070:2074, nrs_larvae, 130,
                                temps_miroc_ssp585,
                                salts_miroc_ssp585, 12,
                                larval_formula, 1994)
grids_nrslarvae13 <- pred_loop(2075:2079, nrs_larvae, 130,
                                temps_miroc_ssp585,
                                salts_miroc_ssp585, 13,
                                larval_formula, 1994)
grids_nrslarvae14 <- pred_loop(2080:2084, nrs_larvae, 130,
                                temps_miroc_ssp585,
                                salts_miroc_ssp585, 14,
                                larval_formula, 1994)
grids_nrslarvae15 <- pred_loop(2085:2089, nrs_larvae, 130,
                                temps_miroc_ssp585,
                                salts_miroc_ssp585, 15,
                                larval_formula, 1994)
grids_nrslarvae16 <- pred_loop(2090:2094, nrs_larvae, 130,
                                temps_miroc_ssp585,
                                salts_miroc_ssp585, 16,
                                larval_formula, 1994)
grids_nrslarvae17 <- pred_loop(2095:2099, nrs_larvae, 130, 
                                temps_miroc_ssp585,
                                salts_miroc_ssp585, 17,
                                larval_formula, 1994)

# Combine into one data frame
df_nrslarvae6 <- list(grids_nrslarvae12[[1]], grids_nrslarvae12[[2]], grids_nrslarvae12[[3]], 
                       grids_nrslarvae12[[4]], grids_nrslarvae12[[5]], grids_nrslarvae13[[1]], 
                       grids_nrslarvae13[[2]], grids_nrslarvae13[[3]], grids_nrslarvae13[[4]],
                       grids_nrslarvae13[[5]], grids_nrslarvae14[[1]], grids_nrslarvae14[[2]], 
                       grids_nrslarvae14[[3]], grids_nrslarvae14[[4]], grids_nrslarvae14[[5]],
                       grids_nrslarvae15[[1]], grids_nrslarvae15[[2]], grids_nrslarvae15[[3]], 
                       grids_nrslarvae15[[4]], grids_nrslarvae15[[5]], grids_nrslarvae16[[1]], 
                       grids_nrslarvae16[[2]], grids_nrslarvae16[[3]], grids_nrslarvae16[[4]],
                       grids_nrslarvae17[[5]], grids_nrslarvae17[[1]], grids_nrslarvae17[[2]],
                       grids_nrslarvae17[[3]], grids_nrslarvae17[[4]], grids_nrslarvae17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_nrslarvae6), fixed = T)
df_nrslarvae_avg6_miroc585 <- data.frame(lat = df_nrslarvae6$lat, 
                                          lon = df_nrslarvae6$lon, 
                                          dist = df_nrslarvae6$dist,
                                          avg_pred = rowSums(df_nrslarvae6[, x])/30)
saveRDS(df_nrslarvae_avg6_miroc585, file = here("data", "df_nrslarvae_avg6_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae_avg6_miroc585, "Forecasted Distribution 2070 - 2099 \n MIROC SSP585")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_miroc_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

rm(temps_miroc_ssp585)
rm(salts_miroc_ssp585)

### Average Predictions ------------------------------------------------------------------------------------------------------------
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
        zlim = c(0, 1000),
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
             zlim = c(0, 1000),
             legend.args = list("Avg. Predicted \n Occurrence",
                                side = 2, cex = 1))
}


#### 2015-2039 ---------------------------------------------------------------------------------------------------------------------
##### Larvae
df_nrslarvae_avg1_cesm126 <- readRDS(here('data', 'df_nrslarvae_avg1_cesm126.rds'))
df_nrslarvae_avg4_cesm585 <- readRDS(here('data', 'df_nrslarvae_avg4_cesm585.rds'))
df_nrslarvae_avg1_gfdl126 <- readRDS(here('data', 'df_nrslarvae_avg1_gfdl126.rds'))
df_nrslarvae_avg4_gfdl585 <- readRDS(here('data', 'df_nrslarvae_avg4_gfdl585.rds'))
df_nrslarvae_avg1_miroc126 <- readRDS(here('data', 'df_nrslarvae_avg1_miroc126.rds'))
df_nrslarvae_avg4_miroc585 <- readRDS(here('data', 'df_nrslarvae_avg4_miroc585.rds'))

df_nrslarvae_merged1 <- list(df_nrslarvae_avg1_cesm126, df_nrslarvae_avg4_cesm585,
                              df_nrslarvae_avg1_gfdl126, df_nrslarvae_avg4_gfdl585,
                              df_nrslarvae_avg1_miroc126, df_nrslarvae_avg4_miroc585) %>%
  reduce(inner_join, by = c("lon", "lat", "dist"))

x <- grepl("pred", names(df_nrslarvae_merged1), fixed = T)
df_nrslarvae_final1 <- data.frame(lat = df_nrslarvae_merged1$lat,
                                   lon = df_nrslarvae_merged1$lon,
                                   avg_pred = (rowSums(df_nrslarvae_merged1[, x])/6))

windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae_final1, "Forecasted Distribution 2015 - 2039")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_avg1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


#### 2040-2069 ---------------------------------------------------------------------------------------------------------------------
##### Larvae
df_nrslarvae_avg2_cesm126 <- readRDS(here('data', 'df_nrslarvae_avg2_cesm126.rds'))
df_nrslarvae_avg5_cesm585 <- readRDS(here('data', 'df_nrslarvae_avg5_cesm585.rds'))
df_nrslarvae_avg2_gfdl126 <- readRDS(here('data', 'df_nrslarvae_avg2_gfdl126.rds'))
df_nrslarvae_avg5_gfdl585 <- readRDS(here('data', 'df_nrslarvae_avg5_gfdl585.rds'))
df_nrslarvae_avg2_miroc126 <- readRDS(here('data', 'df_nrslarvae_avg2_miroc126.rds'))
df_nrslarvae_avg5_miroc585 <- readRDS(here('data', 'df_nrslarvae_avg5_miroc585.rds'))

df_nrslarvae_merged2 <- list(df_nrslarvae_avg2_cesm126, df_nrslarvae_avg5_cesm585,
                              df_nrslarvae_avg2_gfdl126, df_nrslarvae_avg5_gfdl585,
                              df_nrslarvae_avg2_miroc126, df_nrslarvae_avg5_miroc585) %>%
  reduce(inner_join, by = c("lon", "lat", "dist"))

x <- grepl("pred", names(df_nrslarvae_merged2), fixed = T)
df_nrslarvae_final2 <- data.frame(lat = df_nrslarvae_merged2$lat,
                                   lon = df_nrslarvae_merged2$lon,
                                   avg_pred = (rowSums(df_nrslarvae_merged2[, x])/6))

windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae_final2, "Forecasted Distribution 2040 - 2069")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_avg2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


#### 2070-2099----------------------------------------------------------------------------------------------------------------------
df_nrslarvae_avg3_cesm126 <- readRDS(here('data', 'df_nrslarvae_avg3_cesm126.rds'))
df_nrslarvae_avg6_cesm585 <- readRDS(here('data', 'df_nrslarvae_avg6_cesm585.rds'))
df_nrslarvae_avg3_gfdl126 <- readRDS(here('data', 'df_nrslarvae_avg3_gfdl126.rds'))
df_nrslarvae_avg6_gfdl585 <- readRDS(here('data', 'df_nrslarvae_avg6_gfdl585.rds'))
df_nrslarvae_avg3_miroc126 <- readRDS(here('data', 'df_nrslarvae_avg3_miroc126.rds'))
df_nrslarvae_avg6_miroc585 <- readRDS(here('data', 'df_nrslarvae_avg6_miroc585.rds'))

df_nrslarvae_merged3 <- list(df_nrslarvae_avg3_cesm126, df_nrslarvae_avg6_cesm585,
                              df_nrslarvae_avg3_gfdl126, df_nrslarvae_avg6_gfdl585,
                              df_nrslarvae_avg3_miroc126, df_nrslarvae_avg6_miroc585) %>%
  reduce(inner_join, by = c("lon", "lat", "dist"))

x <- grepl("pred", names(df_nrslarvae_merged3), fixed = T)
df_nrslarvae_final3 <- data.frame(lat = df_nrslarvae_merged3$lat,
                                   lon = df_nrslarvae_merged3$lon,
                                   avg_pred = (rowSums(df_nrslarvae_merged3[, x])/6))

windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae_final3, "Forecasted Distribution 2070 - 2099")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_avg3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()