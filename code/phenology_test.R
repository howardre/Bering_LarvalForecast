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
source(here('code/functions', 'doyance_function.R'))


# Load ROMS temperature means and forecast
roms_temps <- readRDS(here('data', 'roms_temps.rds'))

# Load fish data
akp_egg <- as.data.frame(filter((readRDS(here('data', 'akp_egg.rds'))),
                                lat >= 52 & lat <= 62,
                                lon >= -176.5 & lon <= -156.5))
akp_egg$mean_temp <- roms_temps$mean[match(akp_egg$year, roms_temps$year)]
akp_egg$catch <- akp_egg$larvalcatchper10m2 + 1

akp_larvae <- as.data.frame(filter(readRDS(here('data', 'akp_larvae.rds')),
                                   lat >= 52 & lat <= 62,
                                   lon >= -176.5 & lon <= -156.5))
akp_larvae$mean_temp <- roms_temps$mean[match(akp_larvae$year, roms_temps$year)]
akp_larvae$catch <- akp_larvae$larvalcatchper10m2 + 1

# Match ROMS output function
varid_match <- function(data, model_output1, model_output2, list){
  data$roms_date <- NA
  data$roms_temperature <- NA
  data$roms_salinity <- NA
  for (i in 1:nrow(data)) {
    idx_time <- order(abs(model_output1[[list]][[3]] - data$date[i]))[1]
    data$roms_date[i] <- model_output1[[list]][[3]][idx_time]
    idx_grid <- order(doyance_function(data$lat[i],
                                        data$lon[i],
                                        c(model_output1[[list]][[2]]),
                                        c(model_output1[[list]][[1]])))[1]
    data$roms_temperature[i] <- c(model_output1[[list]][[4]][, , idx_time])[idx_grid]
    data$roms_salinity[i] <- c(model_output2[[list]][[4]][, , idx_time])[idx_grid]
  }
  return(data)
}


# Function to get predictions
get_phenology_preds <- function(data, year, date, day,
                              start_date, end_date,
                              temp_output, salt_output,
                              list, formula){
  # Prediction grid
  nlat = 40
  nlon = 60
  latd = seq(min(data$lat), max(data$lat), length.out = nlat)
  lond = seq(min(data$lon), max(data$lon), length.out = nlon)
  
  # Calculate mean temperature
  time_index <- temp_output[[list]][[3]] >= start_date & temp_output[[list]][[3]] <= end_date
  temp_array <- temp_output[[list]][[4]][, , time_index]
  temp_data <- as.data.frame(cbind(lon = as.vector(temp_output[[list]][[1]]), 
                                   lat = as.vector(temp_output[[list]][[2]]), 
                                   temp = as.vector(temp_array)))
  temp_filtered <- temp_data %>% filter(lon >= -170 & lon <= -165, lat >= 56 & lat <= 58)
  mean <- mean(temp_filtered$temp, na.rm = T)

  # Parameterized model
  gam <- formula
  
  # Create phenology grid
  phenology_grid <-  data.frame('lon' = rep(-170, 100),
                                'lat' = rep(57, 100),
                                'doy' = seq(min(data$doy, na.rm = T),
                                            max(data$doy, na.rm = T),
                                            length = 100),
                                'year' = rep(year, 100),
                                'date' = as.Date(day, origin = paste(year, "-01-01", sep = "")))
  phenology_grid <- varid_match(phenology_grid, temp_output, salt_output, list)
  phenology_grid$mean_temp <- mean
  phenology_grid$pred <- predict(gam, newdata = phenology_grid, type = "link")
  phenology_grid$se <- predict(gam, newdata = phenology_grid, se = T, type = "link")[[2]]
  phenology_grid <- phenology_grid[c(1:3, 10, 11)]
  return(phenology_grid)
  
}


egg_formula <- gam(catch ~ s(year) +
                     s(doy, k = 8) +
                     s(lon, lat) +
                     s(roms_temperature, k = 6) +
                     s(roms_salinity, k = 6) +
                     s(doy, by = mean_temp, k = 6),
                   data = akp_egg,
                   family = tw(link = 'log'),
                   method = 'REML')

larval_formula <- gam(catch ~ s(year) +
                        s(doy, k = 8) +
                        s(lon, lat) +
                        s(roms_temperature, k = 6) +
                        s(roms_salinity, k = 6) +
                        s(doy, by = mean_temp, k = 6),
                      data = akp_larvae,
                      family = tw(link = 'log'),
                      method = 'REML')

# Function to loop through years
pred_loop_phenology <- function(range, data, day,
                                temp_output, salt_output,
                                list, formula){
  grids <- list()
  for(j in range) {
    date1 <- paste(j, "-05-17", sep = "")
    date2 <- paste(j, "-02-01", sep = "")
    date3 <- paste(j, "-04-30", sep = "")
    grid <- get_phenology_preds(data, j, date1,
                                day, date2, date3,
                                temp_output, salt_output,
                                list, formula)
    grids[[paste("year", j, sep = "")]] <- grid
  }
  return(grids)
}

# Function to make maps
phenology_curve <- function(grid, title){
  plot(grid$doy,
       grid$avg_pred,
       main = title,
       type = 'l',
       ylab = 'Egg density ln(n/10m2)',
       xlab = 'Day of the year',
       cex.lab = 1.1,
       cex.axis = 1.1,
       cex.main = 1.2,
       xlim = c(min(grid$doy, na.rm = T), max(grid$doy, na.rm = T)),
       ylim = range(c(grid$pred_up, grid$pred_lw)),
       col = 'blue',
       lwd = 2)
  polygon(c(grid$doy, rev(grid$doy)),
          c(grid$pred_lw, rev(grid$pred_up)),
          col = alpha('gray', 0.6),
          lty = 0)
}



### Bias correction testing ----
# Average by month for hindcast and historical found in each fish dataset
# Need to match the forecast to the fish dataset by lat, lon, year, month
# Then subtract the baseline (historical) from the forecast to get the deltas
# Finally add the deltas to the hindcast to get final values



### Plaice Eggs --------------------------------------------------------------------------------------------------------------------------
#### Forecast and average into 3 time periods ---------------------------------------------------------------------------------------------
##### CESM 126 ----------------------------------------------------------------------------------------------------------------------------
temps_cesm_ssp126 <- readRDS(here('data', 'temps_cesm_ssp126.rds'))
salts_cesm_ssp126 <- readRDS(here('data', 'salts_cesm_ssp126.rds'))

## 2015 - 2039
grids_akpegg1 <- pred_loop_phenology(2015:2019, akp_egg, 137, 
                                     temps_cesm_ssp126,
                                     salts_cesm_ssp126, 1,
                                     egg_formula)
grids_akpegg2 <- pred_loop_phenology(2020:2024, akp_egg, 137,
                                     temps_cesm_ssp126,
                                     salts_cesm_ssp126, 2,
                                     egg_formula)
grids_akpegg3 <- pred_loop_phenology(2025:2029, akp_egg, 137,
                                     temps_cesm_ssp126,
                                     salts_cesm_ssp126, 3,
                                     egg_formula)
grids_akpegg4 <- pred_loop_phenology(2030:2034, akp_egg, 137,
                                     temps_cesm_ssp126,
                                     salts_cesm_ssp126, 4,
                                     egg_formula)
grids_akpegg5 <- pred_loop_phenology(2035:2039, akp_egg, 137,
                                     temps_cesm_ssp126,
                                     salts_cesm_ssp126, 5,
                                     egg_formula)

# Combine into one data frame
df_akpegg1 <- list(grids_akpegg1[[1]], grids_akpegg1[[2]], grids_akpegg1[[3]], 
                   grids_akpegg1[[4]], grids_akpegg1[[5]], grids_akpegg2[[1]], 
                   grids_akpegg2[[2]], grids_akpegg2[[3]], grids_akpegg2[[4]],
                   grids_akpegg2[[5]], grids_akpegg3[[1]], grids_akpegg3[[2]], 
                   grids_akpegg3[[3]], grids_akpegg3[[4]], grids_akpegg3[[5]],
                   grids_akpegg4[[1]], grids_akpegg4[[2]], grids_akpegg4[[3]], 
                   grids_akpegg4[[4]], grids_akpegg4[[5]], grids_akpegg5[[1]], 
                   grids_akpegg5[[2]], grids_akpegg5[[3]], grids_akpegg5[[4]],
                   grids_akpegg5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_akpegg1), fixed = T)
y <- grepl("se", names(df_akpegg1), fixed = T)
df_akpegg_avg1_cesm126 <- data.frame(lat = df_akpegg1$lat, 
                                     lon = df_akpegg1$lon, 
                                     doy = df_akpegg1$doy,
                                     avg_pred = rowSums(df_akpegg1[, x])/25,
                                     avg_se = rowSums(df_akpegg1[, y])/25)
df_akpegg_avg1_cesm126$pred_up <- df_akpegg_avg1_cesm126$avg_pred + 1.96 * df_akpegg_avg1_cesm126$avg_se
df_akpegg_avg1_cesm126$pred_lw <- df_akpegg_avg1_cesm126$avg_pred - 1.96 * df_akpegg_avg1_cesm126$avg_se

saveRDS(df_akpegg_avg1_cesm126, file = here("data", "df_akpegg_avg1_cesm126_pheno.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
phenology_curve(df_akpegg_avg1_cesm126, "Forecasted Phenology 2015 - 2039 \n CESM SSP126")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_cesm_ssp126_curve_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
grids_akpegg6 <- pred_loop_phenology(2040:2044, akp_egg, 137, 
                           temps_cesm_ssp126,
                           salts_cesm_ssp126, 6,
                           egg_formula)
grids_akpegg7 <- pred_loop_phenology(2045:2049, akp_egg, 137, 
                           temps_cesm_ssp126,
                           salts_cesm_ssp126, 7,
                           egg_formula)
grids_akpegg8 <- pred_loop_phenology(2050:2054, akp_egg, 137, 
                           temps_cesm_ssp126,
                           salts_cesm_ssp126, 8,
                           egg_formula)
grids_akpegg9 <- pred_loop_phenology(2055:2059, akp_egg, 137, 
                           temps_cesm_ssp126,
                           salts_cesm_ssp126, 9,
                           egg_formula)
grids_akpegg10 <- pred_loop_phenology(2060:2064, akp_egg, 137, 
                            temps_cesm_ssp126,
                            salts_cesm_ssp126, 10,
                            egg_formula)
grids_akpegg11 <- pred_loop_phenology(2065:2069, akp_egg, 137, 
                            temps_cesm_ssp126,
                            salts_cesm_ssp126, 11,
                            egg_formula)

# Combine into one data frame
df_akpegg2 <- list(grids_akpegg6[[1]], grids_akpegg6[[2]], grids_akpegg6[[3]], 
                   grids_akpegg6[[4]], grids_akpegg6[[5]], grids_akpegg7[[1]], 
                   grids_akpegg7[[2]], grids_akpegg7[[3]], grids_akpegg7[[4]],
                   grids_akpegg7[[5]], grids_akpegg8[[1]], grids_akpegg8[[2]], 
                   grids_akpegg8[[3]], grids_akpegg8[[4]], grids_akpegg8[[5]],
                   grids_akpegg9[[1]], grids_akpegg9[[2]], grids_akpegg9[[3]], 
                   grids_akpegg9[[4]], grids_akpegg9[[5]], grids_akpegg10[[1]], 
                   grids_akpegg10[[2]], grids_akpegg10[[3]], grids_akpegg10[[4]],
                   grids_akpegg11[[5]], grids_akpegg11[[1]], grids_akpegg11[[2]],
                   grids_akpegg11[[3]], grids_akpegg11[[4]], grids_akpegg11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_akpegg2), fixed = T)
y <- grepl("se", names(df_akpegg2), fixed = T)
df_akpegg_avg2_cesm126 <- data.frame(lat = df_akpegg2$lat, 
                                     lon = df_akpegg2$lon, 
                                     doy = df_akpegg2$doy,
                                     avg_pred = rowSums(df_akpegg2[, x])/25,
                                     avg_se = rowSums(df_akpegg2[, y])/25)
df_akpegg_avg2_cesm126$pred_up <- df_akpegg_avg2_cesm126$avg_pred + 1.96 * df_akpegg_avg2_cesm126$avg_se
df_akpegg_avg2_cesm126$pred_lw <- df_akpegg_avg2_cesm126$avg_pred - 1.96 * df_akpegg_avg2_cesm126$avg_se

saveRDS(df_akpegg_avg2_cesm126, file = here("data", "df_akpegg_avg2_cesm126_pheno.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
phenology_curve(df_akpegg_avg2_cesm126, "Forecasted Phenology 2025 - 2039 \n CESM SSP126")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_cesm_ssp126_curve_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_akpegg12 <- pred_loop_phenology(2070:2074, akp_egg, 137,
                            temps_cesm_ssp126,
                            salts_cesm_ssp126, 12,
                            egg_formula)
grids_akpegg13 <- pred_loop_phenology(2075:2079, akp_egg, 137,
                            temps_cesm_ssp126,
                            salts_cesm_ssp126, 13,
                            egg_formula)
grids_akpegg14 <- pred_loop_phenology(2080:2084, akp_egg, 137,
                            temps_cesm_ssp126,
                            salts_cesm_ssp126, 14,
                            egg_formula)
grids_akpegg15 <- pred_loop_phenology(2085:2089, akp_egg, 137,
                            temps_cesm_ssp126,
                            salts_cesm_ssp126, 15,
                            egg_formula)
grids_akpegg16 <- pred_loop_phenology(2090:2094, akp_egg, 137,
                            temps_cesm_ssp126,
                            salts_cesm_ssp126, 16,
                            egg_formula)
grids_akpegg17 <- pred_loop_phenology(2095:2099, akp_egg, 137, 
                            temps_cesm_ssp126,
                            salts_cesm_ssp126, 17,
                            egg_formula)

# Combine into one data frame
df_akpegg3 <- list(grids_akpegg12[[1]], grids_akpegg12[[2]], grids_akpegg12[[3]], 
                   grids_akpegg12[[4]], grids_akpegg12[[5]], grids_akpegg13[[1]], 
                   grids_akpegg13[[2]], grids_akpegg13[[3]], grids_akpegg13[[4]],
                   grids_akpegg13[[5]], grids_akpegg14[[1]], grids_akpegg14[[2]], 
                   grids_akpegg14[[3]], grids_akpegg14[[4]], grids_akpegg14[[5]],
                   grids_akpegg15[[1]], grids_akpegg15[[2]], grids_akpegg15[[3]], 
                   grids_akpegg15[[4]], grids_akpegg15[[5]], grids_akpegg16[[1]], 
                   grids_akpegg16[[2]], grids_akpegg16[[3]], grids_akpegg16[[4]],
                   grids_akpegg17[[5]], grids_akpegg17[[1]], grids_akpegg17[[2]],
                   grids_akpegg17[[3]], grids_akpegg17[[4]], grids_akpegg17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_akpegg3), fixed = T)
y <- grepl("se", names(df_akpegg3), fixed = T)
df_akpegg_avg3_cesm126 <- data.frame(lat = df_akpegg3$lat, 
                                     lon = df_akpegg3$lon, 
                                     doy = df_akpegg3$doy,
                                     avg_pred = rowSums(df_akpegg3[, x])/25,
                                     avg_se = rowSums(df_akpegg3[, y])/25)
df_akpegg_avg3_cesm126$pred_up <- df_akpegg_avg3_cesm126$avg_pred + 1.96 * df_akpegg_avg3_cesm126$avg_se
df_akpegg_avg3_cesm126$pred_lw <- df_akpegg_avg3_cesm126$avg_pred - 1.96 * df_akpegg_avg3_cesm126$avg_se

saveRDS(df_akpegg_avg3_cesm126, file = here("data", "df_akpegg_avg3_cesm126_pheno.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
phenology_curve(df_akpegg_avg3_cesm126, "Forecasted Phenology 2035 - 2039 \n CESM SSP126")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_cesm_ssp126_curve_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

##### CESM 585 -----------------------------------------------------------------------------------------------------------------
temps_cesm_ssp585 <- readRDS(here('data', 'temps_cesm_ssp585.rds'))
salts_cesm_ssp585 <- readRDS(here('data', 'salts_cesm_ssp585.rds'))

## 2015 - 2039
grids_akpegg1 <- pred_loop_phenology(2015:2019, akp_egg, 137, 
                           temps_cesm_ssp585,
                           salts_cesm_ssp585, 1,
                           egg_formula)
grids_akpegg2 <- pred_loop_phenology(2020:2024, akp_egg, 137, 
                           temps_cesm_ssp585,
                           salts_cesm_ssp585, 2,
                           egg_formula)
grids_akpegg3 <- pred_loop_phenology(2025:2029, akp_egg, 137,
                           temps_cesm_ssp585,
                           salts_cesm_ssp585, 3,
                           egg_formula)
grids_akpegg4 <- pred_loop_phenology(2030:2034, akp_egg, 137, 
                           temps_cesm_ssp585,
                           salts_cesm_ssp585, 4,
                           egg_formula)
grids_akpegg5 <- pred_loop_phenology(2035:2039, akp_egg, 137, 
                           temps_cesm_ssp585,
                           salts_cesm_ssp585, 5,
                           egg_formula)

# Combine into one data frame
df_akpegg4 <- list(grids_akpegg1[[1]], grids_akpegg1[[2]], grids_akpegg1[[3]], 
                   grids_akpegg1[[4]], grids_akpegg1[[5]], grids_akpegg2[[1]], 
                   grids_akpegg2[[2]], grids_akpegg2[[3]], grids_akpegg2[[4]],
                   grids_akpegg2[[5]], grids_akpegg3[[1]], grids_akpegg3[[2]], 
                   grids_akpegg3[[3]], grids_akpegg3[[4]], grids_akpegg3[[5]],
                   grids_akpegg4[[1]], grids_akpegg4[[2]], grids_akpegg4[[3]], 
                   grids_akpegg4[[4]], grids_akpegg4[[5]], grids_akpegg5[[1]], 
                   grids_akpegg5[[2]], grids_akpegg5[[3]], grids_akpegg5[[4]],
                   grids_akpegg5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_akpegg1), fixed = T)
y <- grepl("se", names(df_akpegg1), fixed = T)
df_akpegg_avg1_cesm585 <- data.frame(lat = df_akpegg1$lat, 
                                     lon = df_akpegg1$lon, 
                                     doy = df_akpegg1$doy,
                                     avg_pred = rowSums(df_akpegg1[, x])/25,
                                     avg_se = rowSums(df_akpegg1[, y])/25)
df_akpegg_avg1_cesm585$pred_up <- df_akpegg_avg1_cesm585$avg_pred + 1.96 * df_akpegg_avg1_cesm585$avg_se
df_akpegg_avg1_cesm585$pred_lw <- df_akpegg_avg1_cesm585$avg_pred - 1.96 * df_akpegg_avg1_cesm585$avg_se

saveRDS(df_akpegg_avg1_cesm585, file = here("data", "df_akpegg_avg1_cesm585_pheno.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
phenology_curve(df_akpegg_avg1_cesm585, "Forecasted Phenology 2015 - 2019 \n CESM SSP585")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_cesm_ssp585_curve_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_akpegg6 <- pred_loop_phenology(2040:2044, akp_egg, 137, 
                           temps_cesm_ssp585,
                           salts_cesm_ssp585, 6,
                           egg_formula)
grids_akpegg7 <- pred_loop_phenology(2045:2049, akp_egg, 137, 
                           temps_cesm_ssp585,
                           salts_cesm_ssp585, 7,
                           egg_formula)
grids_akpegg8 <- pred_loop_phenology(2050:2054, akp_egg, 137, 
                           temps_cesm_ssp585,
                           salts_cesm_ssp585, 8,
                           egg_formula)
grids_akpegg9 <- pred_loop_phenology(2055:2059, akp_egg, 137, 
                           temps_cesm_ssp585,
                           salts_cesm_ssp585, 9,
                           egg_formula)
grids_akpegg10 <- pred_loop_phenology(2060:2064, akp_egg, 137, 
                            temps_cesm_ssp585,
                            salts_cesm_ssp585, 10,
                            egg_formula)
grids_akpegg11 <- pred_loop_phenology(2065:2069, akp_egg, 137, 
                            temps_cesm_ssp585,
                            salts_cesm_ssp585, 11,
                            egg_formula)

# Combine into one data frame
df_akpegg5 <- list(grids_akpegg6[[1]], grids_akpegg6[[2]], grids_akpegg6[[3]], 
                   grids_akpegg6[[4]], grids_akpegg6[[5]], grids_akpegg7[[1]], 
                   grids_akpegg7[[2]], grids_akpegg7[[3]], grids_akpegg7[[4]],
                   grids_akpegg7[[5]], grids_akpegg8[[1]], grids_akpegg8[[2]], 
                   grids_akpegg8[[3]], grids_akpegg8[[4]], grids_akpegg8[[5]],
                   grids_akpegg9[[1]], grids_akpegg9[[2]], grids_akpegg9[[3]], 
                   grids_akpegg9[[4]], grids_akpegg9[[5]], grids_akpegg10[[1]], 
                   grids_akpegg10[[2]], grids_akpegg10[[3]], grids_akpegg10[[4]],
                   grids_akpegg11[[5]], grids_akpegg11[[1]], grids_akpegg11[[2]],
                   grids_akpegg11[[3]], grids_akpegg11[[4]], grids_akpegg11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_akpegg5), fixed = T)
df_akpegg_avg5_cesm585 <- data.frame(lat = df_akpegg5$lat, 
                                     lon = df_akpegg5$lon, 
                                     doy = df_akpegg5$doy,
                                     avg_pred = rowSums(df_akpegg5[, x])/30)
saveRDS(df_akpegg_avg5_cesm585, file = here("data", "df_akpegg_avg5_cesm585_pheno.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akpegg_avg5_cesm585, "Forecasted doyribution 2040 - 2069 \n CESM SSP585")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_cesm_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_akpegg12 <- pred_loop_phenology(2070:2074, akp_egg, 137,
                            temps_cesm_ssp585,
                            salts_cesm_ssp585, 12,
                            egg_formula)
grids_akpegg13 <- pred_loop_phenology(2075:2079, akp_egg, 137,
                            temps_cesm_ssp585,
                            salts_cesm_ssp585, 13,
                            egg_formula)
grids_akpegg14 <- pred_loop_phenology(2080:2084, akp_egg, 137,
                            temps_cesm_ssp585,
                            salts_cesm_ssp585, 14,
                            egg_formula)
grids_akpegg15 <- pred_loop_phenology(2085:2089, akp_egg, 137,
                            temps_cesm_ssp585,
                            salts_cesm_ssp585, 15,
                            egg_formula)
grids_akpegg16 <- pred_loop_phenology(2090:2094, akp_egg, 137,
                            temps_cesm_ssp585,
                            salts_cesm_ssp585, 16,
                            egg_formula)
grids_akpegg17 <- pred_loop_phenology(2095:2099, akp_egg, 137, 
                            temps_cesm_ssp585,
                            salts_cesm_ssp585, 17,
                            egg_formula)

# Combine into one data frame
df_akpegg6 <- list(grids_akpegg12[[1]], grids_akpegg12[[2]], grids_akpegg12[[3]], 
                   grids_akpegg12[[4]], grids_akpegg12[[5]], grids_akpegg13[[1]], 
                   grids_akpegg13[[2]], grids_akpegg13[[3]], grids_akpegg13[[4]],
                   grids_akpegg13[[5]], grids_akpegg14[[1]], grids_akpegg14[[2]], 
                   grids_akpegg14[[3]], grids_akpegg14[[4]], grids_akpegg14[[5]],
                   grids_akpegg15[[1]], grids_akpegg15[[2]], grids_akpegg15[[3]], 
                   grids_akpegg15[[4]], grids_akpegg15[[5]], grids_akpegg16[[1]], 
                   grids_akpegg16[[2]], grids_akpegg16[[3]], grids_akpegg16[[4]],
                   grids_akpegg17[[5]], grids_akpegg17[[1]], grids_akpegg17[[2]],
                   grids_akpegg17[[3]], grids_akpegg17[[4]], grids_akpegg17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_akpegg6), fixed = T)
df_akpegg_avg6_cesm585 <- data.frame(lat = df_akpegg6$lat, 
                                     lon = df_akpegg6$lon, 
                                     doy = df_akpegg6$doy,
                                     avg_pred = rowSums(df_akpegg6[, x])/30)
saveRDS(df_akpegg_avg6_cesm585, file = here("data", "df_akpegg_avg6_cesm585_pheno.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akpegg_avg6_cesm585, "Forecasted doyribution 2070 - 2099 \n CESM SSP585")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_cesm_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### GFDL 126 ----------------------------------------------------------------------------------------------------------------------------
temps_gfdl_ssp126 <- readRDS(here('data', 'temps_gfdl_ssp126.rds'))
salts_gfdl_ssp126 <- readRDS(here('data', 'salts_gfdl_ssp126.rds'))

## 2015 - 2039
grids_akpegg1 <- pred_loop_phenology(2015:2019, akp_egg, 137, 
                           temps_gfdl_ssp126,
                           salts_gfdl_ssp126, 1,
                           egg_formula)
grids_akpegg2 <- pred_loop_phenology(2020:2024, akp_egg, 137, 
                           temps_gfdl_ssp126,
                           salts_gfdl_ssp126, 2,
                           egg_formula)
grids_akpegg3 <- pred_loop_phenology(2025:2029, akp_egg, 137,
                           temps_gfdl_ssp126,
                           salts_gfdl_ssp126, 3,
                           egg_formula)
grids_akpegg4 <- pred_loop_phenology(2030:2034, akp_egg, 137, 
                           temps_gfdl_ssp126,
                           salts_gfdl_ssp126, 4,
                           egg_formula)
grids_akpegg5 <- pred_loop_phenology(2035:2039, akp_egg, 137, 
                           temps_gfdl_ssp126,
                           salts_gfdl_ssp126, 5,
                           egg_formula)

# Combine into one data frame
df_akpegg1 <- list(grids_akpegg1[[1]], grids_akpegg1[[2]], grids_akpegg1[[3]], 
                   grids_akpegg1[[4]], grids_akpegg1[[5]], grids_akpegg2[[1]], 
                   grids_akpegg2[[2]], grids_akpegg2[[3]], grids_akpegg2[[4]],
                   grids_akpegg2[[5]], grids_akpegg3[[1]], grids_akpegg3[[2]], 
                   grids_akpegg3[[3]], grids_akpegg3[[4]], grids_akpegg3[[5]],
                   grids_akpegg4[[1]], grids_akpegg4[[2]], grids_akpegg4[[3]], 
                   grids_akpegg4[[4]], grids_akpegg4[[5]], grids_akpegg5[[1]], 
                   grids_akpegg5[[2]], grids_akpegg5[[3]], grids_akpegg5[[4]],
                   grids_akpegg5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_akpegg1), fixed = T)
df_akpegg_avg1_gfdl126 <- data.frame(lat = df_akpegg1$lat, 
                                     lon = df_akpegg1$lon, 
                                     doy = df_akpegg1$doy,
                                     avg_pred = rowSums(df_akpegg1[, x])/25)
saveRDS(df_akpegg_avg1_gfdl126, file = here("data", "df_akpegg_avg1_gfdl126_pheno.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akpegg_avg1_gfdl126, "Forecasted doyribution 2015 - 2039 \n GFDL SSP126")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_gfdl_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_akpegg6 <- pred_loop_phenology(2040:2044, akp_egg, 137, 
                           temps_gfdl_ssp126,
                           salts_gfdl_ssp126, 6,
                           egg_formula)
grids_akpegg7 <- pred_loop_phenology(2045:2049, akp_egg, 137, 
                           temps_gfdl_ssp126,
                           salts_gfdl_ssp126, 7,
                           egg_formula)
grids_akpegg8 <- pred_loop_phenology(2050:2054, akp_egg, 137, 
                           temps_gfdl_ssp126,
                           salts_gfdl_ssp126, 8,
                           egg_formula)
grids_akpegg9 <- pred_loop_phenology(2055:2059, akp_egg, 137, 
                           temps_gfdl_ssp126,
                           salts_gfdl_ssp126, 9,
                           egg_formula)
grids_akpegg10 <- pred_loop_phenology(2060:2064, akp_egg, 137, 
                            temps_gfdl_ssp126,
                            salts_gfdl_ssp126, 10,
                            egg_formula)
grids_akpegg11 <- pred_loop_phenology(2065:2069, akp_egg, 137, 
                            temps_gfdl_ssp126,
                            salts_gfdl_ssp126, 11,
                            egg_formula)

# Combine into one data frame
df_akpegg2 <- list(grids_akpegg6[[1]], grids_akpegg6[[2]], grids_akpegg6[[3]], 
                   grids_akpegg6[[4]], grids_akpegg6[[5]], grids_akpegg7[[1]], 
                   grids_akpegg7[[2]], grids_akpegg7[[3]], grids_akpegg7[[4]],
                   grids_akpegg7[[5]], grids_akpegg8[[1]], grids_akpegg8[[2]], 
                   grids_akpegg8[[3]], grids_akpegg8[[4]], grids_akpegg8[[5]],
                   grids_akpegg9[[1]], grids_akpegg9[[2]], grids_akpegg9[[3]], 
                   grids_akpegg9[[4]], grids_akpegg9[[5]], grids_akpegg10[[1]], 
                   grids_akpegg10[[2]], grids_akpegg10[[3]], grids_akpegg10[[4]],
                   grids_akpegg11[[5]], grids_akpegg11[[1]], grids_akpegg11[[2]],
                   grids_akpegg11[[3]], grids_akpegg11[[4]], grids_akpegg11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_akpegg2), fixed = T)
df_akpegg_avg2_gfdl126 <- data.frame(lat = df_akpegg2$lat, 
                                     lon = df_akpegg2$lon, 
                                     doy = df_akpegg2$doy,
                                     avg_pred = rowSums(df_akpegg2[, x])/30)
saveRDS(df_akpegg_avg2_gfdl126, file = here("data", "df_akpegg_avg2_gfdl126_pheno.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akpegg_avg2_gfdl126, "Forecasted doyribution 2040 - 2069 \n GFDL SSP126")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_gfdl_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_akpegg12 <- pred_loop_phenology(2070:2074, akp_egg, 137,
                            temps_gfdl_ssp126,
                            salts_gfdl_ssp126, 12,
                            egg_formula)
grids_akpegg13 <- pred_loop_phenology(2075:2079, akp_egg, 137,
                            temps_gfdl_ssp126,
                            salts_gfdl_ssp126, 13,
                            egg_formula)
grids_akpegg14 <- pred_loop_phenology(2080:2084, akp_egg, 137,
                            temps_gfdl_ssp126,
                            salts_gfdl_ssp126, 14,
                            egg_formula)
grids_akpegg15 <- pred_loop_phenology(2085:2089, akp_egg, 137,
                            temps_gfdl_ssp126,
                            salts_gfdl_ssp126, 15,
                            egg_formula)
grids_akpegg16 <- pred_loop_phenology(2090:2094, akp_egg, 137,
                            temps_gfdl_ssp126,
                            salts_gfdl_ssp126, 16,
                            egg_formula)
grids_akpegg17 <- pred_loop_phenology(2095:2099, akp_egg, 137, 
                            temps_gfdl_ssp126,
                            salts_gfdl_ssp126, 17,
                            egg_formula)

# Combine into one data frame
df_akpegg3 <- list(grids_akpegg12[[1]], grids_akpegg12[[2]], grids_akpegg12[[3]], 
                   grids_akpegg12[[4]], grids_akpegg12[[5]], grids_akpegg13[[1]], 
                   grids_akpegg13[[2]], grids_akpegg13[[3]], grids_akpegg13[[4]],
                   grids_akpegg13[[5]], grids_akpegg14[[1]], grids_akpegg14[[2]], 
                   grids_akpegg14[[3]], grids_akpegg14[[4]], grids_akpegg14[[5]],
                   grids_akpegg15[[1]], grids_akpegg15[[2]], grids_akpegg15[[3]], 
                   grids_akpegg15[[4]], grids_akpegg15[[5]], grids_akpegg16[[1]], 
                   grids_akpegg16[[2]], grids_akpegg16[[3]], grids_akpegg16[[4]],
                   grids_akpegg17[[5]], grids_akpegg17[[1]], grids_akpegg17[[2]],
                   grids_akpegg17[[3]], grids_akpegg17[[4]], grids_akpegg17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_akpegg3), fixed = T)
df_akpegg_avg3_gfdl126 <- data.frame(lat = df_akpegg3$lat, 
                                     lon = df_akpegg3$lon, 
                                     doy = df_akpegg3$doy,
                                     avg_pred = rowSums(df_akpegg3[, x])/30)
saveRDS(df_akpegg_avg3_gfdl126, file = here("data", "df_akpegg_avg3_gfdl126_pheno.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akpegg_avg3_gfdl126, "Forecasted doyribution 2070 - 2099 \n GFDL SSP126")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_gfdl_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

##### GFDL 585 -----------------------------------------------------------------------------------------------------------------
temps_gfdl_ssp585 <- readRDS(here('data', 'temps_gfdl_ssp585.rds'))
salts_gfdl_ssp585 <- readRDS(here('data', 'salts_gfdl_ssp585.rds'))

## 2015 - 2039
grids_akpegg1 <- pred_loop_phenology(2015:2019, akp_egg, 137, 
                           temps_gfdl_ssp585,
                           salts_gfdl_ssp585, 1,
                           egg_formula)
grids_akpegg2 <- pred_loop_phenology(2020:2024, akp_egg, 137, 
                           temps_gfdl_ssp585,
                           salts_gfdl_ssp585, 2,
                           egg_formula)
grids_akpegg3 <- pred_loop_phenology(2025:2029, akp_egg, 137,
                           temps_gfdl_ssp585,
                           salts_gfdl_ssp585, 3,
                           egg_formula)
grids_akpegg4 <- pred_loop_phenology(2030:2034, akp_egg, 137, 
                           temps_gfdl_ssp585,
                           salts_gfdl_ssp585, 4,
                           egg_formula)
grids_akpegg5 <- pred_loop_phenology(2035:2039, akp_egg, 137, 
                           temps_gfdl_ssp585,
                           salts_gfdl_ssp585, 5,
                           egg_formula)

# Combine into one data frame
df_akpegg4 <- list(grids_akpegg1[[1]], grids_akpegg1[[2]], grids_akpegg1[[3]], 
                   grids_akpegg1[[4]], grids_akpegg1[[5]], grids_akpegg2[[1]], 
                   grids_akpegg2[[2]], grids_akpegg2[[3]], grids_akpegg2[[4]],
                   grids_akpegg2[[5]], grids_akpegg3[[1]], grids_akpegg3[[2]], 
                   grids_akpegg3[[3]], grids_akpegg3[[4]], grids_akpegg3[[5]],
                   grids_akpegg4[[1]], grids_akpegg4[[2]], grids_akpegg4[[3]], 
                   grids_akpegg4[[4]], grids_akpegg4[[5]], grids_akpegg5[[1]], 
                   grids_akpegg5[[2]], grids_akpegg5[[3]], grids_akpegg5[[4]],
                   grids_akpegg5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_akpegg4), fixed = T)
df_akpegg_avg4_gfdl585 <- data.frame(lat = df_akpegg4$lat, 
                                     lon = df_akpegg4$lon, 
                                     doy = df_akpegg4$doy,
                                     avg_pred = rowSums(df_akpegg4[, x])/25)
saveRDS(df_akpegg_avg4_gfdl585, file = here("data", "df_akpegg_avg4_gfdl585_pheno.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akpegg_avg4_gfdl585, "Forecasted doyribution 2015 - 2039 \n GFDL SSP585")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_gfdl_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_akpegg6 <- pred_loop_phenology(2040:2044, akp_egg, 137, 
                           temps_gfdl_ssp585,
                           salts_gfdl_ssp585, 6,
                           egg_formula)
grids_akpegg7 <- pred_loop_phenology(2045:2049, akp_egg, 137, 
                           temps_gfdl_ssp585,
                           salts_gfdl_ssp585, 7,
                           egg_formula)
grids_akpegg8 <- pred_loop_phenology(2050:2054, akp_egg, 137, 
                           temps_gfdl_ssp585,
                           salts_gfdl_ssp585, 8,
                           egg_formula)
grids_akpegg9 <- pred_loop_phenology(2055:2059, akp_egg, 137, 
                           temps_gfdl_ssp585,
                           salts_gfdl_ssp585, 9,
                           egg_formula)
grids_akpegg10 <- pred_loop_phenology(2060:2064, akp_egg, 137, 
                            temps_gfdl_ssp585,
                            salts_gfdl_ssp585, 10,
                            egg_formula)
grids_akpegg11 <- pred_loop_phenology(2065:2069, akp_egg, 137, 
                            temps_gfdl_ssp585,
                            salts_gfdl_ssp585, 11,
                            egg_formula)

# Combine into one data frame
df_akpegg5 <- list(grids_akpegg6[[1]], grids_akpegg6[[2]], grids_akpegg6[[3]], 
                   grids_akpegg6[[4]], grids_akpegg6[[5]], grids_akpegg7[[1]], 
                   grids_akpegg7[[2]], grids_akpegg7[[3]], grids_akpegg7[[4]],
                   grids_akpegg7[[5]], grids_akpegg8[[1]], grids_akpegg8[[2]], 
                   grids_akpegg8[[3]], grids_akpegg8[[4]], grids_akpegg8[[5]],
                   grids_akpegg9[[1]], grids_akpegg9[[2]], grids_akpegg9[[3]], 
                   grids_akpegg9[[4]], grids_akpegg9[[5]], grids_akpegg10[[1]], 
                   grids_akpegg10[[2]], grids_akpegg10[[3]], grids_akpegg10[[4]],
                   grids_akpegg11[[5]], grids_akpegg11[[1]], grids_akpegg11[[2]],
                   grids_akpegg11[[3]], grids_akpegg11[[4]], grids_akpegg11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_akpegg5), fixed = T)
df_akpegg_avg5_gfdl585 <- data.frame(lat = df_akpegg5$lat, 
                                     lon = df_akpegg5$lon, 
                                     doy = df_akpegg5$doy,
                                     avg_pred = rowSums(df_akpegg5[, x])/30)
saveRDS(df_akpegg_avg5_gfdl585, file = here("data", "df_akpegg_avg5_gfdl585_pheno.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akpegg_avg5_gfdl585, "Forecasted doyribution 2040 - 2069 \n GFDL SSP585")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_gfdl_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_akpegg12 <- pred_loop_phenology(2070:2074, akp_egg, 137,
                            temps_gfdl_ssp585,
                            salts_gfdl_ssp585, 12,
                            egg_formula)
grids_akpegg13 <- pred_loop_phenology(2075:2079, akp_egg, 137,
                            temps_gfdl_ssp585,
                            salts_gfdl_ssp585, 13,
                            egg_formula)
grids_akpegg14 <- pred_loop_phenology(2080:2084, akp_egg, 137,
                            temps_gfdl_ssp585,
                            salts_gfdl_ssp585, 14,
                            egg_formula)
grids_akpegg15 <- pred_loop_phenology(2085:2089, akp_egg, 137,
                            temps_gfdl_ssp585,
                            salts_gfdl_ssp585, 15,
                            egg_formula)
grids_akpegg16 <- pred_loop_phenology(2090:2094, akp_egg, 137,
                            temps_gfdl_ssp585,
                            salts_gfdl_ssp585, 16,
                            egg_formula)
grids_akpegg17 <- pred_loop_phenology(2095:2099, akp_egg, 137, 
                            temps_gfdl_ssp585,
                            salts_gfdl_ssp585, 17,
                            egg_formula)

# Combine into one data frame
df_akpegg6 <- list(grids_akpegg12[[1]], grids_akpegg12[[2]], grids_akpegg12[[3]], 
                   grids_akpegg12[[4]], grids_akpegg12[[5]], grids_akpegg13[[1]], 
                   grids_akpegg13[[2]], grids_akpegg13[[3]], grids_akpegg13[[4]],
                   grids_akpegg13[[5]], grids_akpegg14[[1]], grids_akpegg14[[2]], 
                   grids_akpegg14[[3]], grids_akpegg14[[4]], grids_akpegg14[[5]],
                   grids_akpegg15[[1]], grids_akpegg15[[2]], grids_akpegg15[[3]], 
                   grids_akpegg15[[4]], grids_akpegg15[[5]], grids_akpegg16[[1]], 
                   grids_akpegg16[[2]], grids_akpegg16[[3]], grids_akpegg16[[4]],
                   grids_akpegg17[[5]], grids_akpegg17[[1]], grids_akpegg17[[2]],
                   grids_akpegg17[[3]], grids_akpegg17[[4]], grids_akpegg17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_akpegg6), fixed = T)
df_akpegg_avg6_gfdl585 <- data.frame(lat = df_akpegg6$lat, 
                                     lon = df_akpegg6$lon, 
                                     doy = df_akpegg6$doy,
                                     avg_pred = rowSums(df_akpegg6[, x])/30)
saveRDS(df_akpegg_avg6_gfdl585, file = here("data", "df_akpegg_avg6_gfdl585_pheno.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akpegg_avg6_gfdl585, "Forecasted doyribution 2070 - 2099 \n GFDL SSP585")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_gfdl_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### MIROC 126 ----------------------------------------------------------------------------------------------------------------------------
temps_miroc_ssp126 <- readRDS(here('data', 'temps_miroc_ssp126.rds'))
salts_miroc_ssp126 <- readRDS(here('data', 'salts_miroc_ssp126.rds'))

## 2015 - 2039
grids_akpegg1 <- pred_loop_phenology(2015:2019, akp_egg, 137, 
                           temps_miroc_ssp126,
                           salts_miroc_ssp126, 1,
                           egg_formula)
grids_akpegg2 <- pred_loop_phenology(2020:2024, akp_egg, 137, 
                           temps_miroc_ssp126,
                           salts_miroc_ssp126, 2,
                           egg_formula)
grids_akpegg3 <- pred_loop_phenology(2025:2029, akp_egg, 137,
                           temps_miroc_ssp126,
                           salts_miroc_ssp126, 3,
                           egg_formula)
grids_akpegg4 <- pred_loop_phenology(2030:2034, akp_egg, 137, 
                           temps_miroc_ssp126,
                           salts_miroc_ssp126, 4,
                           egg_formula)
grids_akpegg5 <- pred_loop_phenology(2035:2039, akp_egg, 137, 
                           temps_miroc_ssp126,
                           salts_miroc_ssp126, 5,
                           egg_formula)

# Combine into one data frame
df_akpegg1 <- list(grids_akpegg1[[1]], grids_akpegg1[[2]], grids_akpegg1[[3]], 
                   grids_akpegg1[[4]], grids_akpegg1[[5]], grids_akpegg2[[1]], 
                   grids_akpegg2[[2]], grids_akpegg2[[3]], grids_akpegg2[[4]],
                   grids_akpegg2[[5]], grids_akpegg3[[1]], grids_akpegg3[[2]], 
                   grids_akpegg3[[3]], grids_akpegg3[[4]], grids_akpegg3[[5]],
                   grids_akpegg4[[1]], grids_akpegg4[[2]], grids_akpegg4[[3]], 
                   grids_akpegg4[[4]], grids_akpegg4[[5]], grids_akpegg5[[1]], 
                   grids_akpegg5[[2]], grids_akpegg5[[3]], grids_akpegg5[[4]],
                   grids_akpegg5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_akpegg1), fixed = T)
df_akpegg_avg1_miroc126 <- data.frame(lat = df_akpegg1$lat, 
                                      lon = df_akpegg1$lon, 
                                      doy = df_akpegg1$doy,
                                      avg_pred = rowSums(df_akpegg1[, x])/25)
saveRDS(df_akpegg_avg1_miroc126, file = here("data", "df_akpegg_avg1_miroc126_pheno.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akpegg_avg1_miroc126, "Forecasted doyribution 2015 - 2039 \n MIROC SSP126")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_miroc_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_akpegg6 <- pred_loop_phenology(2040:2044, akp_egg, 137, 
                           temps_miroc_ssp126,
                           salts_miroc_ssp126, 6,
                           egg_formula)
grids_akpegg7 <- pred_loop_phenology(2045:2049, akp_egg, 137, 
                           temps_miroc_ssp126,
                           salts_miroc_ssp126, 7,
                           egg_formula)
grids_akpegg8 <- pred_loop_phenology(2050:2054, akp_egg, 137, 
                           temps_miroc_ssp126,
                           salts_miroc_ssp126, 8,
                           egg_formula)
grids_akpegg9 <- pred_loop_phenology(2055:2059, akp_egg, 137, 
                           temps_miroc_ssp126,
                           salts_miroc_ssp126, 9,
                           egg_formula)
grids_akpegg10 <- pred_loop_phenology(2060:2064, akp_egg, 137, 
                            temps_miroc_ssp126,
                            salts_miroc_ssp126, 10,
                            egg_formula)
grids_akpegg11 <- pred_loop_phenology(2065:2069, akp_egg, 137, 
                            temps_miroc_ssp126,
                            salts_miroc_ssp126, 11,
                            egg_formula)

# Combine into one data frame
df_akpegg2 <- list(grids_akpegg6[[1]], grids_akpegg6[[2]], grids_akpegg6[[3]], 
                   grids_akpegg6[[4]], grids_akpegg6[[5]], grids_akpegg7[[1]], 
                   grids_akpegg7[[2]], grids_akpegg7[[3]], grids_akpegg7[[4]],
                   grids_akpegg7[[5]], grids_akpegg8[[1]], grids_akpegg8[[2]], 
                   grids_akpegg8[[3]], grids_akpegg8[[4]], grids_akpegg8[[5]],
                   grids_akpegg9[[1]], grids_akpegg9[[2]], grids_akpegg9[[3]], 
                   grids_akpegg9[[4]], grids_akpegg9[[5]], grids_akpegg10[[1]], 
                   grids_akpegg10[[2]], grids_akpegg10[[3]], grids_akpegg10[[4]],
                   grids_akpegg11[[5]], grids_akpegg11[[1]], grids_akpegg11[[2]],
                   grids_akpegg11[[3]], grids_akpegg11[[4]], grids_akpegg11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_akpegg2), fixed = T)
df_akpegg_avg2_miroc126 <- data.frame(lat = df_akpegg2$lat, 
                                      lon = df_akpegg2$lon, 
                                      doy = df_akpegg2$doy,
                                      avg_pred = rowSums(df_akpegg2[, x])/30)
saveRDS(df_akpegg_avg2_miroc126, file = here("data", "df_akpegg_avg2_miroc126_pheno.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akpegg_avg2_miroc126, "Forecasted doyribution 2040 - 2069 \n MIROC SSP126")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_miroc_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_akpegg12 <- pred_loop_phenology(2070:2074, akp_egg, 137,
                            temps_miroc_ssp126,
                            salts_miroc_ssp126, 12,
                            egg_formula)
grids_akpegg13 <- pred_loop_phenology(2075:2079, akp_egg, 137,
                            temps_miroc_ssp126,
                            salts_miroc_ssp126, 13,
                            egg_formula)
grids_akpegg14 <- pred_loop_phenology(2080:2084, akp_egg, 137,
                            temps_miroc_ssp126,
                            salts_miroc_ssp126, 14,
                            egg_formula)
grids_akpegg15 <- pred_loop_phenology(2085:2089, akp_egg, 137,
                            temps_miroc_ssp126,
                            salts_miroc_ssp126, 15,
                            egg_formula)
grids_akpegg16 <- pred_loop_phenology(2090:2094, akp_egg, 137,
                            temps_miroc_ssp126,
                            salts_miroc_ssp126, 16,
                            egg_formula)
grids_akpegg17 <- pred_loop_phenology(2095:2099, akp_egg, 137, 
                            temps_miroc_ssp126,
                            salts_miroc_ssp126, 17,
                            egg_formula)

# Combine into one data frame
df_akpegg3 <- list(grids_akpegg12[[1]], grids_akpegg12[[2]], grids_akpegg12[[3]], 
                   grids_akpegg12[[4]], grids_akpegg12[[5]], grids_akpegg13[[1]], 
                   grids_akpegg13[[2]], grids_akpegg13[[3]], grids_akpegg13[[4]],
                   grids_akpegg13[[5]], grids_akpegg14[[1]], grids_akpegg14[[2]], 
                   grids_akpegg14[[3]], grids_akpegg14[[4]], grids_akpegg14[[5]],
                   grids_akpegg15[[1]], grids_akpegg15[[2]], grids_akpegg15[[3]], 
                   grids_akpegg15[[4]], grids_akpegg15[[5]], grids_akpegg16[[1]], 
                   grids_akpegg16[[2]], grids_akpegg16[[3]], grids_akpegg16[[4]],
                   grids_akpegg17[[5]], grids_akpegg17[[1]], grids_akpegg17[[2]],
                   grids_akpegg17[[3]], grids_akpegg17[[4]], grids_akpegg17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_akpegg3), fixed = T)
df_akpegg_avg3_miroc126 <- data.frame(lat = df_akpegg3$lat, 
                                      lon = df_akpegg3$lon, 
                                      doy = df_akpegg3$doy,
                                      avg_pred = rowSums(df_akpegg3[, x])/30)
saveRDS(df_akpegg_avg3_miroc126, file = here("data", "df_akpegg_avg3_miroc126_pheno.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akpegg_avg3_miroc126, "Forecasted doyribution 2070 - 2099 \n MIROC SSP126")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_miroc_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

##### MIROC 585 -----------------------------------------------------------------------------------------------------------------
temps_miroc_ssp585 <- readRDS(here('data', 'temps_miroc_ssp585.rds'))
salts_miroc_ssp585 <- readRDS(here('data', 'salts_miroc_ssp585.rds'))

## 2015 - 2039
grids_akpegg1 <- pred_loop_phenology(2015:2019, akp_egg, 137, 
                           temps_miroc_ssp585,
                           salts_miroc_ssp585, 1,
                           egg_formula)
grids_akpegg2 <- pred_loop_phenology(2020:2024, akp_egg, 137, 
                           temps_miroc_ssp585,
                           salts_miroc_ssp585, 2,
                           egg_formula)
grids_akpegg3 <- pred_loop_phenology(2025:2029, akp_egg, 137,
                           temps_miroc_ssp585,
                           salts_miroc_ssp585, 3,
                           egg_formula)
grids_akpegg4 <- pred_loop_phenology(2030:2034, akp_egg, 137, 
                           temps_miroc_ssp585,
                           salts_miroc_ssp585, 4,
                           egg_formula)
grids_akpegg5 <- pred_loop_phenology(2035:2039, akp_egg, 137, 
                           temps_miroc_ssp585,
                           salts_miroc_ssp585, 5,
                           egg_formula)

# Combine into one data frame
df_akpegg4 <- list(grids_akpegg1[[1]], grids_akpegg1[[2]], grids_akpegg1[[3]], 
                   grids_akpegg1[[4]], grids_akpegg1[[5]], grids_akpegg2[[1]], 
                   grids_akpegg2[[2]], grids_akpegg2[[3]], grids_akpegg2[[4]],
                   grids_akpegg2[[5]], grids_akpegg3[[1]], grids_akpegg3[[2]], 
                   grids_akpegg3[[3]], grids_akpegg3[[4]], grids_akpegg3[[5]],
                   grids_akpegg4[[1]], grids_akpegg4[[2]], grids_akpegg4[[3]], 
                   grids_akpegg4[[4]], grids_akpegg4[[5]], grids_akpegg5[[1]], 
                   grids_akpegg5[[2]], grids_akpegg5[[3]], grids_akpegg5[[4]],
                   grids_akpegg5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_akpegg4), fixed = T)
df_akpegg_avg4_miroc585 <- data.frame(lat = df_akpegg4$lat, 
                                      lon = df_akpegg4$lon, 
                                      doy = df_akpegg4$doy,
                                      avg_pred = rowSums(df_akpegg4[, x])/25)
saveRDS(df_akpegg_avg4_miroc585, file = here("data", "df_akpegg_avg4_miroc585_pheno.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akpegg_avg4_miroc585, "Forecasted doyribution 2015 - 2039 \n MIROC SSP585")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_miroc_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_akpegg6 <- pred_loop_phenology(2040:2044, akp_egg, 137, 
                           temps_miroc_ssp585,
                           salts_miroc_ssp585, 6,
                           egg_formula)
grids_akpegg7 <- pred_loop_phenology(2045:2049, akp_egg, 137, 
                           temps_miroc_ssp585,
                           salts_miroc_ssp585, 7,
                           egg_formula)
grids_akpegg8 <- pred_loop_phenology(2050:2054, akp_egg, 137, 
                           temps_miroc_ssp585,
                           salts_miroc_ssp585, 8,
                           egg_formula)
grids_akpegg9 <- pred_loop_phenology(2055:2059, akp_egg, 137, 
                           temps_miroc_ssp585,
                           salts_miroc_ssp585, 9,
                           egg_formula)
grids_akpegg10 <- pred_loop_phenology(2060:2064, akp_egg, 137, 
                            temps_miroc_ssp585,
                            salts_miroc_ssp585, 10,
                            egg_formula)
grids_akpegg11 <- pred_loop_phenology(2065:2069, akp_egg, 137, 
                            temps_miroc_ssp585,
                            salts_miroc_ssp585, 11,
                            egg_formula)

# Combine into one data frame
df_akpegg5 <- list(grids_akpegg6[[1]], grids_akpegg6[[2]], grids_akpegg6[[3]], 
                   grids_akpegg6[[4]], grids_akpegg6[[5]], grids_akpegg7[[1]], 
                   grids_akpegg7[[2]], grids_akpegg7[[3]], grids_akpegg7[[4]],
                   grids_akpegg7[[5]], grids_akpegg8[[1]], grids_akpegg8[[2]], 
                   grids_akpegg8[[3]], grids_akpegg8[[4]], grids_akpegg8[[5]],
                   grids_akpegg9[[1]], grids_akpegg9[[2]], grids_akpegg9[[3]], 
                   grids_akpegg9[[4]], grids_akpegg9[[5]], grids_akpegg10[[1]], 
                   grids_akpegg10[[2]], grids_akpegg10[[3]], grids_akpegg10[[4]],
                   grids_akpegg11[[5]], grids_akpegg11[[1]], grids_akpegg11[[2]],
                   grids_akpegg11[[3]], grids_akpegg11[[4]], grids_akpegg11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_akpegg5), fixed = T)
df_akpegg_avg5_miroc585 <- data.frame(lat = df_akpegg5$lat, 
                                      lon = df_akpegg5$lon, 
                                      doy = df_akpegg5$doy,
                                      avg_pred = rowSums(df_akpegg5[, x])/30)
saveRDS(df_akpegg_avg5_miroc585, file = here("data", "df_akpegg_avg5_miroc585_pheno.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akpegg_avg5_miroc585, "Forecasted doyribution 2040 - 2069 \n MIROC SSP585")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_miroc_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_akpegg12 <- pred_loop_phenology(2070:2074, akp_egg, 137,
                            temps_miroc_ssp585,
                            salts_miroc_ssp585, 12,
                            egg_formula)
grids_akpegg13 <- pred_loop_phenology(2075:2079, akp_egg, 137,
                            temps_miroc_ssp585,
                            salts_miroc_ssp585, 13,
                            egg_formula)
grids_akpegg14 <- pred_loop_phenology(2080:2084, akp_egg, 137,
                            temps_miroc_ssp585,
                            salts_miroc_ssp585, 14,
                            egg_formula)
grids_akpegg15 <- pred_loop_phenology(2085:2089, akp_egg, 137,
                            temps_miroc_ssp585,
                            salts_miroc_ssp585, 15,
                            egg_formula)
grids_akpegg16 <- pred_loop_phenology(2090:2094, akp_egg, 137,
                            temps_miroc_ssp585,
                            salts_miroc_ssp585, 16,
                            egg_formula)
grids_akpegg17 <- pred_loop_phenology(2095:2099, akp_egg, 137, 
                            temps_miroc_ssp585,
                            salts_miroc_ssp585, 17,
                            egg_formula)

# Combine into one data frame
df_akpegg6 <- list(grids_akpegg12[[1]], grids_akpegg12[[2]], grids_akpegg12[[3]], 
                   grids_akpegg12[[4]], grids_akpegg12[[5]], grids_akpegg13[[1]], 
                   grids_akpegg13[[2]], grids_akpegg13[[3]], grids_akpegg13[[4]],
                   grids_akpegg13[[5]], grids_akpegg14[[1]], grids_akpegg14[[2]], 
                   grids_akpegg14[[3]], grids_akpegg14[[4]], grids_akpegg14[[5]],
                   grids_akpegg15[[1]], grids_akpegg15[[2]], grids_akpegg15[[3]], 
                   grids_akpegg15[[4]], grids_akpegg15[[5]], grids_akpegg16[[1]], 
                   grids_akpegg16[[2]], grids_akpegg16[[3]], grids_akpegg16[[4]],
                   grids_akpegg17[[5]], grids_akpegg17[[1]], grids_akpegg17[[2]],
                   grids_akpegg17[[3]], grids_akpegg17[[4]], grids_akpegg17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_akpegg6), fixed = T)
df_akpegg_avg6_miroc585 <- data.frame(lat = df_akpegg6$lat, 
                                      lon = df_akpegg6$lon, 
                                      doy = df_akpegg6$doy,
                                      avg_pred = rowSums(df_akpegg6[, x])/30)
saveRDS(df_akpegg_avg6_miroc585, file = here("data", "df_akpegg_avg6_miroc585_pheno.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akpegg_avg6_miroc585, "Forecasted doyribution 2070 - 2099 \n MIROC SSP585")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_miroc_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()



### Plaice Larvae --------------------------------------------------------------------------------------------------------------------------
#### Forecast and average into 3 time periods ---------------------------------------------------------------------------------------------
##### CESM 126 ----------------------------------------------------------------------------------------------------------------------------
temps_cesm_ssp126 <- readRDS(here('data', 'temps_cesm_ssp126.rds'))
salts_cesm_ssp126 <- readRDS(here('data', 'salts_cesm_ssp126.rds'))

## 2015 - 2039
grids_akplarvae1 <- pred_loop_phenology(2015:2019, akp_larvae, 137, 
                              temps_cesm_ssp126,
                              salts_cesm_ssp126, 1,
                              larval_formula)
grids_akplarvae2 <- pred_loop_phenology(2020:2024, akp_larvae, 137, 
                              temps_cesm_ssp126,
                              salts_cesm_ssp126, 2,
                              larval_formula)
grids_akplarvae3 <- pred_loop_phenology(2025:2029, akp_larvae, 137,
                              temps_cesm_ssp126,
                              salts_cesm_ssp126, 3,
                              larval_formula)
grids_akplarvae4 <- pred_loop_phenology(2030:2034, akp_larvae, 137, 
                              temps_cesm_ssp126,
                              salts_cesm_ssp126, 4,
                              larval_formula)
grids_akplarvae5 <- pred_loop_phenology(2035:2039, akp_larvae, 137, 
                              temps_cesm_ssp126,
                              salts_cesm_ssp126, 5,
                              larval_formula)

# Combine into one data frame
df_akplarvae1 <- list(grids_akplarvae1[[1]], grids_akplarvae1[[2]], grids_akplarvae1[[3]], 
                      grids_akplarvae1[[4]], grids_akplarvae1[[5]], grids_akplarvae2[[1]], 
                      grids_akplarvae2[[2]], grids_akplarvae2[[3]], grids_akplarvae2[[4]],
                      grids_akplarvae2[[5]], grids_akplarvae3[[1]], grids_akplarvae3[[2]], 
                      grids_akplarvae3[[3]], grids_akplarvae3[[4]], grids_akplarvae3[[5]],
                      grids_akplarvae4[[1]], grids_akplarvae4[[2]], grids_akplarvae4[[3]], 
                      grids_akplarvae4[[4]], grids_akplarvae4[[5]], grids_akplarvae5[[1]], 
                      grids_akplarvae5[[2]], grids_akplarvae5[[3]], grids_akplarvae5[[4]],
                      grids_akplarvae5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_akplarvae1), fixed = T)
df_akplarvae_avg1_cesm126 <- data.frame(lat = df_akplarvae1$lat, 
                                        lon = df_akplarvae1$lon, 
                                        doy = df_akplarvae1$doy,
                                        avg_pred = rowSums(df_akplarvae1[, x])/25)
saveRDS(df_akplarvae_avg1_cesm126, file = here("data", "df_akplarvae_avg1_cesm126_pheno.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae_avg1_cesm126, "Forecasted doyribution 2015 - 2039 \n CESM SSP126")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_cesm_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_akplarvae6 <- pred_loop_phenology(2040:2044, akp_larvae, 137, 
                              temps_cesm_ssp126,
                              salts_cesm_ssp126, 6,
                              larval_formula)
grids_akplarvae7 <- pred_loop_phenology(2045:2049, akp_larvae, 137, 
                              temps_cesm_ssp126,
                              salts_cesm_ssp126, 7,
                              larval_formula)
grids_akplarvae8 <- pred_loop_phenology(2050:2054, akp_larvae, 137, 
                              temps_cesm_ssp126,
                              salts_cesm_ssp126, 8,
                              larval_formula)
grids_akplarvae9 <- pred_loop_phenology(2055:2059, akp_larvae, 137, 
                              temps_cesm_ssp126,
                              salts_cesm_ssp126, 9,
                              larval_formula)
grids_akplarvae10 <- pred_loop_phenology(2060:2064, akp_larvae, 137, 
                               temps_cesm_ssp126,
                               salts_cesm_ssp126, 10,
                               larval_formula)
grids_akplarvae11 <- pred_loop_phenology(2065:2069, akp_larvae, 137, 
                               temps_cesm_ssp126,
                               salts_cesm_ssp126, 11,
                               larval_formula)

# Combine into one data frame
df_akplarvae2 <- list(grids_akplarvae6[[1]], grids_akplarvae6[[2]], grids_akplarvae6[[3]], 
                      grids_akplarvae6[[4]], grids_akplarvae6[[5]], grids_akplarvae7[[1]], 
                      grids_akplarvae7[[2]], grids_akplarvae7[[3]], grids_akplarvae7[[4]],
                      grids_akplarvae7[[5]], grids_akplarvae8[[1]], grids_akplarvae8[[2]], 
                      grids_akplarvae8[[3]], grids_akplarvae8[[4]], grids_akplarvae8[[5]],
                      grids_akplarvae9[[1]], grids_akplarvae9[[2]], grids_akplarvae9[[3]], 
                      grids_akplarvae9[[4]], grids_akplarvae9[[5]], grids_akplarvae10[[1]], 
                      grids_akplarvae10[[2]], grids_akplarvae10[[3]], grids_akplarvae10[[4]],
                      grids_akplarvae11[[5]], grids_akplarvae11[[1]], grids_akplarvae11[[2]],
                      grids_akplarvae11[[3]], grids_akplarvae11[[4]], grids_akplarvae11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_akplarvae2), fixed = T)
df_akplarvae_avg2_cesm126 <- data.frame(lat = df_akplarvae2$lat, 
                                        lon = df_akplarvae2$lon, 
                                        doy = df_akplarvae2$doy,
                                        avg_pred = rowSums(df_akplarvae2[, x])/30)
saveRDS(df_akplarvae_avg2_cesm126, file = here("data", "df_akplarvae_avg2_cesm126_pheno.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae_avg2_cesm126, "Forecasted doyribution 2040 - 2069 \n CESM SSP126")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_cesm_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_akplarvae12 <- pred_loop_phenology(2070:2074, akp_larvae, 137,
                               temps_cesm_ssp126,
                               salts_cesm_ssp126, 12,
                               larval_formula)
grids_akplarvae13 <- pred_loop_phenology(2075:2079, akp_larvae, 137,
                               temps_cesm_ssp126,
                               salts_cesm_ssp126, 13,
                               larval_formula)
grids_akplarvae14 <- pred_loop_phenology(2080:2084, akp_larvae, 137,
                               temps_cesm_ssp126,
                               salts_cesm_ssp126, 14,
                               larval_formula)
grids_akplarvae15 <- pred_loop_phenology(2085:2089, akp_larvae, 137,
                               temps_cesm_ssp126,
                               salts_cesm_ssp126, 15,
                               larval_formula)
grids_akplarvae16 <- pred_loop_phenology(2090:2094, akp_larvae, 137,
                               temps_cesm_ssp126,
                               salts_cesm_ssp126, 16,
                               larval_formula)
grids_akplarvae17 <- pred_loop_phenology(2095:2099, akp_larvae, 137, 
                               temps_cesm_ssp126,
                               salts_cesm_ssp126, 17,
                               larval_formula)

# Combine into one data frame
df_akplarvae3 <- list(grids_akplarvae12[[1]], grids_akplarvae12[[2]], grids_akplarvae12[[3]], 
                      grids_akplarvae12[[4]], grids_akplarvae12[[5]], grids_akplarvae13[[1]], 
                      grids_akplarvae13[[2]], grids_akplarvae13[[3]], grids_akplarvae13[[4]],
                      grids_akplarvae13[[5]], grids_akplarvae14[[1]], grids_akplarvae14[[2]], 
                      grids_akplarvae14[[3]], grids_akplarvae14[[4]], grids_akplarvae14[[5]],
                      grids_akplarvae15[[1]], grids_akplarvae15[[2]], grids_akplarvae15[[3]], 
                      grids_akplarvae15[[4]], grids_akplarvae15[[5]], grids_akplarvae16[[1]], 
                      grids_akplarvae16[[2]], grids_akplarvae16[[3]], grids_akplarvae16[[4]],
                      grids_akplarvae17[[5]], grids_akplarvae17[[1]], grids_akplarvae17[[2]],
                      grids_akplarvae17[[3]], grids_akplarvae17[[4]], grids_akplarvae17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_akplarvae3), fixed = T)
df_akplarvae_avg3_cesm126 <- data.frame(lat = df_akplarvae3$lat, 
                                        lon = df_akplarvae3$lon, 
                                        doy = df_akplarvae3$doy,
                                        avg_pred = rowSums(df_akplarvae3[, x])/30)
saveRDS(df_akplarvae_avg3_cesm126, file = here("data", "df_akplarvae_avg3_cesm126_pheno.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae_avg3_cesm126, "Forecasted doyribution 2070 - 2099 \n CESM SSP126")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_cesm_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

##### CESM 585 -----------------------------------------------------------------------------------------------------------------
temps_cesm_ssp585 <- readRDS(here('data', 'temps_cesm_ssp585.rds'))
salts_cesm_ssp585 <- readRDS(here('data', 'salts_cesm_ssp585.rds'))

## 2015 - 2039
grids_akplarvae1 <- pred_loop_phenology(2015:2019, akp_larvae, 137, 
                              temps_cesm_ssp585,
                              salts_cesm_ssp585, 1,
                              larval_formula)
grids_akplarvae2 <- pred_loop_phenology(2020:2024, akp_larvae, 137, 
                              temps_cesm_ssp585,
                              salts_cesm_ssp585, 2,
                              larval_formula)
grids_akplarvae3 <- pred_loop_phenology(2025:2029, akp_larvae, 137,
                              temps_cesm_ssp585,
                              salts_cesm_ssp585, 3,
                              larval_formula)
grids_akplarvae4 <- pred_loop_phenology(2030:2034, akp_larvae, 137, 
                              temps_cesm_ssp585,
                              salts_cesm_ssp585, 4,
                              larval_formula)
grids_akplarvae5 <- pred_loop_phenology(2035:2039, akp_larvae, 137, 
                              temps_cesm_ssp585,
                              salts_cesm_ssp585, 5,
                              larval_formula)

# Combine into one data frame
df_akplarvae4 <- list(grids_akplarvae1[[1]], grids_akplarvae1[[2]], grids_akplarvae1[[3]], 
                      grids_akplarvae1[[4]], grids_akplarvae1[[5]], grids_akplarvae2[[1]], 
                      grids_akplarvae2[[2]], grids_akplarvae2[[3]], grids_akplarvae2[[4]],
                      grids_akplarvae2[[5]], grids_akplarvae3[[1]], grids_akplarvae3[[2]], 
                      grids_akplarvae3[[3]], grids_akplarvae3[[4]], grids_akplarvae3[[5]],
                      grids_akplarvae4[[1]], grids_akplarvae4[[2]], grids_akplarvae4[[3]], 
                      grids_akplarvae4[[4]], grids_akplarvae4[[5]], grids_akplarvae5[[1]], 
                      grids_akplarvae5[[2]], grids_akplarvae5[[3]], grids_akplarvae5[[4]],
                      grids_akplarvae5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_akplarvae4), fixed = T)
df_akplarvae_avg4_cesm585 <- data.frame(lat = df_akplarvae4$lat, 
                                        lon = df_akplarvae4$lon, 
                                        doy = df_akplarvae4$doy,
                                        avg_pred = rowSums(df_akplarvae4[, x])/25)
saveRDS(df_akplarvae_avg4_cesm585, file = here("data", "df_akplarvae_avg4_cesm585_pheno.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae_avg4_cesm585, "Forecasted doyribution 2015 - 2039 \n CESM SSP585")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_cesm_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_akplarvae6 <- pred_loop_phenology(2040:2044, akp_larvae, 137, 
                              temps_cesm_ssp585,
                              salts_cesm_ssp585, 6,
                              larval_formula)
grids_akplarvae7 <- pred_loop_phenology(2045:2049, akp_larvae, 137, 
                              temps_cesm_ssp585,
                              salts_cesm_ssp585, 7,
                              larval_formula)
grids_akplarvae8 <- pred_loop_phenology(2050:2054, akp_larvae, 137, 
                              temps_cesm_ssp585,
                              salts_cesm_ssp585, 8,
                              larval_formula)
grids_akplarvae9 <- pred_loop_phenology(2055:2059, akp_larvae, 137, 
                              temps_cesm_ssp585,
                              salts_cesm_ssp585, 9,
                              larval_formula)
grids_akplarvae10 <- pred_loop_phenology(2060:2064, akp_larvae, 137, 
                               temps_cesm_ssp585,
                               salts_cesm_ssp585, 10,
                               larval_formula)
grids_akplarvae11 <- pred_loop_phenology(2065:2069, akp_larvae, 137, 
                               temps_cesm_ssp585,
                               salts_cesm_ssp585, 11,
                               larval_formula)

# Combine into one data frame
df_akplarvae5 <- list(grids_akplarvae6[[1]], grids_akplarvae6[[2]], grids_akplarvae6[[3]], 
                      grids_akplarvae6[[4]], grids_akplarvae6[[5]], grids_akplarvae7[[1]], 
                      grids_akplarvae7[[2]], grids_akplarvae7[[3]], grids_akplarvae7[[4]],
                      grids_akplarvae7[[5]], grids_akplarvae8[[1]], grids_akplarvae8[[2]], 
                      grids_akplarvae8[[3]], grids_akplarvae8[[4]], grids_akplarvae8[[5]],
                      grids_akplarvae9[[1]], grids_akplarvae9[[2]], grids_akplarvae9[[3]], 
                      grids_akplarvae9[[4]], grids_akplarvae9[[5]], grids_akplarvae10[[1]], 
                      grids_akplarvae10[[2]], grids_akplarvae10[[3]], grids_akplarvae10[[4]],
                      grids_akplarvae11[[5]], grids_akplarvae11[[1]], grids_akplarvae11[[2]],
                      grids_akplarvae11[[3]], grids_akplarvae11[[4]], grids_akplarvae11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_akplarvae5), fixed = T)
df_akplarvae_avg5_cesm585 <- data.frame(lat = df_akplarvae5$lat, 
                                        lon = df_akplarvae5$lon, 
                                        doy = df_akplarvae5$doy,
                                        avg_pred = rowSums(df_akplarvae5[, x])/30)
saveRDS(df_akplarvae_avg5_cesm585, file = here("data", "df_akplarvae_avg5_cesm585_pheno.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae_avg5_cesm585, "Forecasted doyribution 2040 - 2069 \n CESM SSP585")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_cesm_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_akplarvae12 <- pred_loop_phenology(2070:2074, akp_larvae, 137,
                               temps_cesm_ssp585,
                               salts_cesm_ssp585, 12,
                               larval_formula)
grids_akplarvae13 <- pred_loop_phenology(2075:2079, akp_larvae, 137,
                               temps_cesm_ssp585,
                               salts_cesm_ssp585, 13,
                               larval_formula)
grids_akplarvae14 <- pred_loop_phenology(2080:2084, akp_larvae, 137,
                               temps_cesm_ssp585,
                               salts_cesm_ssp585, 14,
                               larval_formula)
grids_akplarvae15 <- pred_loop_phenology(2085:2089, akp_larvae, 137,
                               temps_cesm_ssp585,
                               salts_cesm_ssp585, 15,
                               larval_formula)
grids_akplarvae16 <- pred_loop_phenology(2090:2094, akp_larvae, 137,
                               temps_cesm_ssp585,
                               salts_cesm_ssp585, 16,
                               larval_formula)
grids_akplarvae17 <- pred_loop_phenology(2095:2099, akp_larvae, 137, 
                               temps_cesm_ssp585,
                               salts_cesm_ssp585, 17,
                               larval_formula)

# Combine into one data frame
df_akplarvae6 <- list(grids_akplarvae12[[1]], grids_akplarvae12[[2]], grids_akplarvae12[[3]], 
                      grids_akplarvae12[[4]], grids_akplarvae12[[5]], grids_akplarvae13[[1]], 
                      grids_akplarvae13[[2]], grids_akplarvae13[[3]], grids_akplarvae13[[4]],
                      grids_akplarvae13[[5]], grids_akplarvae14[[1]], grids_akplarvae14[[2]], 
                      grids_akplarvae14[[3]], grids_akplarvae14[[4]], grids_akplarvae14[[5]],
                      grids_akplarvae15[[1]], grids_akplarvae15[[2]], grids_akplarvae15[[3]], 
                      grids_akplarvae15[[4]], grids_akplarvae15[[5]], grids_akplarvae16[[1]], 
                      grids_akplarvae16[[2]], grids_akplarvae16[[3]], grids_akplarvae16[[4]],
                      grids_akplarvae17[[5]], grids_akplarvae17[[1]], grids_akplarvae17[[2]],
                      grids_akplarvae17[[3]], grids_akplarvae17[[4]], grids_akplarvae17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_akplarvae6), fixed = T)
df_akplarvae_avg6_cesm585 <- data.frame(lat = df_akplarvae6$lat, 
                                        lon = df_akplarvae6$lon, 
                                        doy = df_akplarvae6$doy,
                                        avg_pred = rowSums(df_akplarvae6[, x])/30)
saveRDS(df_akplarvae_avg6_cesm585, file = here("data", "df_akplarvae_avg6_cesm585_pheno.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae_avg6_cesm585, "Forecasted doyribution 2070 - 2099 \n CESM SSP585")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_cesm_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### GFDL 126 ----------------------------------------------------------------------------------------------------------------------------
temps_gfdl_ssp126 <- readRDS(here('data', 'temps_gfdl_ssp126.rds'))
salts_gfdl_ssp126 <- readRDS(here('data', 'salts_gfdl_ssp126.rds'))

## 2015 - 2039
grids_akplarvae1 <- pred_loop_phenology(2015:2019, akp_larvae, 137, 
                              temps_gfdl_ssp126,
                              salts_gfdl_ssp126, 1,
                              larval_formula)
grids_akplarvae2 <- pred_loop_phenology(2020:2024, akp_larvae, 137, 
                              temps_gfdl_ssp126,
                              salts_gfdl_ssp126, 2,
                              larval_formula)
grids_akplarvae3 <- pred_loop_phenology(2025:2029, akp_larvae, 137,
                              temps_gfdl_ssp126,
                              salts_gfdl_ssp126, 3,
                              larval_formula)
grids_akplarvae4 <- pred_loop_phenology(2030:2034, akp_larvae, 137, 
                              temps_gfdl_ssp126,
                              salts_gfdl_ssp126, 4,
                              larval_formula)
grids_akplarvae5 <- pred_loop_phenology(2035:2039, akp_larvae, 137, 
                              temps_gfdl_ssp126,
                              salts_gfdl_ssp126, 5,
                              larval_formula)

# Combine into one data frame
df_akplarvae1 <- list(grids_akplarvae1[[1]], grids_akplarvae1[[2]], grids_akplarvae1[[3]], 
                      grids_akplarvae1[[4]], grids_akplarvae1[[5]], grids_akplarvae2[[1]], 
                      grids_akplarvae2[[2]], grids_akplarvae2[[3]], grids_akplarvae2[[4]],
                      grids_akplarvae2[[5]], grids_akplarvae3[[1]], grids_akplarvae3[[2]], 
                      grids_akplarvae3[[3]], grids_akplarvae3[[4]], grids_akplarvae3[[5]],
                      grids_akplarvae4[[1]], grids_akplarvae4[[2]], grids_akplarvae4[[3]], 
                      grids_akplarvae4[[4]], grids_akplarvae4[[5]], grids_akplarvae5[[1]], 
                      grids_akplarvae5[[2]], grids_akplarvae5[[3]], grids_akplarvae5[[4]],
                      grids_akplarvae5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_akplarvae1), fixed = T)
df_akplarvae_avg1_gfdl126 <- data.frame(lat = df_akplarvae1$lat, 
                                        lon = df_akplarvae1$lon, 
                                        doy = df_akplarvae1$doy,
                                        avg_pred = rowSums(df_akplarvae1[, x])/25)
saveRDS(df_akplarvae_avg1_gfdl126, file = here("data", "df_akplarvae_avg1_gfdl126_pheno.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae_avg1_gfdl126, "Forecasted doyribution 2015 - 2039 \n GFDL SSP126")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_gfdl_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_akplarvae6 <- pred_loop_phenology(2040:2044, akp_larvae, 137, 
                              temps_gfdl_ssp126,
                              salts_gfdl_ssp126, 6,
                              larval_formula)
grids_akplarvae7 <- pred_loop_phenology(2045:2049, akp_larvae, 137, 
                              temps_gfdl_ssp126,
                              salts_gfdl_ssp126, 7,
                              larval_formula)
grids_akplarvae8 <- pred_loop_phenology(2050:2054, akp_larvae, 137, 
                              temps_gfdl_ssp126,
                              salts_gfdl_ssp126, 8,
                              larval_formula)
grids_akplarvae9 <- pred_loop_phenology(2055:2059, akp_larvae, 137, 
                              temps_gfdl_ssp126,
                              salts_gfdl_ssp126, 9,
                              larval_formula)
grids_akplarvae10 <- pred_loop_phenology(2060:2064, akp_larvae, 137, 
                               temps_gfdl_ssp126,
                               salts_gfdl_ssp126, 10,
                               larval_formula)
grids_akplarvae11 <- pred_loop_phenology(2065:2069, akp_larvae, 137, 
                               temps_gfdl_ssp126,
                               salts_gfdl_ssp126, 11,
                               larval_formula)

# Combine into one data frame
df_akplarvae2 <- list(grids_akplarvae6[[1]], grids_akplarvae6[[2]], grids_akplarvae6[[3]], 
                      grids_akplarvae6[[4]], grids_akplarvae6[[5]], grids_akplarvae7[[1]], 
                      grids_akplarvae7[[2]], grids_akplarvae7[[3]], grids_akplarvae7[[4]],
                      grids_akplarvae7[[5]], grids_akplarvae8[[1]], grids_akplarvae8[[2]], 
                      grids_akplarvae8[[3]], grids_akplarvae8[[4]], grids_akplarvae8[[5]],
                      grids_akplarvae9[[1]], grids_akplarvae9[[2]], grids_akplarvae9[[3]], 
                      grids_akplarvae9[[4]], grids_akplarvae9[[5]], grids_akplarvae10[[1]], 
                      grids_akplarvae10[[2]], grids_akplarvae10[[3]], grids_akplarvae10[[4]],
                      grids_akplarvae11[[5]], grids_akplarvae11[[1]], grids_akplarvae11[[2]],
                      grids_akplarvae11[[3]], grids_akplarvae11[[4]], grids_akplarvae11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_akplarvae2), fixed = T)
df_akplarvae_avg2_gfdl126 <- data.frame(lat = df_akplarvae2$lat, 
                                        lon = df_akplarvae2$lon, 
                                        doy = df_akplarvae2$doy,
                                        avg_pred = rowSums(df_akplarvae2[, x])/30)
saveRDS(df_akplarvae_avg2_gfdl126, file = here("data", "df_akplarvae_avg2_gfdl126_pheno.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae_avg2_gfdl126, "Forecasted doyribution 2040 - 2069 \n GFDL SSP126")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_gfdl_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_akplarvae12 <- pred_loop_phenology(2070:2074, akp_larvae, 137,
                               temps_gfdl_ssp126,
                               salts_gfdl_ssp126, 12,
                               larval_formula)
grids_akplarvae13 <- pred_loop_phenology(2075:2079, akp_larvae, 137,
                               temps_gfdl_ssp126,
                               salts_gfdl_ssp126, 13,
                               larval_formula)
grids_akplarvae14 <- pred_loop_phenology(2080:2084, akp_larvae, 137,
                               temps_gfdl_ssp126,
                               salts_gfdl_ssp126, 14,
                               larval_formula)
grids_akplarvae15 <- pred_loop_phenology(2085:2089, akp_larvae, 137,
                               temps_gfdl_ssp126,
                               salts_gfdl_ssp126, 15,
                               larval_formula)
grids_akplarvae16 <- pred_loop_phenology(2090:2094, akp_larvae, 137,
                               temps_gfdl_ssp126,
                               salts_gfdl_ssp126, 16,
                               larval_formula)
grids_akplarvae17 <- pred_loop_phenology(2095:2099, akp_larvae, 137, 
                               temps_gfdl_ssp126,
                               salts_gfdl_ssp126, 17,
                               larval_formula)

# Combine into one data frame
df_akplarvae3 <- list(grids_akplarvae12[[1]], grids_akplarvae12[[2]], grids_akplarvae12[[3]], 
                      grids_akplarvae12[[4]], grids_akplarvae12[[5]], grids_akplarvae13[[1]], 
                      grids_akplarvae13[[2]], grids_akplarvae13[[3]], grids_akplarvae13[[4]],
                      grids_akplarvae13[[5]], grids_akplarvae14[[1]], grids_akplarvae14[[2]], 
                      grids_akplarvae14[[3]], grids_akplarvae14[[4]], grids_akplarvae14[[5]],
                      grids_akplarvae15[[1]], grids_akplarvae15[[2]], grids_akplarvae15[[3]], 
                      grids_akplarvae15[[4]], grids_akplarvae15[[5]], grids_akplarvae16[[1]], 
                      grids_akplarvae16[[2]], grids_akplarvae16[[3]], grids_akplarvae16[[4]],
                      grids_akplarvae17[[5]], grids_akplarvae17[[1]], grids_akplarvae17[[2]],
                      grids_akplarvae17[[3]], grids_akplarvae17[[4]], grids_akplarvae17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_akplarvae3), fixed = T)
df_akplarvae_avg3_gfdl126 <- data.frame(lat = df_akplarvae3$lat, 
                                        lon = df_akplarvae3$lon, 
                                        doy = df_akplarvae3$doy,
                                        avg_pred = rowSums(df_akplarvae3[, x])/30)
saveRDS(df_akplarvae_avg3_gfdl126, file = here("data", "df_akplarvae_avg3_gfdl126_pheno.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae_avg3_gfdl126, "Forecasted doyribution 2070 - 2099 \n GFDL SSP126")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_gfdl_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

##### GFDL 585 -----------------------------------------------------------------------------------------------------------------
temps_gfdl_ssp585 <- readRDS(here('data', 'temps_gfdl_ssp585.rds'))
salts_gfdl_ssp585 <- readRDS(here('data', 'salts_gfdl_ssp585.rds'))

## 2015 - 2039
grids_akplarvae1 <- pred_loop_phenology(2015:2019, akp_larvae, 137, 
                              temps_gfdl_ssp585,
                              salts_gfdl_ssp585, 1,
                              larval_formula)
grids_akplarvae2 <- pred_loop_phenology(2020:2024, akp_larvae, 137, 
                              temps_gfdl_ssp585,
                              salts_gfdl_ssp585, 2,
                              larval_formula)
grids_akplarvae3 <- pred_loop_phenology(2025:2029, akp_larvae, 137,
                              temps_gfdl_ssp585,
                              salts_gfdl_ssp585, 3,
                              larval_formula)
grids_akplarvae4 <- pred_loop_phenology(2030:2034, akp_larvae, 137, 
                              temps_gfdl_ssp585,
                              salts_gfdl_ssp585, 4,
                              larval_formula)
grids_akplarvae5 <- pred_loop_phenology(2035:2039, akp_larvae, 137, 
                              temps_gfdl_ssp585,
                              salts_gfdl_ssp585, 5,
                              larval_formula)

# Combine into one data frame
df_akplarvae4 <- list(grids_akplarvae1[[1]], grids_akplarvae1[[2]], grids_akplarvae1[[3]], 
                      grids_akplarvae1[[4]], grids_akplarvae1[[5]], grids_akplarvae2[[1]], 
                      grids_akplarvae2[[2]], grids_akplarvae2[[3]], grids_akplarvae2[[4]],
                      grids_akplarvae2[[5]], grids_akplarvae3[[1]], grids_akplarvae3[[2]], 
                      grids_akplarvae3[[3]], grids_akplarvae3[[4]], grids_akplarvae3[[5]],
                      grids_akplarvae4[[1]], grids_akplarvae4[[2]], grids_akplarvae4[[3]], 
                      grids_akplarvae4[[4]], grids_akplarvae4[[5]], grids_akplarvae5[[1]], 
                      grids_akplarvae5[[2]], grids_akplarvae5[[3]], grids_akplarvae5[[4]],
                      grids_akplarvae5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_akplarvae4), fixed = T)
df_akplarvae_avg4_gfdl585 <- data.frame(lat = df_akplarvae4$lat, 
                                        lon = df_akplarvae4$lon, 
                                        doy = df_akplarvae4$doy,
                                        avg_pred = rowSums(df_akplarvae4[, x])/25)
saveRDS(df_akplarvae_avg4_gfdl585, file = here("data", "df_akplarvae_avg4_gfdl585_pheno.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae_avg4_gfdl585, "Forecasted doyribution 2015 - 2039 \n GFDL SSP585")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_gfdl_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_akplarvae6 <- pred_loop_phenology(2040:2044, akp_larvae, 137, 
                              temps_gfdl_ssp585,
                              salts_gfdl_ssp585, 6,
                              larval_formula)
grids_akplarvae7 <- pred_loop_phenology(2045:2049, akp_larvae, 137, 
                              temps_gfdl_ssp585,
                              salts_gfdl_ssp585, 7,
                              larval_formula)
grids_akplarvae8 <- pred_loop_phenology(2050:2054, akp_larvae, 137, 
                              temps_gfdl_ssp585,
                              salts_gfdl_ssp585, 8,
                              larval_formula)
grids_akplarvae9 <- pred_loop_phenology(2055:2059, akp_larvae, 137, 
                              temps_gfdl_ssp585,
                              salts_gfdl_ssp585, 9,
                              larval_formula)
grids_akplarvae10 <- pred_loop_phenology(2060:2064, akp_larvae, 137, 
                               temps_gfdl_ssp585,
                               salts_gfdl_ssp585, 10,
                               larval_formula)
grids_akplarvae11 <- pred_loop_phenology(2065:2069, akp_larvae, 137, 
                               temps_gfdl_ssp585,
                               salts_gfdl_ssp585, 11,
                               larval_formula)

# Combine into one data frame
df_akplarvae5 <- list(grids_akplarvae6[[1]], grids_akplarvae6[[2]], grids_akplarvae6[[3]], 
                      grids_akplarvae6[[4]], grids_akplarvae6[[5]], grids_akplarvae7[[1]], 
                      grids_akplarvae7[[2]], grids_akplarvae7[[3]], grids_akplarvae7[[4]],
                      grids_akplarvae7[[5]], grids_akplarvae8[[1]], grids_akplarvae8[[2]], 
                      grids_akplarvae8[[3]], grids_akplarvae8[[4]], grids_akplarvae8[[5]],
                      grids_akplarvae9[[1]], grids_akplarvae9[[2]], grids_akplarvae9[[3]], 
                      grids_akplarvae9[[4]], grids_akplarvae9[[5]], grids_akplarvae10[[1]], 
                      grids_akplarvae10[[2]], grids_akplarvae10[[3]], grids_akplarvae10[[4]],
                      grids_akplarvae11[[5]], grids_akplarvae11[[1]], grids_akplarvae11[[2]],
                      grids_akplarvae11[[3]], grids_akplarvae11[[4]], grids_akplarvae11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_akplarvae5), fixed = T)
df_akplarvae_avg5_gfdl585 <- data.frame(lat = df_akplarvae5$lat, 
                                        lon = df_akplarvae5$lon, 
                                        doy = df_akplarvae5$doy,
                                        avg_pred = rowSums(df_akplarvae5[, x])/30)
saveRDS(df_akplarvae_avg5_gfdl585, file = here("data", "df_akplarvae_avg5_gfdl585_pheno.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae_avg5_gfdl585, "Forecasted doyribution 2040 - 2069 \n GFDL SSP585")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_gfdl_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_akplarvae12 <- pred_loop_phenology(2070:2074, akp_larvae, 137,
                               temps_gfdl_ssp585,
                               salts_gfdl_ssp585, 12,
                               larval_formula)
grids_akplarvae13 <- pred_loop_phenology(2075:2079, akp_larvae, 137,
                               temps_gfdl_ssp585,
                               salts_gfdl_ssp585, 13,
                               larval_formula)
grids_akplarvae14 <- pred_loop_phenology(2080:2084, akp_larvae, 137,
                               temps_gfdl_ssp585,
                               salts_gfdl_ssp585, 14,
                               larval_formula)
grids_akplarvae15 <- pred_loop_phenology(2085:2089, akp_larvae, 137,
                               temps_gfdl_ssp585,
                               salts_gfdl_ssp585, 15,
                               larval_formula)
grids_akplarvae16 <- pred_loop_phenology(2090:2094, akp_larvae, 137,
                               temps_gfdl_ssp585,
                               salts_gfdl_ssp585, 16,
                               larval_formula)
grids_akplarvae17 <- pred_loop_phenology(2095:2099, akp_larvae, 137, 
                               temps_gfdl_ssp585,
                               salts_gfdl_ssp585, 17,
                               larval_formula)

# Combine into one data frame
df_akplarvae6 <- list(grids_akplarvae12[[1]], grids_akplarvae12[[2]], grids_akplarvae12[[3]], 
                      grids_akplarvae12[[4]], grids_akplarvae12[[5]], grids_akplarvae13[[1]], 
                      grids_akplarvae13[[2]], grids_akplarvae13[[3]], grids_akplarvae13[[4]],
                      grids_akplarvae13[[5]], grids_akplarvae14[[1]], grids_akplarvae14[[2]], 
                      grids_akplarvae14[[3]], grids_akplarvae14[[4]], grids_akplarvae14[[5]],
                      grids_akplarvae15[[1]], grids_akplarvae15[[2]], grids_akplarvae15[[3]], 
                      grids_akplarvae15[[4]], grids_akplarvae15[[5]], grids_akplarvae16[[1]], 
                      grids_akplarvae16[[2]], grids_akplarvae16[[3]], grids_akplarvae16[[4]],
                      grids_akplarvae17[[5]], grids_akplarvae17[[1]], grids_akplarvae17[[2]],
                      grids_akplarvae17[[3]], grids_akplarvae17[[4]], grids_akplarvae17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_akplarvae6), fixed = T)
df_akplarvae_avg6_gfdl585 <- data.frame(lat = df_akplarvae6$lat, 
                                        lon = df_akplarvae6$lon, 
                                        doy = df_akplarvae6$doy,
                                        avg_pred = rowSums(df_akplarvae6[, x])/30)
saveRDS(df_akplarvae_avg6_gfdl585, file = here("data", "df_akplarvae_avg6_gfdl585_pheno.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae_avg6_gfdl585, "Forecasted doyribution 2070 - 2099 \n GFDL SSP585")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_gfdl_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### MIROC 126 ----------------------------------------------------------------------------------------------------------------------------
temps_miroc_ssp126 <- readRDS(here('data', 'temps_miroc_ssp126.rds'))
salts_miroc_ssp126 <- readRDS(here('data', 'salts_miroc_ssp126.rds'))

## 2015 - 2039
grids_akplarvae1 <- pred_loop_phenology(2015:2019, akp_larvae, 137, 
                              temps_miroc_ssp126,
                              salts_miroc_ssp126, 1,
                              larval_formula)
grids_akplarvae2 <- pred_loop_phenology(2020:2024, akp_larvae, 137, 
                              temps_miroc_ssp126,
                              salts_miroc_ssp126, 2,
                              larval_formula)
grids_akplarvae3 <- pred_loop_phenology(2025:2029, akp_larvae, 137,
                              temps_miroc_ssp126,
                              salts_miroc_ssp126, 3,
                              larval_formula)
grids_akplarvae4 <- pred_loop_phenology(2030:2034, akp_larvae, 137, 
                              temps_miroc_ssp126,
                              salts_miroc_ssp126, 4,
                              larval_formula)
grids_akplarvae5 <- pred_loop_phenology(2035:2039, akp_larvae, 137, 
                              temps_miroc_ssp126,
                              salts_miroc_ssp126, 5,
                              larval_formula)

# Combine into one data frame
df_akplarvae1 <- list(grids_akplarvae1[[1]], grids_akplarvae1[[2]], grids_akplarvae1[[3]], 
                      grids_akplarvae1[[4]], grids_akplarvae1[[5]], grids_akplarvae2[[1]], 
                      grids_akplarvae2[[2]], grids_akplarvae2[[3]], grids_akplarvae2[[4]],
                      grids_akplarvae2[[5]], grids_akplarvae3[[1]], grids_akplarvae3[[2]], 
                      grids_akplarvae3[[3]], grids_akplarvae3[[4]], grids_akplarvae3[[5]],
                      grids_akplarvae4[[1]], grids_akplarvae4[[2]], grids_akplarvae4[[3]], 
                      grids_akplarvae4[[4]], grids_akplarvae4[[5]], grids_akplarvae5[[1]], 
                      grids_akplarvae5[[2]], grids_akplarvae5[[3]], grids_akplarvae5[[4]],
                      grids_akplarvae5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_akplarvae1), fixed = T)
df_akplarvae_avg1_miroc126 <- data.frame(lat = df_akplarvae1$lat, 
                                         lon = df_akplarvae1$lon, 
                                         doy = df_akplarvae1$doy,
                                         avg_pred = rowSums(df_akplarvae1[, x])/25)
saveRDS(df_akplarvae_avg1_miroc126, file = here("data", "df_akplarvae_avg1_miroc126_pheno.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae_avg1_miroc126, "Forecasted doyribution 2015 - 2039 \n MIROC SSP126")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_miroc_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_akplarvae6 <- pred_loop_phenology(2040:2044, akp_larvae, 137, 
                              temps_miroc_ssp126,
                              salts_miroc_ssp126, 6,
                              larval_formula)
grids_akplarvae7 <- pred_loop_phenology(2045:2049, akp_larvae, 137, 
                              temps_miroc_ssp126,
                              salts_miroc_ssp126, 7,
                              larval_formula)
grids_akplarvae8 <- pred_loop_phenology(2050:2054, akp_larvae, 137, 
                              temps_miroc_ssp126,
                              salts_miroc_ssp126, 8,
                              larval_formula)
grids_akplarvae9 <- pred_loop_phenology(2055:2059, akp_larvae, 137, 
                              temps_miroc_ssp126,
                              salts_miroc_ssp126, 9,
                              larval_formula)
grids_akplarvae10 <- pred_loop_phenology(2060:2064, akp_larvae, 137, 
                               temps_miroc_ssp126,
                               salts_miroc_ssp126, 10,
                               larval_formula)
grids_akplarvae11 <- pred_loop_phenology(2065:2069, akp_larvae, 137, 
                               temps_miroc_ssp126,
                               salts_miroc_ssp126, 11,
                               larval_formula)

# Combine into one data frame
df_akplarvae2 <- list(grids_akplarvae6[[1]], grids_akplarvae6[[2]], grids_akplarvae6[[3]], 
                      grids_akplarvae6[[4]], grids_akplarvae6[[5]], grids_akplarvae7[[1]], 
                      grids_akplarvae7[[2]], grids_akplarvae7[[3]], grids_akplarvae7[[4]],
                      grids_akplarvae7[[5]], grids_akplarvae8[[1]], grids_akplarvae8[[2]], 
                      grids_akplarvae8[[3]], grids_akplarvae8[[4]], grids_akplarvae8[[5]],
                      grids_akplarvae9[[1]], grids_akplarvae9[[2]], grids_akplarvae9[[3]], 
                      grids_akplarvae9[[4]], grids_akplarvae9[[5]], grids_akplarvae10[[1]], 
                      grids_akplarvae10[[2]], grids_akplarvae10[[3]], grids_akplarvae10[[4]],
                      grids_akplarvae11[[5]], grids_akplarvae11[[1]], grids_akplarvae11[[2]],
                      grids_akplarvae11[[3]], grids_akplarvae11[[4]], grids_akplarvae11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_akplarvae2), fixed = T)
df_akplarvae_avg2_miroc126 <- data.frame(lat = df_akplarvae2$lat, 
                                         lon = df_akplarvae2$lon, 
                                         doy = df_akplarvae2$doy,
                                         avg_pred = rowSums(df_akplarvae2[, x])/30)
saveRDS(df_akplarvae_avg2_miroc126, file = here("data", "df_akplarvae_avg2_miroc126_pheno.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae_avg2_miroc126, "Forecasted doyribution 2040 - 2069 \n MIROC SSP126")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_miroc_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_akplarvae12 <- pred_loop_phenology(2070:2074, akp_larvae, 137,
                               temps_miroc_ssp126,
                               salts_miroc_ssp126, 12,
                               larval_formula)
grids_akplarvae13 <- pred_loop_phenology(2075:2079, akp_larvae, 137,
                               temps_miroc_ssp126,
                               salts_miroc_ssp126, 13,
                               larval_formula)
grids_akplarvae14 <- pred_loop_phenology(2080:2084, akp_larvae, 137,
                               temps_miroc_ssp126,
                               salts_miroc_ssp126, 14,
                               larval_formula)
grids_akplarvae15 <- pred_loop_phenology(2085:2089, akp_larvae, 137,
                               temps_miroc_ssp126,
                               salts_miroc_ssp126, 15,
                               larval_formula)
grids_akplarvae16 <- pred_loop_phenology(2090:2094, akp_larvae, 137,
                               temps_miroc_ssp126,
                               salts_miroc_ssp126, 16,
                               larval_formula)
grids_akplarvae17 <- pred_loop_phenology(2095:2099, akp_larvae, 137, 
                               temps_miroc_ssp126,
                               salts_miroc_ssp126, 17,
                               larval_formula)

# Combine into one data frame
df_akplarvae3 <- list(grids_akplarvae12[[1]], grids_akplarvae12[[2]], grids_akplarvae12[[3]], 
                      grids_akplarvae12[[4]], grids_akplarvae12[[5]], grids_akplarvae13[[1]], 
                      grids_akplarvae13[[2]], grids_akplarvae13[[3]], grids_akplarvae13[[4]],
                      grids_akplarvae13[[5]], grids_akplarvae14[[1]], grids_akplarvae14[[2]], 
                      grids_akplarvae14[[3]], grids_akplarvae14[[4]], grids_akplarvae14[[5]],
                      grids_akplarvae15[[1]], grids_akplarvae15[[2]], grids_akplarvae15[[3]], 
                      grids_akplarvae15[[4]], grids_akplarvae15[[5]], grids_akplarvae16[[1]], 
                      grids_akplarvae16[[2]], grids_akplarvae16[[3]], grids_akplarvae16[[4]],
                      grids_akplarvae17[[5]], grids_akplarvae17[[1]], grids_akplarvae17[[2]],
                      grids_akplarvae17[[3]], grids_akplarvae17[[4]], grids_akplarvae17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_akplarvae3), fixed = T)
df_akplarvae_avg3_miroc126 <- data.frame(lat = df_akplarvae3$lat, 
                                         lon = df_akplarvae3$lon, 
                                         doy = df_akplarvae3$doy,
                                         avg_pred = rowSums(df_akplarvae3[, x])/30)
saveRDS(df_akplarvae_avg3_miroc126, file = here("data", "df_akplarvae_avg3_miroc126_pheno.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae_avg3_miroc126, "Forecasted doyribution 2070 - 2099 \n MIROC SSP126")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_miroc_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

##### MIROC 585 -----------------------------------------------------------------------------------------------------------------
temps_miroc_ssp585 <- readRDS(here('data', 'temps_miroc_ssp585.rds'))
salts_miroc_ssp585 <- readRDS(here('data', 'salts_miroc_ssp585.rds'))

## 2015 - 2039
grids_akplarvae1 <- pred_loop_phenology(2015:2019, akp_larvae, 137, 
                              temps_miroc_ssp585,
                              salts_miroc_ssp585, 1,
                              larval_formula)
grids_akplarvae2 <- pred_loop_phenology(2020:2024, akp_larvae, 137, 
                              temps_miroc_ssp585,
                              salts_miroc_ssp585, 2,
                              larval_formula)
grids_akplarvae3 <- pred_loop_phenology(2025:2029, akp_larvae, 137,
                              temps_miroc_ssp585,
                              salts_miroc_ssp585, 3,
                              larval_formula)
grids_akplarvae4 <- pred_loop_phenology(2030:2034, akp_larvae, 137, 
                              temps_miroc_ssp585,
                              salts_miroc_ssp585, 4,
                              larval_formula)
grids_akplarvae5 <- pred_loop_phenology(2035:2039, akp_larvae, 137, 
                              temps_miroc_ssp585,
                              salts_miroc_ssp585, 5,
                              larval_formula)

# Combine into one data frame
df_akplarvae4 <- list(grids_akplarvae1[[1]], grids_akplarvae1[[2]], grids_akplarvae1[[3]], 
                      grids_akplarvae1[[4]], grids_akplarvae1[[5]], grids_akplarvae2[[1]], 
                      grids_akplarvae2[[2]], grids_akplarvae2[[3]], grids_akplarvae2[[4]],
                      grids_akplarvae2[[5]], grids_akplarvae3[[1]], grids_akplarvae3[[2]], 
                      grids_akplarvae3[[3]], grids_akplarvae3[[4]], grids_akplarvae3[[5]],
                      grids_akplarvae4[[1]], grids_akplarvae4[[2]], grids_akplarvae4[[3]], 
                      grids_akplarvae4[[4]], grids_akplarvae4[[5]], grids_akplarvae5[[1]], 
                      grids_akplarvae5[[2]], grids_akplarvae5[[3]], grids_akplarvae5[[4]],
                      grids_akplarvae5[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_akplarvae4), fixed = T)
df_akplarvae_avg4_miroc585 <- data.frame(lat = df_akplarvae4$lat, 
                                         lon = df_akplarvae4$lon, 
                                         doy = df_akplarvae4$doy,
                                         avg_pred = rowSums(df_akplarvae4[, x])/25)
saveRDS(df_akplarvae_avg4_miroc585, file = here("data", "df_akplarvae_avg4_miroc585_pheno.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae_avg4_miroc585, "Forecasted doyribution 2015 - 2039 \n MIROC SSP585")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_miroc_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
grids_akplarvae6 <- pred_loop_phenology(2040:2044, akp_larvae, 137, 
                              temps_miroc_ssp585,
                              salts_miroc_ssp585, 6,
                              larval_formula)
grids_akplarvae7 <- pred_loop_phenology(2045:2049, akp_larvae, 137, 
                              temps_miroc_ssp585,
                              salts_miroc_ssp585, 7,
                              larval_formula)
grids_akplarvae8 <- pred_loop_phenology(2050:2054, akp_larvae, 137, 
                              temps_miroc_ssp585,
                              salts_miroc_ssp585, 8,
                              larval_formula)
grids_akplarvae9 <- pred_loop_phenology(2055:2059, akp_larvae, 137, 
                              temps_miroc_ssp585,
                              salts_miroc_ssp585, 9,
                              larval_formula)
grids_akplarvae10 <- pred_loop_phenology(2060:2064, akp_larvae, 137, 
                               temps_miroc_ssp585,
                               salts_miroc_ssp585, 10,
                               larval_formula)
grids_akplarvae11 <- pred_loop_phenology(2065:2069, akp_larvae, 137, 
                               temps_miroc_ssp585,
                               salts_miroc_ssp585, 11,
                               larval_formula)

# Combine into one data frame
df_akplarvae5 <- list(grids_akplarvae6[[1]], grids_akplarvae6[[2]], grids_akplarvae6[[3]], 
                      grids_akplarvae6[[4]], grids_akplarvae6[[5]], grids_akplarvae7[[1]], 
                      grids_akplarvae7[[2]], grids_akplarvae7[[3]], grids_akplarvae7[[4]],
                      grids_akplarvae7[[5]], grids_akplarvae8[[1]], grids_akplarvae8[[2]], 
                      grids_akplarvae8[[3]], grids_akplarvae8[[4]], grids_akplarvae8[[5]],
                      grids_akplarvae9[[1]], grids_akplarvae9[[2]], grids_akplarvae9[[3]], 
                      grids_akplarvae9[[4]], grids_akplarvae9[[5]], grids_akplarvae10[[1]], 
                      grids_akplarvae10[[2]], grids_akplarvae10[[3]], grids_akplarvae10[[4]],
                      grids_akplarvae11[[5]], grids_akplarvae11[[1]], grids_akplarvae11[[2]],
                      grids_akplarvae11[[3]], grids_akplarvae11[[4]], grids_akplarvae11[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_akplarvae5), fixed = T)
df_akplarvae_avg5_miroc585 <- data.frame(lat = df_akplarvae5$lat, 
                                         lon = df_akplarvae5$lon, 
                                         doy = df_akplarvae5$doy,
                                         avg_pred = rowSums(df_akplarvae5[, x])/30)
saveRDS(df_akplarvae_avg5_miroc585, file = here("data", "df_akplarvae_avg5_miroc585_pheno.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae_avg5_miroc585, "Forecasted doyribution 2040 - 2069 \n MIROC SSP585")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_miroc_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
grids_akplarvae12 <- pred_loop_phenology(2070:2074, akp_larvae, 137,
                               temps_miroc_ssp585,
                               salts_miroc_ssp585, 12,
                               larval_formula)
grids_akplarvae13 <- pred_loop_phenology(2075:2079, akp_larvae, 137,
                               temps_miroc_ssp585,
                               salts_miroc_ssp585, 13,
                               larval_formula)
grids_akplarvae14 <- pred_loop_phenology(2080:2084, akp_larvae, 137,
                               temps_miroc_ssp585,
                               salts_miroc_ssp585, 14,
                               larval_formula)
grids_akplarvae15 <- pred_loop_phenology(2085:2089, akp_larvae, 137,
                               temps_miroc_ssp585,
                               salts_miroc_ssp585, 15,
                               larval_formula)
grids_akplarvae16 <- pred_loop_phenology(2090:2094, akp_larvae, 137,
                               temps_miroc_ssp585,
                               salts_miroc_ssp585, 16,
                               larval_formula)
grids_akplarvae17 <- pred_loop_phenology(2095:2099, akp_larvae, 137, 
                               temps_miroc_ssp585,
                               salts_miroc_ssp585, 17,
                               larval_formula)

# Combine into one data frame
df_akplarvae6 <- list(grids_akplarvae12[[1]], grids_akplarvae12[[2]], grids_akplarvae12[[3]], 
                      grids_akplarvae12[[4]], grids_akplarvae12[[5]], grids_akplarvae13[[1]], 
                      grids_akplarvae13[[2]], grids_akplarvae13[[3]], grids_akplarvae13[[4]],
                      grids_akplarvae13[[5]], grids_akplarvae14[[1]], grids_akplarvae14[[2]], 
                      grids_akplarvae14[[3]], grids_akplarvae14[[4]], grids_akplarvae14[[5]],
                      grids_akplarvae15[[1]], grids_akplarvae15[[2]], grids_akplarvae15[[3]], 
                      grids_akplarvae15[[4]], grids_akplarvae15[[5]], grids_akplarvae16[[1]], 
                      grids_akplarvae16[[2]], grids_akplarvae16[[3]], grids_akplarvae16[[4]],
                      grids_akplarvae17[[5]], grids_akplarvae17[[1]], grids_akplarvae17[[2]],
                      grids_akplarvae17[[3]], grids_akplarvae17[[4]], grids_akplarvae17[[5]]) %>%
  reduce(inner_join, by = c("lon", "lat", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_akplarvae6), fixed = T)
df_akplarvae_avg6_miroc585 <- data.frame(lat = df_akplarvae6$lat, 
                                         lon = df_akplarvae6$lon, 
                                         doy = df_akplarvae6$doy,
                                         avg_pred = rowSums(df_akplarvae6[, x])/30)
saveRDS(df_akplarvae_avg6_miroc585, file = here("data", "df_akplarvae_avg6_miroc585_pheno.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae_avg6_miroc585, "Forecasted doyribution 2070 - 2099 \n MIROC SSP585")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_miroc_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

### Average Predictions ------------------------------------------------------------------------------------------------------------
#### 2015-2039 ---------------------------------------------------------------------------------------------------------------------
##### Eggs
df_akpegg_avg1_cesm126 <- readRDS(here('data', 'df_akpegg_avg1_cesm126_pheno.rds'))
df_akpegg_avg4_cesm585 <- readRDS(here('data', 'df_akpegg_avg4_cesm585_pheno.rds'))
df_akpegg_avg1_gfdl126 <- readRDS(here('data', 'df_akpegg_avg1_gfdl126_pheno.rds'))
df_akpegg_avg4_gfdl585 <- readRDS(here('data', 'df_akpegg_avg4_gfdl585_pheno.rds'))
df_akpegg_avg1_miroc126 <- readRDS(here('data', 'df_akpegg_avg1_miroc126_pheno.rds'))
df_akpegg_avg4_miroc585 <- readRDS(here('data', 'df_akpegg_avg4_miroc585_pheno.rds'))

df_akpegg_merged1 <- list(df_akpegg_avg1_cesm126, df_akpegg_avg4_cesm585,
                          df_akpegg_avg1_gfdl126, df_akpegg_avg4_gfdl585,
                          df_akpegg_avg1_miroc126, df_akpegg_avg4_miroc585) %>%
  reduce(inner_join, by = c("lon", "lat", "doy"))

x <- grepl("pred", names(df_akpegg_merged1), fixed = T)
df_akpegg_final1 <- data.frame(lat = df_akpegg_merged1$lat,
                               lon = df_akpegg_merged1$lon,
                               avg_pred = (rowSums(df_akpegg_merged1[, x])/6))

windows(width = 6, height = 6, family = "serif")
grid_predict(df_akpegg_final1, "Forecasted doyribution 2015 - 2039")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_avg1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

##### Larvae
df_akplarvae_avg1_cesm126 <- readRDS(here('data', 'df_akplarvae_avg1_cesm126_pheno.rds'))
df_akplarvae_avg4_cesm585 <- readRDS(here('data', 'df_akplarvae_avg4_cesm585_pheno.rds'))
df_akplarvae_avg1_gfdl126 <- readRDS(here('data', 'df_akplarvae_avg1_gfdl126_pheno.rds'))
df_akplarvae_avg4_gfdl585 <- readRDS(here('data', 'df_akplarvae_avg4_gfdl585_pheno.rds'))
df_akplarvae_avg1_miroc126 <- readRDS(here('data', 'df_akplarvae_avg1_miroc126_pheno.rds'))
df_akplarvae_avg4_miroc585 <- readRDS(here('data', 'df_akplarvae_avg4_miroc585_pheno.rds'))

df_akplarvae_merged1 <- list(df_akplarvae_avg1_cesm126, df_akplarvae_avg4_cesm585,
                             df_akplarvae_avg1_gfdl126, df_akplarvae_avg4_gfdl585,
                             df_akplarvae_avg1_miroc126, df_akplarvae_avg4_miroc585) %>%
  reduce(inner_join, by = c("lon", "lat", "doy"))

x <- grepl("pred", names(df_akplarvae_merged1), fixed = T)
df_akplarvae_final1 <- data.frame(lat = df_akplarvae_merged1$lat,
                                  lon = df_akplarvae_merged1$lon,
                                  avg_pred = (rowSums(df_akplarvae_merged1[, x])/6))

windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae_final1, "Forecasted doyribution 2015 - 2039")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_avg1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


#### 2040-2069 ---------------------------------------------------------------------------------------------------------------------
##### Eggs
df_akpegg_avg2_cesm126 <- readRDS(here('data', 'df_akpegg_avg2_cesm126_pheno.rds'))
df_akpegg_avg5_cesm585 <- readRDS(here('data', 'df_akpegg_avg5_cesm585_pheno.rds'))
df_akpegg_avg2_gfdl126 <- readRDS(here('data', 'df_akpegg_avg2_gfdl126_pheno.rds'))
df_akpegg_avg5_gfdl585 <- readRDS(here('data', 'df_akpegg_avg5_gfdl585_pheno.rds'))
df_akpegg_avg2_miroc126 <- readRDS(here('data', 'df_akpegg_avg2_miroc126_pheno.rds'))
df_akpegg_avg5_miroc585 <- readRDS(here('data', 'df_akpegg_avg5_miroc585_pheno.rds'))

df_akpegg_merged2 <- list(df_akpegg_avg2_cesm126, df_akpegg_avg5_cesm585,
                          df_akpegg_avg2_gfdl126, df_akpegg_avg5_gfdl585,
                          df_akpegg_avg2_miroc126, df_akpegg_avg5_miroc585) %>%
  reduce(inner_join, by = c("lon", "lat", "doy"))

x <- grepl("pred", names(df_akpegg_merged2), fixed = T)
df_akpegg_final2 <- data.frame(lat = df_akpegg_merged2$lat,
                               lon = df_akpegg_merged2$lon,
                               avg_pred = (rowSums(df_akpegg_merged2[, x])/6))

windows(width = 6, height = 6, family = "serif")
grid_predict(df_akpegg_final2, "Forecasted doyribution 2040 - 2069")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_avg2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

##### Larvae
df_akplarvae_avg2_cesm126 <- readRDS(here('data', 'df_akplarvae_avg2_cesm126_pheno.rds'))
df_akplarvae_avg5_cesm585 <- readRDS(here('data', 'df_akplarvae_avg5_cesm585_pheno.rds'))
df_akplarvae_avg2_gfdl126 <- readRDS(here('data', 'df_akplarvae_avg2_gfdl126_pheno.rds'))
df_akplarvae_avg5_gfdl585 <- readRDS(here('data', 'df_akplarvae_avg5_gfdl585_pheno.rds'))
df_akplarvae_avg2_miroc126 <- readRDS(here('data', 'df_akplarvae_avg2_miroc126_pheno.rds'))
df_akplarvae_avg5_miroc585 <- readRDS(here('data', 'df_akplarvae_avg5_miroc585_pheno.rds'))

df_akplarvae_merged2 <- list(df_akplarvae_avg2_cesm126, df_akplarvae_avg5_cesm585,
                             df_akplarvae_avg2_gfdl126, df_akplarvae_avg5_gfdl585,
                             df_akplarvae_avg2_miroc126, df_akplarvae_avg5_miroc585) %>%
  reduce(inner_join, by = c("lon", "lat", "doy"))

x <- grepl("pred", names(df_akplarvae_merged2), fixed = T)
df_akplarvae_final2 <- data.frame(lat = df_akplarvae_merged2$lat,
                                  lon = df_akplarvae_merged2$lon,
                                  avg_pred = (rowSums(df_akplarvae_merged2[, x])/6))

windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae_final2, "Forecasted doyribution 2040 - 2069")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_avg2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


#### 2070-2099----------------------------------------------------------------------------------------------------------------------
df_akpegg_avg3_cesm126 <- readRDS(here('data', 'df_akpegg_avg3_cesm126_pheno.rds'))
df_akpegg_avg6_cesm585 <- readRDS(here('data', 'df_akpegg_avg6_cesm585_pheno.rds'))
df_akpegg_avg3_gfdl126 <- readRDS(here('data', 'df_akpegg_avg3_gfdl126_pheno.rds'))
df_akpegg_avg6_gfdl585 <- readRDS(here('data', 'df_akpegg_avg6_gfdl585_pheno.rds'))
df_akpegg_avg3_miroc126 <- readRDS(here('data', 'df_akpegg_avg3_miroc126_pheno.rds'))
df_akpegg_avg6_miroc585 <- readRDS(here('data', 'df_akpegg_avg6_miroc585_pheno.rds'))

df_akpegg_merged3 <- list(df_akpegg_avg3_cesm126, df_akpegg_avg6_cesm585,
                          df_akpegg_avg3_gfdl126, df_akpegg_avg6_gfdl585,
                          df_akpegg_avg3_miroc126, df_akpegg_avg6_miroc585) %>%
  reduce(inner_join, by = c("lon", "lat", "doy"))

x <- grepl("pred", names(df_akpegg_merged3), fixed = T)
df_akpegg_final3 <- data.frame(lat = df_akpegg_merged3$lat,
                               lon = df_akpegg_merged3$lon,
                               avg_pred = (rowSums(df_akpegg_merged3[, x])/6))

windows(width = 6, height = 6, family = "serif")
grid_predict(df_akpegg_final3, "Forecasted doyribution 2070 - 2099")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_avg3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

df_akplarvae_avg3_cesm126 <- readRDS(here('data', 'df_akplarvae_avg3_cesm126_pheno.rds'))
df_akplarvae_avg6_cesm585 <- readRDS(here('data', 'df_akplarvae_avg6_cesm585_pheno.rds'))
df_akplarvae_avg3_gfdl126 <- readRDS(here('data', 'df_akplarvae_avg3_gfdl126_pheno.rds'))
df_akplarvae_avg6_gfdl585 <- readRDS(here('data', 'df_akplarvae_avg6_gfdl585_pheno.rds'))
df_akplarvae_avg3_miroc126 <- readRDS(here('data', 'df_akplarvae_avg3_miroc126_pheno.rds'))
df_akplarvae_avg6_miroc585 <- readRDS(here('data', 'df_akplarvae_avg6_miroc585_pheno.rds'))

df_akplarvae_merged3 <- list(df_akplarvae_avg3_cesm126, df_akplarvae_avg6_cesm585,
                             df_akplarvae_avg3_gfdl126, df_akplarvae_avg6_gfdl585,
                             df_akplarvae_avg3_miroc126, df_akplarvae_avg6_miroc585) %>%
  reduce(inner_join, by = c("lon", "lat", "doy"))

x <- grepl("pred", names(df_akplarvae_merged3), fixed = T)
df_akplarvae_final3 <- data.frame(lat = df_akplarvae_merged3$lat,
                                  lon = df_akplarvae_merged3$lon,
                                  avg_pred = (rowSums(df_akplarvae_merged3[, x])/6))

windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae_final3, "Forecasted doyribution 2070 - 2099")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_avg3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()