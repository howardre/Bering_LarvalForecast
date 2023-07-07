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

# Formulas
egg_formula <- gam(catch ~ s(year, bs = 're') +
                     s(doy, k = 8) +
                     s(lon, lat) +
                     s(roms_temperature, k = 6) +
                     s(roms_salinity, k = 6) +
                     s(doy, by = mean_temp, k = 6),
                   data = fhs_egg,
                   family = tw(link = 'log'),
                   method = 'REML')

larval_formula <- gam(catch ~ s(year, bs = 're') + 
                        s(doy, k = 8) +
                        s(lon, lat) +
                        s(roms_temperature, k = 6) +
                        s(roms_salinity, k = 6) +
                        s(doy, by = mean_temp, k = 6),
                      data = fhs_larvae,
                      family = tw(link = 'log'),
                      method = 'REML')

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
  #temporary
  # message(paste(the_year))
  # message(paste("nnidx max", max(grid_extent$nn.idx)))
  # message(paste("nnidx min", min(grid_extent$nn.idx)))
  # end temporary
  grid_extent$roms_temperature <- bc_temps[c(grid_extent$nn.idx), 10] # Match nearest temp
  grid_extent$roms_salinity <- bc_salts[c(grid_extent$nn.idx), 10] # Match nearest temp
  grid_extent <- grid_extent[-c(6:8)] # remove extra columns before predicting
  
  # Calculate mean temperature
  
  temp_filtered <- temp_output %>% filter(lon >= -170 & lon <= -165, 
                                          lat >= 56 & lat <= 58,
                                          month >= 2 & month <= 4,
                                          year == the_year)
  #temporary
  # message(paste('temp_output pre', nrow(temp_output)))
  # message(paste("temp filtered", nrow(temp_filtered)))
  # end temporary
  mean <- mean(temp_filtered$bc, na.rm = T)
  grid_extent$mean_temp <- mean
  
  # Parameterized model
  gam <- formula
  
  # Predict on forecasted output
  grid_extent$pred <- exp(predict(gam,
                                  newdata = grid_extent,
                                  type = "link",
                                  exclude = "s(year)"))
  grid_extent$pred[grid_extent$dist > 30000] <- NA
  return(grid_extent)
}


# New function to loop through years
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


## Use if any issues arise in get_preds function
# nlat = 40
# nlon = 60
# latd = seq(min(fhs_egg$lat), max(fhs_egg$lat), length.out = nlat)
# lond = seq(min(fhs_egg$lon), max(fhs_egg$lon), length.out = nlon)
# grid_extent <- expand.grid(lond, latd)
# names(grid_extent) <- c('lon', 'lat')
# 
# # Calculate distance of each grid point to closest 'positive observation'
# grid_extent$dist <- NA
# for (k in 1:nrow(grid_extent)) {
#   dist <- distance_function(grid_extent$lat[k],
#                             grid_extent$lon[k],
#                             fhs_egg$lat,
#                             fhs_egg$lon)
#   grid_extent$dist[k] <- min(dist)
# }

# Assign a within sample year and doy to the grid fhs_egg
# grid_extent$year <- 2040
# grid_extent$doy <- rep(130, length(grid_extent))
# grid_extent$month <- 5
# 
# cesm_temps2$lon <- -cesm_temps2$lon
# cesm_salts2$lon <- -cesm_salts2$lon
# # Use RANN package to match nearest temperature value based on lat, lon, month, year
# bc_temps <- cesm_temps2 %>% filter(month == 5 & year == 2040 & projection == 'ssp126')
# bc_salts <- cesm_salts2 %>% filter(month == 5 & year == 2040 & projection == 'ssp126')
# grid_extent[, c(7, 8)] <- as.data.frame(RANN::nn2(bc_temps[, c('lat', 'lon')],
#                                                   grid_extent[, c('lat', 'lon')],
#                                                   k = 1))
# grid_extent$roms_temperature <- bc_temps[c(grid_extent$nn.idx), 10] # Match nearest temp
# grid_extent$roms_salinity <- bc_salts[c(grid_extent$nn.idx), 10] # Match nearest temp
# grid_extent <- grid_extent[-c(6:8)] # remove extra columns before predicting
# temp_filtered <- cesm_temps1 %>% filter(lon >= -170 & lon <= -165, 
#                                         lat >= 56 & lat <= 58,
#                                         month >= 2 & month <= 4,
#                                         year == year)
# mean <- mean(temp_filtered$bc, na.rm = T)
# grid_extent$mean_temp <- mean
# gam <- egg_formula
# 
# # Predict on forecasted output
# grid_extent$pred <- exp(predict(gam,
#                                 newdata = grid_extent,
#                                 type = "link",
#                                 exclude = "s(year)"))
# grid_extent$pred[grid_extent$dist > 30000] <- NA


### Flathead Eggs --------------------------------------------------------------------------------------------------------------------------
#### Forecast and average into 3 time periods ---------------------------------------------------------------------------------------------
##### CESM 126 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
cesm_temps1 <- readRDS(here('data', 'cesm_forecast_temp1.rds'))
cesm_salts1 <- readRDS(here('data', 'cesm_forecast_salt1.rds'))

preds_fhsegg1_cesm126 <- pred_loop(2015:2039, fhs_egg, 130,
                                   5, 'ssp126', cesm_temps1, 
                                   cesm_salts1, egg_formula)

# Combine into one data frame
df_fhsegg1_cesm126 <- list(preds_fhsegg1_cesm126[[1]], preds_fhsegg1_cesm126[[2]],
                           preds_fhsegg1_cesm126[[3]], preds_fhsegg1_cesm126[[4]],
                           preds_fhsegg1_cesm126[[5]], preds_fhsegg1_cesm126[[6]],
                           preds_fhsegg1_cesm126[[7]], preds_fhsegg1_cesm126[[8]],
                           preds_fhsegg1_cesm126[[9]], preds_fhsegg1_cesm126[[10]],
                           preds_fhsegg1_cesm126[[11]], preds_fhsegg1_cesm126[[12]],
                           preds_fhsegg1_cesm126[[13]], preds_fhsegg1_cesm126[[14]],
                           preds_fhsegg1_cesm126[[15]], preds_fhsegg1_cesm126[[16]],
                           preds_fhsegg1_cesm126[[17]], preds_fhsegg1_cesm126[[18]],
                           preds_fhsegg1_cesm126[[19]], preds_fhsegg1_cesm126[[20]],
                           preds_fhsegg1_cesm126[[21]], preds_fhsegg1_cesm126[[22]],
                           preds_fhsegg1_cesm126[[23]], preds_fhsegg1_cesm126[[24]],
                           preds_fhsegg1_cesm126[[25]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhsegg1_cesm126), fixed = T)
df_fhsegg_avg1_cesm126 <- data.frame(lat = df_fhsegg1_cesm126$lat, 
                                     lon = df_fhsegg1_cesm126$lon, 
                                     dist = df_fhsegg1_cesm126$dist,
                                     avg_pred = rowSums(df_fhsegg1_cesm126[, x])/25)
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
cesm_temps2 <- readRDS(here('data', 'cesm_forecast_temp2.rds'))
cesm_salts2 <- readRDS(here('data', 'cesm_forecast_salt2.rds'))

preds_fhsegg2_cesm126 <- pred_loop(2040:2069, fhs_egg, 130,
                                   5, 'ssp126', cesm_temps2, 
                                   cesm_salts2, egg_formula)

# Combine into one data frame
df_fhsegg2_cesm126 <- list(preds_fhsegg2_cesm126[[1]], preds_fhsegg2_cesm126[[2]],
                           preds_fhsegg2_cesm126[[3]], preds_fhsegg2_cesm126[[4]],
                           preds_fhsegg2_cesm126[[5]], preds_fhsegg2_cesm126[[6]],
                           preds_fhsegg2_cesm126[[7]], preds_fhsegg2_cesm126[[8]],
                           preds_fhsegg2_cesm126[[9]], preds_fhsegg2_cesm126[[10]],
                           preds_fhsegg2_cesm126[[11]], preds_fhsegg2_cesm126[[12]],
                           preds_fhsegg2_cesm126[[13]], preds_fhsegg2_cesm126[[14]],
                           preds_fhsegg2_cesm126[[15]], preds_fhsegg2_cesm126[[16]],
                           preds_fhsegg2_cesm126[[17]], preds_fhsegg2_cesm126[[18]],
                           preds_fhsegg2_cesm126[[19]], preds_fhsegg2_cesm126[[20]],
                           preds_fhsegg2_cesm126[[21]], preds_fhsegg2_cesm126[[22]],
                           preds_fhsegg2_cesm126[[23]], preds_fhsegg2_cesm126[[24]],
                           preds_fhsegg2_cesm126[[25]], preds_fhsegg2_cesm126[[26]],
                           preds_fhsegg2_cesm126[[27]], preds_fhsegg2_cesm126[[28]],
                           preds_fhsegg2_cesm126[[29]], preds_fhsegg2_cesm126[[30]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhsegg2_cesm126), fixed = T)
df_fhsegg_avg2_cesm126 <- data.frame(lat = df_fhsegg2_cesm126$lat, 
                                     lon = df_fhsegg2_cesm126$lon, 
                                     dist = df_fhsegg2_cesm126$dist,
                                     avg_pred = rowSums(df_fhsegg2_cesm126[, x])/30)
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
cesm_temps3 <- readRDS(here('data', 'cesm_forecast_temp3.rds'))
cesm_salts3 <- readRDS(here('data', 'cesm_forecast_salt3.rds'))

preds_fhsegg3_cesm126 <- pred_loop(2070:2099, fhs_egg, 130,
                                   5, 'ssp126', cesm_temps3, 
                                   cesm_salts3, egg_formula)

# Combine into one data frame
df_fhsegg3_cesm126 <- list(preds_fhsegg3_cesm126[[1]], preds_fhsegg3_cesm126[[2]],
                           preds_fhsegg3_cesm126[[3]], preds_fhsegg3_cesm126[[4]],
                           preds_fhsegg3_cesm126[[5]], preds_fhsegg3_cesm126[[6]],
                           preds_fhsegg3_cesm126[[7]], preds_fhsegg3_cesm126[[8]],
                           preds_fhsegg3_cesm126[[9]], preds_fhsegg3_cesm126[[10]],
                           preds_fhsegg3_cesm126[[11]], preds_fhsegg3_cesm126[[12]],
                           preds_fhsegg3_cesm126[[13]], preds_fhsegg3_cesm126[[14]],
                           preds_fhsegg3_cesm126[[15]], preds_fhsegg3_cesm126[[16]],
                           preds_fhsegg3_cesm126[[17]], preds_fhsegg3_cesm126[[18]],
                           preds_fhsegg3_cesm126[[19]], preds_fhsegg3_cesm126[[20]],
                           preds_fhsegg3_cesm126[[21]], preds_fhsegg3_cesm126[[22]],
                           preds_fhsegg3_cesm126[[23]], preds_fhsegg3_cesm126[[24]],
                           preds_fhsegg3_cesm126[[25]], preds_fhsegg3_cesm126[[26]], 
                           preds_fhsegg3_cesm126[[27]], preds_fhsegg3_cesm126[[28]], 
                           preds_fhsegg3_cesm126[[29]], preds_fhsegg3_cesm126[[30]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhsegg3_cesm126), fixed = T)
df_fhsegg_avg3_cesm126 <- data.frame(lat = df_fhsegg3_cesm126$lat, 
                                     lon = df_fhsegg3_cesm126$lon, 
                                     dist = df_fhsegg3_cesm126$dist,
                                     avg_pred = rowSums(df_fhsegg3_cesm126[, x])/30)
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


rm(cesm_temps1, cesm_temps2, cesm_temps3,
   cesm_salts1, cesm_salts2, cesm_salts3)

##### CESM 585 -----------------------------------------------------------------------------------------------------------------
## 2015 - 2039
cesm_temps1 <- readRDS(here('data', 'cesm_forecast_temp1.rds'))
cesm_salts1 <- readRDS(here('data', 'cesm_forecast_salt1.rds'))

preds_fhsegg1_cesm585 <- pred_loop(2015:2039, fhs_egg, 130,
                                   5, 'ssp585', cesm_temps1, 
                                   cesm_salts1, egg_formula)

# Combine into one data frame
df_fhsegg1_cesm585 <- list(preds_fhsegg1_cesm585[[1]], preds_fhsegg1_cesm585[[2]],
                           preds_fhsegg1_cesm585[[3]], preds_fhsegg1_cesm585[[4]],
                           preds_fhsegg1_cesm585[[5]], preds_fhsegg1_cesm585[[6]],
                           preds_fhsegg1_cesm585[[7]], preds_fhsegg1_cesm585[[8]],
                           preds_fhsegg1_cesm585[[9]], preds_fhsegg1_cesm585[[10]],
                           preds_fhsegg1_cesm585[[11]], preds_fhsegg1_cesm585[[12]],
                           preds_fhsegg1_cesm585[[13]], preds_fhsegg1_cesm585[[14]],
                           preds_fhsegg1_cesm585[[15]], preds_fhsegg1_cesm585[[16]],
                           preds_fhsegg1_cesm585[[17]], preds_fhsegg1_cesm585[[18]],
                           preds_fhsegg1_cesm585[[19]], preds_fhsegg1_cesm585[[20]],
                           preds_fhsegg1_cesm585[[21]], preds_fhsegg1_cesm585[[22]],
                           preds_fhsegg1_cesm585[[23]], preds_fhsegg1_cesm585[[24]],
                           preds_fhsegg1_cesm585[[25]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhsegg1_cesm585), fixed = T)
df_fhsegg_avg1_cesm585 <- data.frame(lat = df_fhsegg1_cesm585$lat, 
                                     lon = df_fhsegg1_cesm585$lon, 
                                     dist = df_fhsegg1_cesm585$dist,
                                     avg_pred = rowSums(df_fhsegg1_cesm585[, x])/25)
saveRDS(df_fhsegg_avg1_cesm585, file = here("data", "df_fhsegg_avg1_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg_avg1_cesm585, "Forecasted Distribution 2015 - 2039 \n CESM SSP585")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_cesm_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
cesm_temps2 <- readRDS(here('data', 'cesm_forecast_temp2.rds'))
cesm_salts2 <- readRDS(here('data', 'cesm_forecast_salt2.rds'))

preds_fhsegg2_cesm585 <- pred_loop(2040:2069, fhs_egg, 130,
                                   5, 'ssp585', cesm_temps2, 
                                   cesm_salts2, egg_formula)

# Combine into one data frame
df_fhsegg2_cesm585 <- list(preds_fhsegg2_cesm585[[1]], preds_fhsegg2_cesm585[[2]],
                           preds_fhsegg2_cesm585[[3]], preds_fhsegg2_cesm585[[4]],
                           preds_fhsegg2_cesm585[[5]], preds_fhsegg2_cesm585[[6]],
                           preds_fhsegg2_cesm585[[7]], preds_fhsegg2_cesm585[[8]],
                           preds_fhsegg2_cesm585[[9]], preds_fhsegg2_cesm585[[10]],
                           preds_fhsegg2_cesm585[[11]], preds_fhsegg2_cesm585[[12]],
                           preds_fhsegg2_cesm585[[13]], preds_fhsegg2_cesm585[[14]],
                           preds_fhsegg2_cesm585[[15]], preds_fhsegg2_cesm585[[16]],
                           preds_fhsegg2_cesm585[[17]], preds_fhsegg2_cesm585[[18]],
                           preds_fhsegg2_cesm585[[19]], preds_fhsegg2_cesm585[[20]],
                           preds_fhsegg2_cesm585[[21]], preds_fhsegg2_cesm585[[22]],
                           preds_fhsegg2_cesm585[[23]], preds_fhsegg2_cesm585[[24]],
                           preds_fhsegg2_cesm585[[25]], preds_fhsegg2_cesm585[[26]],
                           preds_fhsegg2_cesm585[[27]], preds_fhsegg2_cesm585[[28]],
                           preds_fhsegg2_cesm585[[29]], preds_fhsegg2_cesm585[[30]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhsegg2_cesm585), fixed = T)
df_fhsegg_avg2_cesm585 <- data.frame(lat = df_fhsegg2_cesm585$lat, 
                                     lon = df_fhsegg2_cesm585$lon, 
                                     dist = df_fhsegg2_cesm585$dist,
                                     avg_pred = rowSums(df_fhsegg2_cesm585[, x])/30)
saveRDS(df_fhsegg_avg2_cesm585, file = here("data", "df_fhsegg_avg2_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg_avg2_cesm585, "Forecasted Distribution 2040 - 2069 \n CESM SSP585")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_cesm_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
cesm_temps3 <- readRDS(here('data', 'cesm_forecast_temp3.rds'))
cesm_salts3 <- readRDS(here('data', 'cesm_forecast_salt3.rds'))

preds_fhsegg3_cesm585 <- pred_loop(2070:2099, fhs_egg, 130,
                                   5, 'ssp585', cesm_temps3, 
                                   cesm_salts3, egg_formula)

# Combine into one data frame
df_fhsegg3_cesm585 <- list(preds_fhsegg3_cesm585[[1]], preds_fhsegg3_cesm585[[2]],
                           preds_fhsegg3_cesm585[[3]], preds_fhsegg3_cesm585[[4]],
                           preds_fhsegg3_cesm585[[5]], preds_fhsegg3_cesm585[[6]],
                           preds_fhsegg3_cesm585[[7]], preds_fhsegg3_cesm585[[8]],
                           preds_fhsegg3_cesm585[[9]], preds_fhsegg3_cesm585[[10]],
                           preds_fhsegg3_cesm585[[11]], preds_fhsegg3_cesm585[[12]],
                           preds_fhsegg3_cesm585[[13]], preds_fhsegg3_cesm585[[14]],
                           preds_fhsegg3_cesm585[[15]], preds_fhsegg3_cesm585[[16]],
                           preds_fhsegg3_cesm585[[17]], preds_fhsegg3_cesm585[[18]],
                           preds_fhsegg3_cesm585[[19]], preds_fhsegg3_cesm585[[20]],
                           preds_fhsegg3_cesm585[[21]], preds_fhsegg3_cesm585[[22]],
                           preds_fhsegg3_cesm585[[23]], preds_fhsegg3_cesm585[[24]],
                           preds_fhsegg3_cesm585[[25]], preds_fhsegg3_cesm585[[26]], 
                           preds_fhsegg3_cesm585[[27]], preds_fhsegg3_cesm585[[28]], 
                           preds_fhsegg3_cesm585[[29]], preds_fhsegg3_cesm585[[30]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhsegg3_cesm585), fixed = T)
df_fhsegg_avg3_cesm585 <- data.frame(lat = df_fhsegg3_cesm585$lat, 
                                     lon = df_fhsegg3_cesm585$lon, 
                                     dist = df_fhsegg3_cesm585$dist,
                                     avg_pred = rowSums(df_fhsegg3_cesm585[, x])/30)
saveRDS(df_fhsegg_avg3_cesm585, file = here("data", "df_fhsegg_avg3_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg_avg3_cesm585, "Forecasted Distribution 2070 - 2099 \n CESM SSP585")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_cesm_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


rm(cesm_temps1, cesm_temps2, cesm_temps3,
   cesm_salts1, cesm_salts2, cesm_salts3)

##### GFDL 126 ----------------------------------------------------------------------------------------------------------------------------
gfdl_temps1 <- readRDS(here('data', 'gfdl_forecast_temp1.rds'))
gfdl_salts1 <- readRDS(here('data', 'gfdl_forecast_salt1.rds'))

preds_fhsegg1_gfdl126 <- pred_loop(2015:2039, fhs_egg, 130,
                                   5, 'ssp126', gfdl_temps1, 
                                   gfdl_salts1, egg_formula)

# Combine into one data frame
df_fhsegg1_gfdl126 <- list(preds_fhsegg1_gfdl126[[1]], preds_fhsegg1_gfdl126[[2]],
                           preds_fhsegg1_gfdl126[[3]], preds_fhsegg1_gfdl126[[4]],
                           preds_fhsegg1_gfdl126[[5]], preds_fhsegg1_gfdl126[[6]],
                           preds_fhsegg1_gfdl126[[7]], preds_fhsegg1_gfdl126[[8]],
                           preds_fhsegg1_gfdl126[[9]], preds_fhsegg1_gfdl126[[10]],
                           preds_fhsegg1_gfdl126[[11]], preds_fhsegg1_gfdl126[[12]],
                           preds_fhsegg1_gfdl126[[13]], preds_fhsegg1_gfdl126[[14]],
                           preds_fhsegg1_gfdl126[[15]], preds_fhsegg1_gfdl126[[16]],
                           preds_fhsegg1_gfdl126[[17]], preds_fhsegg1_gfdl126[[18]],
                           preds_fhsegg1_gfdl126[[19]], preds_fhsegg1_gfdl126[[20]],
                           preds_fhsegg1_gfdl126[[21]], preds_fhsegg1_gfdl126[[22]],
                           preds_fhsegg1_gfdl126[[23]], preds_fhsegg1_gfdl126[[24]],
                           preds_fhsegg1_gfdl126[[25]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhsegg1_gfdl126), fixed = T)
df_fhsegg_avg1_gfdl126 <- data.frame(lat = df_fhsegg1_gfdl126$lat, 
                                     lon = df_fhsegg1_gfdl126$lon, 
                                     dist = df_fhsegg1_gfdl126$dist,
                                     avg_pred = rowSums(df_fhsegg1_gfdl126[, x])/25)
saveRDS(df_fhsegg_avg1_gfdl126, file = here("data", "df_fhsegg_avg1_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg_avg1_gfdl126, "Forecasted Distribution 2015 - 2039 \n gfdl SSP126")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_gfdl_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
gfdl_temps2 <- readRDS(here('data', 'gfdl_forecast_temp2.rds'))
gfdl_salts2 <- readRDS(here('data', 'gfdl_forecast_salt2.rds'))

preds_fhsegg2_gfdl126 <- pred_loop(2040:2069, fhs_egg, 130,
                                   5, 'ssp126', gfdl_temps2, 
                                   gfdl_salts2, egg_formula)

# Combine into one data frame
df_fhsegg2_gfdl126 <- list(preds_fhsegg2_gfdl126[[1]], preds_fhsegg2_gfdl126[[2]],
                           preds_fhsegg2_gfdl126[[3]], preds_fhsegg2_gfdl126[[4]],
                           preds_fhsegg2_gfdl126[[5]], preds_fhsegg2_gfdl126[[6]],
                           preds_fhsegg2_gfdl126[[7]], preds_fhsegg2_gfdl126[[8]],
                           preds_fhsegg2_gfdl126[[9]], preds_fhsegg2_gfdl126[[10]],
                           preds_fhsegg2_gfdl126[[11]], preds_fhsegg2_gfdl126[[12]],
                           preds_fhsegg2_gfdl126[[13]], preds_fhsegg2_gfdl126[[14]],
                           preds_fhsegg2_gfdl126[[15]], preds_fhsegg2_gfdl126[[16]],
                           preds_fhsegg2_gfdl126[[17]], preds_fhsegg2_gfdl126[[18]],
                           preds_fhsegg2_gfdl126[[19]], preds_fhsegg2_gfdl126[[20]],
                           preds_fhsegg2_gfdl126[[21]], preds_fhsegg2_gfdl126[[22]],
                           preds_fhsegg2_gfdl126[[23]], preds_fhsegg2_gfdl126[[24]],
                           preds_fhsegg2_gfdl126[[25]], preds_fhsegg2_gfdl126[[26]],
                           preds_fhsegg2_gfdl126[[27]], preds_fhsegg2_gfdl126[[28]],
                           preds_fhsegg2_gfdl126[[29]], preds_fhsegg2_gfdl126[[30]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhsegg2_gfdl126), fixed = T)
df_fhsegg_avg2_gfdl126 <- data.frame(lat = df_fhsegg2_gfdl126$lat, 
                                     lon = df_fhsegg2_gfdl126$lon, 
                                     dist = df_fhsegg2_gfdl126$dist,
                                     avg_pred = rowSums(df_fhsegg2_gfdl126[, x])/30)
saveRDS(df_fhsegg_avg2_gfdl126, file = here("data", "df_fhsegg_avg2_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg_avg2_gfdl126, "Forecasted Distribution 2040 - 2069 \n gfdl SSP126")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_gfdl_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
gfdl_temps3 <- readRDS(here('data', 'gfdl_forecast_temp3.rds'))
gfdl_salts3 <- readRDS(here('data', 'gfdl_forecast_salt3.rds'))

preds_fhsegg3_gfdl126 <- pred_loop(2070:2099, fhs_egg, 130,
                                   5, 'ssp126', gfdl_temps3, 
                                   gfdl_salts3, egg_formula)

# Combine into one data frame
df_fhsegg3_gfdl126 <- list(preds_fhsegg3_gfdl126[[1]], preds_fhsegg3_gfdl126[[2]],
                           preds_fhsegg3_gfdl126[[3]], preds_fhsegg3_gfdl126[[4]],
                           preds_fhsegg3_gfdl126[[5]], preds_fhsegg3_gfdl126[[6]],
                           preds_fhsegg3_gfdl126[[7]], preds_fhsegg3_gfdl126[[8]],
                           preds_fhsegg3_gfdl126[[9]], preds_fhsegg3_gfdl126[[10]],
                           preds_fhsegg3_gfdl126[[11]], preds_fhsegg3_gfdl126[[12]],
                           preds_fhsegg3_gfdl126[[13]], preds_fhsegg3_gfdl126[[14]],
                           preds_fhsegg3_gfdl126[[15]], preds_fhsegg3_gfdl126[[16]],
                           preds_fhsegg3_gfdl126[[17]], preds_fhsegg3_gfdl126[[18]],
                           preds_fhsegg3_gfdl126[[19]], preds_fhsegg3_gfdl126[[20]],
                           preds_fhsegg3_gfdl126[[21]], preds_fhsegg3_gfdl126[[22]],
                           preds_fhsegg3_gfdl126[[23]], preds_fhsegg3_gfdl126[[24]],
                           preds_fhsegg3_gfdl126[[25]], preds_fhsegg3_gfdl126[[26]], 
                           preds_fhsegg3_gfdl126[[27]], preds_fhsegg3_gfdl126[[28]], 
                           preds_fhsegg3_gfdl126[[29]], preds_fhsegg3_gfdl126[[30]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhsegg3_gfdl126), fixed = T)
df_fhsegg_avg3_gfdl126 <- data.frame(lat = df_fhsegg3_gfdl126$lat, 
                                     lon = df_fhsegg3_gfdl126$lon, 
                                     dist = df_fhsegg3_gfdl126$dist,
                                     avg_pred = rowSums(df_fhsegg3_gfdl126[, x])/30)
saveRDS(df_fhsegg_avg3_gfdl126, file = here("data", "df_fhsegg_avg3_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg_avg3_gfdl126, "Forecasted Distribution 2070 - 2099 \n gfdl SSP126")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_gfdl_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


rm(gfdl_temps1, gfdl_temps2, gfdl_temps3,
   gfdl_salts1, gfdl_salts2, gfdl_salts3)

##### GFDL 585 -----------------------------------------------------------------------------------------------------------------
## 2015 - 2039
gfdl_temps1 <- readRDS(here('data', 'gfdl_forecast_temp1.rds'))
gfdl_salts1 <- readRDS(here('data', 'gfdl_forecast_salt1.rds'))

preds_fhsegg1_gfdl585 <- pred_loop(2015:2039, fhs_egg, 130,
                                   5, 'ssp585', gfdl_temps1, 
                                   gfdl_salts1, egg_formula)

# Combine into one data frame
df_fhsegg1_gfdl585 <- list(preds_fhsegg1_gfdl585[[1]], preds_fhsegg1_gfdl585[[2]],
                           preds_fhsegg1_gfdl585[[3]], preds_fhsegg1_gfdl585[[4]],
                           preds_fhsegg1_gfdl585[[5]], preds_fhsegg1_gfdl585[[6]],
                           preds_fhsegg1_gfdl585[[7]], preds_fhsegg1_gfdl585[[8]],
                           preds_fhsegg1_gfdl585[[9]], preds_fhsegg1_gfdl585[[10]],
                           preds_fhsegg1_gfdl585[[11]], preds_fhsegg1_gfdl585[[12]],
                           preds_fhsegg1_gfdl585[[13]], preds_fhsegg1_gfdl585[[14]],
                           preds_fhsegg1_gfdl585[[15]], preds_fhsegg1_gfdl585[[16]],
                           preds_fhsegg1_gfdl585[[17]], preds_fhsegg1_gfdl585[[18]],
                           preds_fhsegg1_gfdl585[[19]], preds_fhsegg1_gfdl585[[20]],
                           preds_fhsegg1_gfdl585[[21]], preds_fhsegg1_gfdl585[[22]],
                           preds_fhsegg1_gfdl585[[23]], preds_fhsegg1_gfdl585[[24]],
                           preds_fhsegg1_gfdl585[[25]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhsegg1_gfdl585), fixed = T)
df_fhsegg_avg1_gfdl585 <- data.frame(lat = df_fhsegg1_gfdl585$lat, 
                                     lon = df_fhsegg1_gfdl585$lon, 
                                     dist = df_fhsegg1_gfdl585$dist,
                                     avg_pred = rowSums(df_fhsegg1_gfdl585[, x])/25)
saveRDS(df_fhsegg_avg1_gfdl585, file = here("data", "df_fhsegg_avg1_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg_avg1_gfdl585, "Forecasted Distribution 2015 - 2039 \n gfdl SSP585")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_gfdl_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
gfdl_temps2 <- readRDS(here('data', 'gfdl_forecast_temp2.rds'))
gfdl_salts2 <- readRDS(here('data', 'gfdl_forecast_salt2.rds'))

preds_fhsegg2_gfdl585 <- pred_loop(2040:2069, fhs_egg, 130,
                                   5, 'ssp585', gfdl_temps2, 
                                   gfdl_salts2, egg_formula)

# Combine into one data frame
df_fhsegg2_gfdl585 <- list(preds_fhsegg2_gfdl585[[1]], preds_fhsegg2_gfdl585[[2]],
                           preds_fhsegg2_gfdl585[[3]], preds_fhsegg2_gfdl585[[4]],
                           preds_fhsegg2_gfdl585[[5]], preds_fhsegg2_gfdl585[[6]],
                           preds_fhsegg2_gfdl585[[7]], preds_fhsegg2_gfdl585[[8]],
                           preds_fhsegg2_gfdl585[[9]], preds_fhsegg2_gfdl585[[10]],
                           preds_fhsegg2_gfdl585[[11]], preds_fhsegg2_gfdl585[[12]],
                           preds_fhsegg2_gfdl585[[13]], preds_fhsegg2_gfdl585[[14]],
                           preds_fhsegg2_gfdl585[[15]], preds_fhsegg2_gfdl585[[16]],
                           preds_fhsegg2_gfdl585[[17]], preds_fhsegg2_gfdl585[[18]],
                           preds_fhsegg2_gfdl585[[19]], preds_fhsegg2_gfdl585[[20]],
                           preds_fhsegg2_gfdl585[[21]], preds_fhsegg2_gfdl585[[22]],
                           preds_fhsegg2_gfdl585[[23]], preds_fhsegg2_gfdl585[[24]],
                           preds_fhsegg2_gfdl585[[25]], preds_fhsegg2_gfdl585[[26]],
                           preds_fhsegg2_gfdl585[[27]], preds_fhsegg2_gfdl585[[28]],
                           preds_fhsegg2_gfdl585[[29]], preds_fhsegg2_gfdl585[[30]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhsegg2_gfdl585), fixed = T)
df_fhsegg_avg2_gfdl585 <- data.frame(lat = df_fhsegg2_gfdl585$lat, 
                                     lon = df_fhsegg2_gfdl585$lon, 
                                     dist = df_fhsegg2_gfdl585$dist,
                                     avg_pred = rowSums(df_fhsegg2_gfdl585[, x])/30)
saveRDS(df_fhsegg_avg2_gfdl585, file = here("data", "df_fhsegg_avg2_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg_avg2_gfdl585, "Forecasted Distribution 2040 - 2069 \n gfdl SSP585")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_gfdl_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
gfdl_temps3 <- readRDS(here('data', 'gfdl_forecast_temp3.rds'))
gfdl_salts3 <- readRDS(here('data', 'gfdl_forecast_salt3.rds'))

preds_fhsegg3_gfdl585 <- pred_loop(2070:2099, fhs_egg, 130,
                                   5, 'ssp585', gfdl_temps3, 
                                   gfdl_salts3, egg_formula)

# Combine into one data frame
df_fhsegg3_gfdl585 <- list(preds_fhsegg3_gfdl585[[1]], preds_fhsegg3_gfdl585[[2]],
                           preds_fhsegg3_gfdl585[[3]], preds_fhsegg3_gfdl585[[4]],
                           preds_fhsegg3_gfdl585[[5]], preds_fhsegg3_gfdl585[[6]],
                           preds_fhsegg3_gfdl585[[7]], preds_fhsegg3_gfdl585[[8]],
                           preds_fhsegg3_gfdl585[[9]], preds_fhsegg3_gfdl585[[10]],
                           preds_fhsegg3_gfdl585[[11]], preds_fhsegg3_gfdl585[[12]],
                           preds_fhsegg3_gfdl585[[13]], preds_fhsegg3_gfdl585[[14]],
                           preds_fhsegg3_gfdl585[[15]], preds_fhsegg3_gfdl585[[16]],
                           preds_fhsegg3_gfdl585[[17]], preds_fhsegg3_gfdl585[[18]],
                           preds_fhsegg3_gfdl585[[19]], preds_fhsegg3_gfdl585[[20]],
                           preds_fhsegg3_gfdl585[[21]], preds_fhsegg3_gfdl585[[22]],
                           preds_fhsegg3_gfdl585[[23]], preds_fhsegg3_gfdl585[[24]],
                           preds_fhsegg3_gfdl585[[25]], preds_fhsegg3_gfdl585[[26]], 
                           preds_fhsegg3_gfdl585[[27]], preds_fhsegg3_gfdl585[[28]], 
                           preds_fhsegg3_gfdl585[[29]], preds_fhsegg3_gfdl585[[30]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhsegg3_gfdl585), fixed = T)
df_fhsegg_avg3_gfdl585 <- data.frame(lat = df_fhsegg3_gfdl585$lat, 
                                     lon = df_fhsegg3_gfdl585$lon, 
                                     dist = df_fhsegg3_gfdl585$dist,
                                     avg_pred = rowSums(df_fhsegg3_gfdl585[, x])/30)
saveRDS(df_fhsegg_avg3_gfdl585, file = here("data", "df_fhsegg_avg3_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg_avg3_gfdl585, "Forecasted Distribution 2070 - 2099 \n gfdl SSP585")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_gfdl_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


rm(gfdl_temps1, gfdl_temps2, gfdl_temps3,
   gfdl_salts1, gfdl_salts2, gfdl_salts3)


##### MIROC 126 ----------------------------------------------------------------------------------------------------------------------------
miroc_temps1 <- readRDS(here('data', 'miroc_forecast_temp1.rds'))
miroc_salts1 <- readRDS(here('data', 'miroc_forecast_salt1.rds'))

preds_fhsegg1_miroc126 <- pred_loop(2015:2039, fhs_egg, 130,
                                    5, 'ssp126', miroc_temps1, 
                                    miroc_salts1, egg_formula)

# Combine into one data frame
df_fhsegg1_miroc126 <- list(preds_fhsegg1_miroc126[[1]], preds_fhsegg1_miroc126[[2]],
                            preds_fhsegg1_miroc126[[3]], preds_fhsegg1_miroc126[[4]],
                            preds_fhsegg1_miroc126[[5]], preds_fhsegg1_miroc126[[6]],
                            preds_fhsegg1_miroc126[[7]], preds_fhsegg1_miroc126[[8]],
                            preds_fhsegg1_miroc126[[9]], preds_fhsegg1_miroc126[[10]],
                            preds_fhsegg1_miroc126[[11]], preds_fhsegg1_miroc126[[12]],
                            preds_fhsegg1_miroc126[[13]], preds_fhsegg1_miroc126[[14]],
                            preds_fhsegg1_miroc126[[15]], preds_fhsegg1_miroc126[[16]],
                            preds_fhsegg1_miroc126[[17]], preds_fhsegg1_miroc126[[18]],
                            preds_fhsegg1_miroc126[[19]], preds_fhsegg1_miroc126[[20]],
                            preds_fhsegg1_miroc126[[21]], preds_fhsegg1_miroc126[[22]],
                            preds_fhsegg1_miroc126[[23]], preds_fhsegg1_miroc126[[24]],
                            preds_fhsegg1_miroc126[[25]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhsegg1_miroc126), fixed = T)
df_fhsegg_avg1_miroc126 <- data.frame(lat = df_fhsegg1_miroc126$lat, 
                                      lon = df_fhsegg1_miroc126$lon, 
                                      dist = df_fhsegg1_miroc126$dist,
                                      avg_pred = rowSums(df_fhsegg1_miroc126[, x])/25)
saveRDS(df_fhsegg_avg1_miroc126, file = here("data", "df_fhsegg_avg1_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg_avg1_miroc126, "Forecasted Distribution 2015 - 2039 \n miroc SSP126")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_miroc_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
miroc_temps2 <- readRDS(here('data', 'miroc_forecast_temp2.rds'))
miroc_salts2 <- readRDS(here('data', 'miroc_forecast_salt2.rds'))

preds_fhsegg2_miroc126 <- pred_loop(2040:2069, fhs_egg, 130,
                                    5, 'ssp126', miroc_temps2, 
                                    miroc_salts2, egg_formula)

# Combine into one data frame
df_fhsegg2_miroc126 <- list(preds_fhsegg2_miroc126[[1]], preds_fhsegg2_miroc126[[2]],
                            preds_fhsegg2_miroc126[[3]], preds_fhsegg2_miroc126[[4]],
                            preds_fhsegg2_miroc126[[5]], preds_fhsegg2_miroc126[[6]],
                            preds_fhsegg2_miroc126[[7]], preds_fhsegg2_miroc126[[8]],
                            preds_fhsegg2_miroc126[[9]], preds_fhsegg2_miroc126[[10]],
                            preds_fhsegg2_miroc126[[11]], preds_fhsegg2_miroc126[[12]],
                            preds_fhsegg2_miroc126[[13]], preds_fhsegg2_miroc126[[14]],
                            preds_fhsegg2_miroc126[[15]], preds_fhsegg2_miroc126[[16]],
                            preds_fhsegg2_miroc126[[17]], preds_fhsegg2_miroc126[[18]],
                            preds_fhsegg2_miroc126[[19]], preds_fhsegg2_miroc126[[20]],
                            preds_fhsegg2_miroc126[[21]], preds_fhsegg2_miroc126[[22]],
                            preds_fhsegg2_miroc126[[23]], preds_fhsegg2_miroc126[[24]],
                            preds_fhsegg2_miroc126[[25]], preds_fhsegg2_miroc126[[26]],
                            preds_fhsegg2_miroc126[[27]], preds_fhsegg2_miroc126[[28]],
                            preds_fhsegg2_miroc126[[29]], preds_fhsegg2_miroc126[[30]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhsegg2_miroc126), fixed = T)
df_fhsegg_avg2_miroc126 <- data.frame(lat = df_fhsegg2_miroc126$lat, 
                                      lon = df_fhsegg2_miroc126$lon, 
                                      dist = df_fhsegg2_miroc126$dist,
                                      avg_pred = rowSums(df_fhsegg2_miroc126[, x])/30)
saveRDS(df_fhsegg_avg2_miroc126, file = here("data", "df_fhsegg_avg2_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg_avg2_miroc126, "Forecasted Distribution 2040 - 2069 \n miroc SSP126")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_miroc_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
miroc_temps3 <- readRDS(here('data', 'miroc_forecast_temp3.rds'))
miroc_salts3 <- readRDS(here('data', 'miroc_forecast_salt3.rds'))

preds_fhsegg3_miroc126 <- pred_loop(2070:2099, fhs_egg, 130,
                                    5, 'ssp126', miroc_temps3, 
                                    miroc_salts3, egg_formula)

# Combine into one data frame
df_fhsegg3_miroc126 <- list(preds_fhsegg3_miroc126[[1]], preds_fhsegg3_miroc126[[2]],
                            preds_fhsegg3_miroc126[[3]], preds_fhsegg3_miroc126[[4]],
                            preds_fhsegg3_miroc126[[5]], preds_fhsegg3_miroc126[[6]],
                            preds_fhsegg3_miroc126[[7]], preds_fhsegg3_miroc126[[8]],
                            preds_fhsegg3_miroc126[[9]], preds_fhsegg3_miroc126[[10]],
                            preds_fhsegg3_miroc126[[11]], preds_fhsegg3_miroc126[[12]],
                            preds_fhsegg3_miroc126[[13]], preds_fhsegg3_miroc126[[14]],
                            preds_fhsegg3_miroc126[[15]], preds_fhsegg3_miroc126[[16]],
                            preds_fhsegg3_miroc126[[17]], preds_fhsegg3_miroc126[[18]],
                            preds_fhsegg3_miroc126[[19]], preds_fhsegg3_miroc126[[20]],
                            preds_fhsegg3_miroc126[[21]], preds_fhsegg3_miroc126[[22]],
                            preds_fhsegg3_miroc126[[23]], preds_fhsegg3_miroc126[[24]],
                            preds_fhsegg3_miroc126[[25]], preds_fhsegg3_miroc126[[26]], 
                            preds_fhsegg3_miroc126[[27]], preds_fhsegg3_miroc126[[28]], 
                            preds_fhsegg3_miroc126[[29]], preds_fhsegg3_miroc126[[30]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhsegg3_miroc126), fixed = T)
df_fhsegg_avg3_miroc126 <- data.frame(lat = df_fhsegg3_miroc126$lat, 
                                      lon = df_fhsegg3_miroc126$lon, 
                                      dist = df_fhsegg3_miroc126$dist,
                                      avg_pred = rowSums(df_fhsegg3_miroc126[, x])/30)
saveRDS(df_fhsegg_avg3_miroc126, file = here("data", "df_fhsegg_avg3_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg_avg3_miroc126, "Forecasted Distribution 2070 - 2099 \n miroc SSP126")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_miroc_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


rm(miroc_temps1, miroc_temps2, miroc_temps3,
   miroc_salts1, miroc_salts2, miroc_salts3)

##### MIROC 585 -----------------------------------------------------------------------------------------------------------------
## 2015 - 2039
miroc_temps1 <- readRDS(here('data', 'miroc_forecast_temp1.rds'))
miroc_salts1 <- readRDS(here('data', 'miroc_forecast_salt1.rds'))

preds_fhsegg1_miroc585 <- pred_loop(2015:2039, fhs_egg, 130,
                                    5, 'ssp585', miroc_temps1, 
                                    miroc_salts1, egg_formula)

# Combine into one data frame
df_fhsegg1_miroc585 <- list(preds_fhsegg1_miroc585[[1]], preds_fhsegg1_miroc585[[2]],
                            preds_fhsegg1_miroc585[[3]], preds_fhsegg1_miroc585[[4]],
                            preds_fhsegg1_miroc585[[5]], preds_fhsegg1_miroc585[[6]],
                            preds_fhsegg1_miroc585[[7]], preds_fhsegg1_miroc585[[8]],
                            preds_fhsegg1_miroc585[[9]], preds_fhsegg1_miroc585[[10]],
                            preds_fhsegg1_miroc585[[11]], preds_fhsegg1_miroc585[[12]],
                            preds_fhsegg1_miroc585[[13]], preds_fhsegg1_miroc585[[14]],
                            preds_fhsegg1_miroc585[[15]], preds_fhsegg1_miroc585[[16]],
                            preds_fhsegg1_miroc585[[17]], preds_fhsegg1_miroc585[[18]],
                            preds_fhsegg1_miroc585[[19]], preds_fhsegg1_miroc585[[20]],
                            preds_fhsegg1_miroc585[[21]], preds_fhsegg1_miroc585[[22]],
                            preds_fhsegg1_miroc585[[23]], preds_fhsegg1_miroc585[[24]],
                            preds_fhsegg1_miroc585[[25]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhsegg1_miroc585), fixed = T)
df_fhsegg_avg1_miroc585 <- data.frame(lat = df_fhsegg1_miroc585$lat, 
                                      lon = df_fhsegg1_miroc585$lon, 
                                      dist = df_fhsegg1_miroc585$dist,
                                      avg_pred = rowSums(df_fhsegg1_miroc585[, x])/25)
saveRDS(df_fhsegg_avg1_miroc585, file = here("data", "df_fhsegg_avg1_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg_avg1_miroc585, "Forecasted Distribution 2015 - 2039 \n miroc SSP585")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_miroc_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
miroc_temps2 <- readRDS(here('data', 'miroc_forecast_temp2.rds'))
miroc_salts2 <- readRDS(here('data', 'miroc_forecast_salt2.rds'))

preds_fhsegg2_miroc585 <- pred_loop(2040:2069, fhs_egg, 130,
                                    5, 'ssp585', miroc_temps2, 
                                    miroc_salts2, egg_formula)

# Combine into one data frame
df_fhsegg2_miroc585 <- list(preds_fhsegg2_miroc585[[1]], preds_fhsegg2_miroc585[[2]],
                            preds_fhsegg2_miroc585[[3]], preds_fhsegg2_miroc585[[4]],
                            preds_fhsegg2_miroc585[[5]], preds_fhsegg2_miroc585[[6]],
                            preds_fhsegg2_miroc585[[7]], preds_fhsegg2_miroc585[[8]],
                            preds_fhsegg2_miroc585[[9]], preds_fhsegg2_miroc585[[10]],
                            preds_fhsegg2_miroc585[[11]], preds_fhsegg2_miroc585[[12]],
                            preds_fhsegg2_miroc585[[13]], preds_fhsegg2_miroc585[[14]],
                            preds_fhsegg2_miroc585[[15]], preds_fhsegg2_miroc585[[16]],
                            preds_fhsegg2_miroc585[[17]], preds_fhsegg2_miroc585[[18]],
                            preds_fhsegg2_miroc585[[19]], preds_fhsegg2_miroc585[[20]],
                            preds_fhsegg2_miroc585[[21]], preds_fhsegg2_miroc585[[22]],
                            preds_fhsegg2_miroc585[[23]], preds_fhsegg2_miroc585[[24]],
                            preds_fhsegg2_miroc585[[25]], preds_fhsegg2_miroc585[[26]],
                            preds_fhsegg2_miroc585[[27]], preds_fhsegg2_miroc585[[28]],
                            preds_fhsegg2_miroc585[[29]], preds_fhsegg2_miroc585[[30]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhsegg2_miroc585), fixed = T)
df_fhsegg_avg2_miroc585 <- data.frame(lat = df_fhsegg2_miroc585$lat, 
                                      lon = df_fhsegg2_miroc585$lon, 
                                      dist = df_fhsegg2_miroc585$dist,
                                      avg_pred = rowSums(df_fhsegg2_miroc585[, x])/30)
saveRDS(df_fhsegg_avg2_miroc585, file = here("data", "df_fhsegg_avg2_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg_avg2_miroc585, "Forecasted Distribution 2040 - 2069 \n miroc SSP585")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_miroc_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
miroc_temps3 <- readRDS(here('data', 'miroc_forecast_temp3.rds'))
miroc_salts3 <- readRDS(here('data', 'miroc_forecast_salt3.rds'))

preds_fhsegg3_miroc585 <- pred_loop(2070:2099, fhs_egg, 130,
                                    5, 'ssp585', miroc_temps3, 
                                    miroc_salts3, egg_formula)

# Combine into one data frame
df_fhsegg3_miroc585 <- list(preds_fhsegg3_miroc585[[1]], preds_fhsegg3_miroc585[[2]],
                            preds_fhsegg3_miroc585[[3]], preds_fhsegg3_miroc585[[4]],
                            preds_fhsegg3_miroc585[[5]], preds_fhsegg3_miroc585[[6]],
                            preds_fhsegg3_miroc585[[7]], preds_fhsegg3_miroc585[[8]],
                            preds_fhsegg3_miroc585[[9]], preds_fhsegg3_miroc585[[10]],
                            preds_fhsegg3_miroc585[[11]], preds_fhsegg3_miroc585[[12]],
                            preds_fhsegg3_miroc585[[13]], preds_fhsegg3_miroc585[[14]],
                            preds_fhsegg3_miroc585[[15]], preds_fhsegg3_miroc585[[16]],
                            preds_fhsegg3_miroc585[[17]], preds_fhsegg3_miroc585[[18]],
                            preds_fhsegg3_miroc585[[19]], preds_fhsegg3_miroc585[[20]],
                            preds_fhsegg3_miroc585[[21]], preds_fhsegg3_miroc585[[22]],
                            preds_fhsegg3_miroc585[[23]], preds_fhsegg3_miroc585[[24]],
                            preds_fhsegg3_miroc585[[25]], preds_fhsegg3_miroc585[[26]], 
                            preds_fhsegg3_miroc585[[27]], preds_fhsegg3_miroc585[[28]], 
                            preds_fhsegg3_miroc585[[29]], preds_fhsegg3_miroc585[[30]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhsegg3_miroc585), fixed = T)
df_fhsegg_avg3_miroc585 <- data.frame(lat = df_fhsegg3_miroc585$lat, 
                                      lon = df_fhsegg3_miroc585$lon, 
                                      dist = df_fhsegg3_miroc585$dist,
                                      avg_pred = rowSums(df_fhsegg3_miroc585[, x])/30)
saveRDS(df_fhsegg_avg3_miroc585, file = here("data", "df_fhsegg_avg3_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg_avg3_miroc585, "Forecasted Distribution 2070 - 2099 \n miroc SSP585")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_miroc_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


rm(miroc_temps1, miroc_temps2, miroc_temps3,
   miroc_salts1, miroc_salts2, miroc_salts3)


### Flathead Larvae --------------------------------------------------------------------------------------------------------------------------
#### Forecast and average into 3 time periods ---------------------------------------------------------------------------------------------
##### CESM 126 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
cesm_temps1 <- readRDS(here('data', 'cesm_forecast_temp1.rds'))
cesm_salts1 <- readRDS(here('data', 'cesm_forecast_salt1.rds'))

preds_fhslarvae1_cesm126 <- pred_loop(2015:2039, fhs_larvae, 130,
                                      5, 'ssp126', cesm_temps1, 
                                      cesm_salts1, larval_formula)

# Combine into one data frame
df_fhslarvae1_cesm126 <- list(preds_fhslarvae1_cesm126[[1]], preds_fhslarvae1_cesm126[[2]],
                              preds_fhslarvae1_cesm126[[3]], preds_fhslarvae1_cesm126[[4]],
                              preds_fhslarvae1_cesm126[[5]], preds_fhslarvae1_cesm126[[6]],
                              preds_fhslarvae1_cesm126[[7]], preds_fhslarvae1_cesm126[[8]],
                              preds_fhslarvae1_cesm126[[9]], preds_fhslarvae1_cesm126[[10]],
                              preds_fhslarvae1_cesm126[[11]], preds_fhslarvae1_cesm126[[12]],
                              preds_fhslarvae1_cesm126[[13]], preds_fhslarvae1_cesm126[[14]],
                              preds_fhslarvae1_cesm126[[15]], preds_fhslarvae1_cesm126[[16]],
                              preds_fhslarvae1_cesm126[[17]], preds_fhslarvae1_cesm126[[18]],
                              preds_fhslarvae1_cesm126[[19]], preds_fhslarvae1_cesm126[[20]],
                              preds_fhslarvae1_cesm126[[21]], preds_fhslarvae1_cesm126[[22]],
                              preds_fhslarvae1_cesm126[[23]], preds_fhslarvae1_cesm126[[24]],
                              preds_fhslarvae1_cesm126[[25]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhslarvae1_cesm126), fixed = T)
df_fhslarvae_avg1_cesm126 <- data.frame(lat = df_fhslarvae1_cesm126$lat, 
                                        lon = df_fhslarvae1_cesm126$lon, 
                                        dist = df_fhslarvae1_cesm126$dist,
                                        avg_pred = rowSums(df_fhslarvae1_cesm126[, x])/25)
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
cesm_temps2 <- readRDS(here('data', 'cesm_forecast_temp2.rds'))
cesm_salts2 <- readRDS(here('data', 'cesm_forecast_salt2.rds'))

preds_fhslarvae2_cesm126 <- pred_loop(2040:2069, fhs_larvae, 130,
                                      5, 'ssp126', cesm_temps2, 
                                      cesm_salts2, larval_formula)

# Combine into one data frame
df_fhslarvae2_cesm126 <- list(preds_fhslarvae2_cesm126[[1]], preds_fhslarvae2_cesm126[[2]],
                              preds_fhslarvae2_cesm126[[3]], preds_fhslarvae2_cesm126[[4]],
                              preds_fhslarvae2_cesm126[[5]], preds_fhslarvae2_cesm126[[6]],
                              preds_fhslarvae2_cesm126[[7]], preds_fhslarvae2_cesm126[[8]],
                              preds_fhslarvae2_cesm126[[9]], preds_fhslarvae2_cesm126[[10]],
                              preds_fhslarvae2_cesm126[[11]], preds_fhslarvae2_cesm126[[12]],
                              preds_fhslarvae2_cesm126[[13]], preds_fhslarvae2_cesm126[[14]],
                              preds_fhslarvae2_cesm126[[15]], preds_fhslarvae2_cesm126[[16]],
                              preds_fhslarvae2_cesm126[[17]], preds_fhslarvae2_cesm126[[18]],
                              preds_fhslarvae2_cesm126[[19]], preds_fhslarvae2_cesm126[[20]],
                              preds_fhslarvae2_cesm126[[21]], preds_fhslarvae2_cesm126[[22]],
                              preds_fhslarvae2_cesm126[[23]], preds_fhslarvae2_cesm126[[24]],
                              preds_fhslarvae2_cesm126[[25]], preds_fhslarvae2_cesm126[[26]],
                              preds_fhslarvae2_cesm126[[27]], preds_fhslarvae2_cesm126[[28]],
                              preds_fhslarvae2_cesm126[[29]], preds_fhslarvae2_cesm126[[30]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhslarvae2_cesm126), fixed = T)
df_fhslarvae_avg2_cesm126 <- data.frame(lat = df_fhslarvae2_cesm126$lat, 
                                        lon = df_fhslarvae2_cesm126$lon, 
                                        dist = df_fhslarvae2_cesm126$dist,
                                        avg_pred = rowSums(df_fhslarvae2_cesm126[, x])/30)
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
cesm_temps3 <- readRDS(here('data', 'cesm_forecast_temp3.rds'))
cesm_salts3 <- readRDS(here('data', 'cesm_forecast_salt3.rds'))

preds_fhslarvae3_cesm126 <- pred_loop(2070:2099, fhs_larvae, 130,
                                      5, 'ssp126', cesm_temps3, 
                                      cesm_salts3, larval_formula)

# Combine into one data frame
df_fhslarvae3_cesm126 <- list(preds_fhslarvae3_cesm126[[1]], preds_fhslarvae3_cesm126[[2]],
                              preds_fhslarvae3_cesm126[[3]], preds_fhslarvae3_cesm126[[4]],
                              preds_fhslarvae3_cesm126[[5]], preds_fhslarvae3_cesm126[[6]],
                              preds_fhslarvae3_cesm126[[7]], preds_fhslarvae3_cesm126[[8]],
                              preds_fhslarvae3_cesm126[[9]], preds_fhslarvae3_cesm126[[10]],
                              preds_fhslarvae3_cesm126[[11]], preds_fhslarvae3_cesm126[[12]],
                              preds_fhslarvae3_cesm126[[13]], preds_fhslarvae3_cesm126[[14]],
                              preds_fhslarvae3_cesm126[[15]], preds_fhslarvae3_cesm126[[16]],
                              preds_fhslarvae3_cesm126[[17]], preds_fhslarvae3_cesm126[[18]],
                              preds_fhslarvae3_cesm126[[19]], preds_fhslarvae3_cesm126[[20]],
                              preds_fhslarvae3_cesm126[[21]], preds_fhslarvae3_cesm126[[22]],
                              preds_fhslarvae3_cesm126[[23]], preds_fhslarvae3_cesm126[[24]],
                              preds_fhslarvae3_cesm126[[25]], preds_fhslarvae3_cesm126[[26]], 
                              preds_fhslarvae3_cesm126[[27]], preds_fhslarvae3_cesm126[[28]], 
                              preds_fhslarvae3_cesm126[[29]], preds_fhslarvae3_cesm126[[30]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhslarvae3_cesm126), fixed = T)
df_fhslarvae_avg3_cesm126 <- data.frame(lat = df_fhslarvae3_cesm126$lat, 
                                        lon = df_fhslarvae3_cesm126$lon, 
                                        dist = df_fhslarvae3_cesm126$dist,
                                        avg_pred = rowSums(df_fhslarvae3_cesm126[, x])/30)
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


rm(cesm_temps1, cesm_temps2, cesm_temps3,
   cesm_salts1, cesm_salts2, cesm_salts3)

##### CESM 585 -----------------------------------------------------------------------------------------------------------------
## 2015 - 2039
cesm_temps1 <- readRDS(here('data', 'cesm_forecast_temp1.rds'))
cesm_salts1 <- readRDS(here('data', 'cesm_forecast_salt1.rds'))

preds_fhslarvae1_cesm585 <- pred_loop(2015:2039, fhs_larvae, 130,
                                      5, 'ssp585', cesm_temps1, 
                                      cesm_salts1, larval_formula)

# Combine into one data frame
df_fhslarvae1_cesm585 <- list(preds_fhslarvae1_cesm585[[1]], preds_fhslarvae1_cesm585[[2]],
                              preds_fhslarvae1_cesm585[[3]], preds_fhslarvae1_cesm585[[4]],
                              preds_fhslarvae1_cesm585[[5]], preds_fhslarvae1_cesm585[[6]],
                              preds_fhslarvae1_cesm585[[7]], preds_fhslarvae1_cesm585[[8]],
                              preds_fhslarvae1_cesm585[[9]], preds_fhslarvae1_cesm585[[10]],
                              preds_fhslarvae1_cesm585[[11]], preds_fhslarvae1_cesm585[[12]],
                              preds_fhslarvae1_cesm585[[13]], preds_fhslarvae1_cesm585[[14]],
                              preds_fhslarvae1_cesm585[[15]], preds_fhslarvae1_cesm585[[16]],
                              preds_fhslarvae1_cesm585[[17]], preds_fhslarvae1_cesm585[[18]],
                              preds_fhslarvae1_cesm585[[19]], preds_fhslarvae1_cesm585[[20]],
                              preds_fhslarvae1_cesm585[[21]], preds_fhslarvae1_cesm585[[22]],
                              preds_fhslarvae1_cesm585[[23]], preds_fhslarvae1_cesm585[[24]],
                              preds_fhslarvae1_cesm585[[25]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhslarvae1_cesm585), fixed = T)
df_fhslarvae_avg1_cesm585 <- data.frame(lat = df_fhslarvae1_cesm585$lat, 
                                        lon = df_fhslarvae1_cesm585$lon, 
                                        dist = df_fhslarvae1_cesm585$dist,
                                        avg_pred = rowSums(df_fhslarvae1_cesm585[, x])/25)
saveRDS(df_fhslarvae_avg1_cesm585, file = here("data", "df_fhslarvae_avg1_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae_avg1_cesm585, "Forecasted Distribution 2015 - 2039 \n CESM SSP585")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_cesm_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
cesm_temps2 <- readRDS(here('data', 'cesm_forecast_temp2.rds'))
cesm_salts2 <- readRDS(here('data', 'cesm_forecast_salt2.rds'))

preds_fhslarvae2_cesm585 <- pred_loop(2040:2069, fhs_larvae, 130,
                                      5, 'ssp585', cesm_temps2, 
                                      cesm_salts2, larval_formula)

# Combine into one data frame
df_fhslarvae2_cesm585 <- list(preds_fhslarvae2_cesm585[[1]], preds_fhslarvae2_cesm585[[2]],
                              preds_fhslarvae2_cesm585[[3]], preds_fhslarvae2_cesm585[[4]],
                              preds_fhslarvae2_cesm585[[5]], preds_fhslarvae2_cesm585[[6]],
                              preds_fhslarvae2_cesm585[[7]], preds_fhslarvae2_cesm585[[8]],
                              preds_fhslarvae2_cesm585[[9]], preds_fhslarvae2_cesm585[[10]],
                              preds_fhslarvae2_cesm585[[11]], preds_fhslarvae2_cesm585[[12]],
                              preds_fhslarvae2_cesm585[[13]], preds_fhslarvae2_cesm585[[14]],
                              preds_fhslarvae2_cesm585[[15]], preds_fhslarvae2_cesm585[[16]],
                              preds_fhslarvae2_cesm585[[17]], preds_fhslarvae2_cesm585[[18]],
                              preds_fhslarvae2_cesm585[[19]], preds_fhslarvae2_cesm585[[20]],
                              preds_fhslarvae2_cesm585[[21]], preds_fhslarvae2_cesm585[[22]],
                              preds_fhslarvae2_cesm585[[23]], preds_fhslarvae2_cesm585[[24]],
                              preds_fhslarvae2_cesm585[[25]], preds_fhslarvae2_cesm585[[26]],
                              preds_fhslarvae2_cesm585[[27]], preds_fhslarvae2_cesm585[[28]],
                              preds_fhslarvae2_cesm585[[29]], preds_fhslarvae2_cesm585[[30]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhslarvae2_cesm585), fixed = T)
df_fhslarvae_avg2_cesm585 <- data.frame(lat = df_fhslarvae2_cesm585$lat, 
                                        lon = df_fhslarvae2_cesm585$lon, 
                                        dist = df_fhslarvae2_cesm585$dist,
                                        avg_pred = rowSums(df_fhslarvae2_cesm585[, x])/30)
saveRDS(df_fhslarvae_avg2_cesm585, file = here("data", "df_fhslarvae_avg2_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae_avg2_cesm585, "Forecasted Distribution 2040 - 2069 \n CESM SSP585")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_cesm_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
cesm_temps3 <- readRDS(here('data', 'cesm_forecast_temp3.rds'))
cesm_salts3 <- readRDS(here('data', 'cesm_forecast_salt3.rds'))

preds_fhslarvae3_cesm585 <- pred_loop(2070:2099, fhs_larvae, 130,
                                      5, 'ssp585', cesm_temps3, 
                                      cesm_salts3, larval_formula)

# Combine into one data frame
df_fhslarvae3_cesm585 <- list(preds_fhslarvae3_cesm585[[1]], preds_fhslarvae3_cesm585[[2]],
                              preds_fhslarvae3_cesm585[[3]], preds_fhslarvae3_cesm585[[4]],
                              preds_fhslarvae3_cesm585[[5]], preds_fhslarvae3_cesm585[[6]],
                              preds_fhslarvae3_cesm585[[7]], preds_fhslarvae3_cesm585[[8]],
                              preds_fhslarvae3_cesm585[[9]], preds_fhslarvae3_cesm585[[10]],
                              preds_fhslarvae3_cesm585[[11]], preds_fhslarvae3_cesm585[[12]],
                              preds_fhslarvae3_cesm585[[13]], preds_fhslarvae3_cesm585[[14]],
                              preds_fhslarvae3_cesm585[[15]], preds_fhslarvae3_cesm585[[16]],
                              preds_fhslarvae3_cesm585[[17]], preds_fhslarvae3_cesm585[[18]],
                              preds_fhslarvae3_cesm585[[19]], preds_fhslarvae3_cesm585[[20]],
                              preds_fhslarvae3_cesm585[[21]], preds_fhslarvae3_cesm585[[22]],
                              preds_fhslarvae3_cesm585[[23]], preds_fhslarvae3_cesm585[[24]],
                              preds_fhslarvae3_cesm585[[25]], preds_fhslarvae3_cesm585[[26]], 
                              preds_fhslarvae3_cesm585[[27]], preds_fhslarvae3_cesm585[[28]], 
                              preds_fhslarvae3_cesm585[[29]], preds_fhslarvae3_cesm585[[30]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhslarvae3_cesm585), fixed = T)
df_fhslarvae_avg3_cesm585 <- data.frame(lat = df_fhslarvae3_cesm585$lat, 
                                        lon = df_fhslarvae3_cesm585$lon, 
                                        dist = df_fhslarvae3_cesm585$dist,
                                        avg_pred = rowSums(df_fhslarvae3_cesm585[, x])/30)
saveRDS(df_fhslarvae_avg3_cesm585, file = here("data", "df_fhslarvae_avg3_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae_avg3_cesm585, "Forecasted Distribution 2070 - 2099 \n CESM SSP585")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_cesm_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


rm(cesm_temps1, cesm_temps2, cesm_temps3,
   cesm_salts1, cesm_salts2, cesm_salts3)

##### GFDL 126 ----------------------------------------------------------------------------------------------------------------------------
gfdl_temps1 <- readRDS(here('data', 'gfdl_forecast_temp1.rds'))
gfdl_salts1 <- readRDS(here('data', 'gfdl_forecast_salt1.rds'))

preds_fhslarvae1_gfdl126 <- pred_loop(2015:2039, fhs_larvae, 130,
                                      5, 'ssp126', gfdl_temps1, 
                                      gfdl_salts1, larval_formula)

# Combine into one data frame
df_fhslarvae1_gfdl126 <- list(preds_fhslarvae1_gfdl126[[1]], preds_fhslarvae1_gfdl126[[2]],
                              preds_fhslarvae1_gfdl126[[3]], preds_fhslarvae1_gfdl126[[4]],
                              preds_fhslarvae1_gfdl126[[5]], preds_fhslarvae1_gfdl126[[6]],
                              preds_fhslarvae1_gfdl126[[7]], preds_fhslarvae1_gfdl126[[8]],
                              preds_fhslarvae1_gfdl126[[9]], preds_fhslarvae1_gfdl126[[10]],
                              preds_fhslarvae1_gfdl126[[11]], preds_fhslarvae1_gfdl126[[12]],
                              preds_fhslarvae1_gfdl126[[13]], preds_fhslarvae1_gfdl126[[14]],
                              preds_fhslarvae1_gfdl126[[15]], preds_fhslarvae1_gfdl126[[16]],
                              preds_fhslarvae1_gfdl126[[17]], preds_fhslarvae1_gfdl126[[18]],
                              preds_fhslarvae1_gfdl126[[19]], preds_fhslarvae1_gfdl126[[20]],
                              preds_fhslarvae1_gfdl126[[21]], preds_fhslarvae1_gfdl126[[22]],
                              preds_fhslarvae1_gfdl126[[23]], preds_fhslarvae1_gfdl126[[24]],
                              preds_fhslarvae1_gfdl126[[25]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhslarvae1_gfdl126), fixed = T)
df_fhslarvae_avg1_gfdl126 <- data.frame(lat = df_fhslarvae1_gfdl126$lat, 
                                        lon = df_fhslarvae1_gfdl126$lon, 
                                        dist = df_fhslarvae1_gfdl126$dist,
                                        avg_pred = rowSums(df_fhslarvae1_gfdl126[, x])/25)
saveRDS(df_fhslarvae_avg1_gfdl126, file = here("data", "df_fhslarvae_avg1_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae_avg1_gfdl126, "Forecasted Distribution 2015 - 2039 \n gfdl SSP126")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_gfdl_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
gfdl_temps2 <- readRDS(here('data', 'gfdl_forecast_temp2.rds'))
gfdl_salts2 <- readRDS(here('data', 'gfdl_forecast_salt2.rds'))

preds_fhslarvae2_gfdl126 <- pred_loop(2040:2069, fhs_larvae, 130,
                                      5, 'ssp126', gfdl_temps2, 
                                      gfdl_salts2, larval_formula)

# Combine into one data frame
df_fhslarvae2_gfdl126 <- list(preds_fhslarvae2_gfdl126[[1]], preds_fhslarvae2_gfdl126[[2]],
                              preds_fhslarvae2_gfdl126[[3]], preds_fhslarvae2_gfdl126[[4]],
                              preds_fhslarvae2_gfdl126[[5]], preds_fhslarvae2_gfdl126[[6]],
                              preds_fhslarvae2_gfdl126[[7]], preds_fhslarvae2_gfdl126[[8]],
                              preds_fhslarvae2_gfdl126[[9]], preds_fhslarvae2_gfdl126[[10]],
                              preds_fhslarvae2_gfdl126[[11]], preds_fhslarvae2_gfdl126[[12]],
                              preds_fhslarvae2_gfdl126[[13]], preds_fhslarvae2_gfdl126[[14]],
                              preds_fhslarvae2_gfdl126[[15]], preds_fhslarvae2_gfdl126[[16]],
                              preds_fhslarvae2_gfdl126[[17]], preds_fhslarvae2_gfdl126[[18]],
                              preds_fhslarvae2_gfdl126[[19]], preds_fhslarvae2_gfdl126[[20]],
                              preds_fhslarvae2_gfdl126[[21]], preds_fhslarvae2_gfdl126[[22]],
                              preds_fhslarvae2_gfdl126[[23]], preds_fhslarvae2_gfdl126[[24]],
                              preds_fhslarvae2_gfdl126[[25]], preds_fhslarvae2_gfdl126[[26]],
                              preds_fhslarvae2_gfdl126[[27]], preds_fhslarvae2_gfdl126[[28]],
                              preds_fhslarvae2_gfdl126[[29]], preds_fhslarvae2_gfdl126[[30]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhslarvae2_gfdl126), fixed = T)
df_fhslarvae_avg2_gfdl126 <- data.frame(lat = df_fhslarvae2_gfdl126$lat, 
                                        lon = df_fhslarvae2_gfdl126$lon, 
                                        dist = df_fhslarvae2_gfdl126$dist,
                                        avg_pred = rowSums(df_fhslarvae2_gfdl126[, x])/30)
saveRDS(df_fhslarvae_avg2_gfdl126, file = here("data", "df_fhslarvae_avg2_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae_avg2_gfdl126, "Forecasted Distribution 2040 - 2069 \n gfdl SSP126")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_gfdl_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
gfdl_temps3 <- readRDS(here('data', 'gfdl_forecast_temp3.rds'))
gfdl_salts3 <- readRDS(here('data', 'gfdl_forecast_salt3.rds'))

preds_fhslarvae3_gfdl126 <- pred_loop(2070:2099, fhs_larvae, 130,
                                      5, 'ssp126', gfdl_temps3, 
                                      gfdl_salts3, larval_formula)

# Combine into one data frame
df_fhslarvae3_gfdl126 <- list(preds_fhslarvae3_gfdl126[[1]], preds_fhslarvae3_gfdl126[[2]],
                              preds_fhslarvae3_gfdl126[[3]], preds_fhslarvae3_gfdl126[[4]],
                              preds_fhslarvae3_gfdl126[[5]], preds_fhslarvae3_gfdl126[[6]],
                              preds_fhslarvae3_gfdl126[[7]], preds_fhslarvae3_gfdl126[[8]],
                              preds_fhslarvae3_gfdl126[[9]], preds_fhslarvae3_gfdl126[[10]],
                              preds_fhslarvae3_gfdl126[[11]], preds_fhslarvae3_gfdl126[[12]],
                              preds_fhslarvae3_gfdl126[[13]], preds_fhslarvae3_gfdl126[[14]],
                              preds_fhslarvae3_gfdl126[[15]], preds_fhslarvae3_gfdl126[[16]],
                              preds_fhslarvae3_gfdl126[[17]], preds_fhslarvae3_gfdl126[[18]],
                              preds_fhslarvae3_gfdl126[[19]], preds_fhslarvae3_gfdl126[[20]],
                              preds_fhslarvae3_gfdl126[[21]], preds_fhslarvae3_gfdl126[[22]],
                              preds_fhslarvae3_gfdl126[[23]], preds_fhslarvae3_gfdl126[[24]],
                              preds_fhslarvae3_gfdl126[[25]], preds_fhslarvae3_gfdl126[[26]], 
                              preds_fhslarvae3_gfdl126[[27]], preds_fhslarvae3_gfdl126[[28]], 
                              preds_fhslarvae3_gfdl126[[29]], preds_fhslarvae3_gfdl126[[30]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhslarvae3_gfdl126), fixed = T)
df_fhslarvae_avg3_gfdl126 <- data.frame(lat = df_fhslarvae3_gfdl126$lat, 
                                        lon = df_fhslarvae3_gfdl126$lon, 
                                        dist = df_fhslarvae3_gfdl126$dist,
                                        avg_pred = rowSums(df_fhslarvae3_gfdl126[, x])/30)
saveRDS(df_fhslarvae_avg3_gfdl126, file = here("data", "df_fhslarvae_avg3_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae_avg3_gfdl126, "Forecasted Distribution 2070 - 2099 \n gfdl SSP126")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_gfdl_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


rm(gfdl_temps1, gfdl_temps2, gfdl_temps3,
   gfdl_salts1, gfdl_salts2, gfdl_salts3)

##### GFDL 585 -----------------------------------------------------------------------------------------------------------------
## 2015 - 2039
gfdl_temps1 <- readRDS(here('data', 'gfdl_forecast_temp1.rds'))
gfdl_salts1 <- readRDS(here('data', 'gfdl_forecast_salt1.rds'))

preds_fhslarvae1_gfdl585 <- pred_loop(2015:2039, fhs_larvae, 130,
                                      5, 'ssp585', gfdl_temps1, 
                                      gfdl_salts1, larval_formula)

# Combine into one data frame
df_fhslarvae1_gfdl585 <- list(preds_fhslarvae1_gfdl585[[1]], preds_fhslarvae1_gfdl585[[2]],
                              preds_fhslarvae1_gfdl585[[3]], preds_fhslarvae1_gfdl585[[4]],
                              preds_fhslarvae1_gfdl585[[5]], preds_fhslarvae1_gfdl585[[6]],
                              preds_fhslarvae1_gfdl585[[7]], preds_fhslarvae1_gfdl585[[8]],
                              preds_fhslarvae1_gfdl585[[9]], preds_fhslarvae1_gfdl585[[10]],
                              preds_fhslarvae1_gfdl585[[11]], preds_fhslarvae1_gfdl585[[12]],
                              preds_fhslarvae1_gfdl585[[13]], preds_fhslarvae1_gfdl585[[14]],
                              preds_fhslarvae1_gfdl585[[15]], preds_fhslarvae1_gfdl585[[16]],
                              preds_fhslarvae1_gfdl585[[17]], preds_fhslarvae1_gfdl585[[18]],
                              preds_fhslarvae1_gfdl585[[19]], preds_fhslarvae1_gfdl585[[20]],
                              preds_fhslarvae1_gfdl585[[21]], preds_fhslarvae1_gfdl585[[22]],
                              preds_fhslarvae1_gfdl585[[23]], preds_fhslarvae1_gfdl585[[24]],
                              preds_fhslarvae1_gfdl585[[25]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhslarvae1_gfdl585), fixed = T)
df_fhslarvae_avg1_gfdl585 <- data.frame(lat = df_fhslarvae1_gfdl585$lat, 
                                        lon = df_fhslarvae1_gfdl585$lon, 
                                        dist = df_fhslarvae1_gfdl585$dist,
                                        avg_pred = rowSums(df_fhslarvae1_gfdl585[, x])/25)
saveRDS(df_fhslarvae_avg1_gfdl585, file = here("data", "df_fhslarvae_avg1_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae_avg1_gfdl585, "Forecasted Distribution 2015 - 2039 \n gfdl SSP585")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_gfdl_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
gfdl_temps2 <- readRDS(here('data', 'gfdl_forecast_temp2.rds'))
gfdl_salts2 <- readRDS(here('data', 'gfdl_forecast_salt2.rds'))

preds_fhslarvae2_gfdl585 <- pred_loop(2040:2069, fhs_larvae, 130,
                                      5, 'ssp585', gfdl_temps2, 
                                      gfdl_salts2, larval_formula)

# Combine into one data frame
df_fhslarvae2_gfdl585 <- list(preds_fhslarvae2_gfdl585[[1]], preds_fhslarvae2_gfdl585[[2]],
                              preds_fhslarvae2_gfdl585[[3]], preds_fhslarvae2_gfdl585[[4]],
                              preds_fhslarvae2_gfdl585[[5]], preds_fhslarvae2_gfdl585[[6]],
                              preds_fhslarvae2_gfdl585[[7]], preds_fhslarvae2_gfdl585[[8]],
                              preds_fhslarvae2_gfdl585[[9]], preds_fhslarvae2_gfdl585[[10]],
                              preds_fhslarvae2_gfdl585[[11]], preds_fhslarvae2_gfdl585[[12]],
                              preds_fhslarvae2_gfdl585[[13]], preds_fhslarvae2_gfdl585[[14]],
                              preds_fhslarvae2_gfdl585[[15]], preds_fhslarvae2_gfdl585[[16]],
                              preds_fhslarvae2_gfdl585[[17]], preds_fhslarvae2_gfdl585[[18]],
                              preds_fhslarvae2_gfdl585[[19]], preds_fhslarvae2_gfdl585[[20]],
                              preds_fhslarvae2_gfdl585[[21]], preds_fhslarvae2_gfdl585[[22]],
                              preds_fhslarvae2_gfdl585[[23]], preds_fhslarvae2_gfdl585[[24]],
                              preds_fhslarvae2_gfdl585[[25]], preds_fhslarvae2_gfdl585[[26]],
                              preds_fhslarvae2_gfdl585[[27]], preds_fhslarvae2_gfdl585[[28]],
                              preds_fhslarvae2_gfdl585[[29]], preds_fhslarvae2_gfdl585[[30]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhslarvae2_gfdl585), fixed = T)
df_fhslarvae_avg2_gfdl585 <- data.frame(lat = df_fhslarvae2_gfdl585$lat, 
                                        lon = df_fhslarvae2_gfdl585$lon, 
                                        dist = df_fhslarvae2_gfdl585$dist,
                                        avg_pred = rowSums(df_fhslarvae2_gfdl585[, x])/30)
saveRDS(df_fhslarvae_avg2_gfdl585, file = here("data", "df_fhslarvae_avg2_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae_avg2_gfdl585, "Forecasted Distribution 2040 - 2069 \n gfdl SSP585")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_gfdl_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
gfdl_temps3 <- readRDS(here('data', 'gfdl_forecast_temp3.rds'))
gfdl_salts3 <- readRDS(here('data', 'gfdl_forecast_salt3.rds'))

preds_fhslarvae3_gfdl585 <- pred_loop(2070:2099, fhs_larvae, 130,
                                      5, 'ssp585', gfdl_temps3, 
                                      gfdl_salts3, larval_formula)

# Combine into one data frame
df_fhslarvae3_gfdl585 <- list(preds_fhslarvae3_gfdl585[[1]], preds_fhslarvae3_gfdl585[[2]],
                              preds_fhslarvae3_gfdl585[[3]], preds_fhslarvae3_gfdl585[[4]],
                              preds_fhslarvae3_gfdl585[[5]], preds_fhslarvae3_gfdl585[[6]],
                              preds_fhslarvae3_gfdl585[[7]], preds_fhslarvae3_gfdl585[[8]],
                              preds_fhslarvae3_gfdl585[[9]], preds_fhslarvae3_gfdl585[[10]],
                              preds_fhslarvae3_gfdl585[[11]], preds_fhslarvae3_gfdl585[[12]],
                              preds_fhslarvae3_gfdl585[[13]], preds_fhslarvae3_gfdl585[[14]],
                              preds_fhslarvae3_gfdl585[[15]], preds_fhslarvae3_gfdl585[[16]],
                              preds_fhslarvae3_gfdl585[[17]], preds_fhslarvae3_gfdl585[[18]],
                              preds_fhslarvae3_gfdl585[[19]], preds_fhslarvae3_gfdl585[[20]],
                              preds_fhslarvae3_gfdl585[[21]], preds_fhslarvae3_gfdl585[[22]],
                              preds_fhslarvae3_gfdl585[[23]], preds_fhslarvae3_gfdl585[[24]],
                              preds_fhslarvae3_gfdl585[[25]], preds_fhslarvae3_gfdl585[[26]], 
                              preds_fhslarvae3_gfdl585[[27]], preds_fhslarvae3_gfdl585[[28]], 
                              preds_fhslarvae3_gfdl585[[29]], preds_fhslarvae3_gfdl585[[30]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhslarvae3_gfdl585), fixed = T)
df_fhslarvae_avg3_gfdl585 <- data.frame(lat = df_fhslarvae3_gfdl585$lat, 
                                        lon = df_fhslarvae3_gfdl585$lon, 
                                        dist = df_fhslarvae3_gfdl585$dist,
                                        avg_pred = rowSums(df_fhslarvae3_gfdl585[, x])/30)
saveRDS(df_fhslarvae_avg3_gfdl585, file = here("data", "df_fhslarvae_avg3_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae_avg3_gfdl585, "Forecasted Distribution 2070 - 2099 \n gfdl SSP585")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_gfdl_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


rm(gfdl_temps1, gfdl_temps2, gfdl_temps3,
   gfdl_salts1, gfdl_salts2, gfdl_salts3)


##### MIROC 126 ----------------------------------------------------------------------------------------------------------------------------
miroc_temps1 <- readRDS(here('data', 'miroc_forecast_temp1.rds'))
miroc_salts1 <- readRDS(here('data', 'miroc_forecast_salt1.rds'))

preds_fhslarvae1_miroc126 <- pred_loop(2015:2039, fhs_larvae, 130,
                                       5, 'ssp126', miroc_temps1, 
                                       miroc_salts1, larval_formula)

# Combine into one data frame
df_fhslarvae1_miroc126 <- list(preds_fhslarvae1_miroc126[[1]], preds_fhslarvae1_miroc126[[2]],
                               preds_fhslarvae1_miroc126[[3]], preds_fhslarvae1_miroc126[[4]],
                               preds_fhslarvae1_miroc126[[5]], preds_fhslarvae1_miroc126[[6]],
                               preds_fhslarvae1_miroc126[[7]], preds_fhslarvae1_miroc126[[8]],
                               preds_fhslarvae1_miroc126[[9]], preds_fhslarvae1_miroc126[[10]],
                               preds_fhslarvae1_miroc126[[11]], preds_fhslarvae1_miroc126[[12]],
                               preds_fhslarvae1_miroc126[[13]], preds_fhslarvae1_miroc126[[14]],
                               preds_fhslarvae1_miroc126[[15]], preds_fhslarvae1_miroc126[[16]],
                               preds_fhslarvae1_miroc126[[17]], preds_fhslarvae1_miroc126[[18]],
                               preds_fhslarvae1_miroc126[[19]], preds_fhslarvae1_miroc126[[20]],
                               preds_fhslarvae1_miroc126[[21]], preds_fhslarvae1_miroc126[[22]],
                               preds_fhslarvae1_miroc126[[23]], preds_fhslarvae1_miroc126[[24]],
                               preds_fhslarvae1_miroc126[[25]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhslarvae1_miroc126), fixed = T)
df_fhslarvae_avg1_miroc126 <- data.frame(lat = df_fhslarvae1_miroc126$lat, 
                                         lon = df_fhslarvae1_miroc126$lon, 
                                         dist = df_fhslarvae1_miroc126$dist,
                                         avg_pred = rowSums(df_fhslarvae1_miroc126[, x])/25)
saveRDS(df_fhslarvae_avg1_miroc126, file = here("data", "df_fhslarvae_avg1_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae_avg1_miroc126, "Forecasted Distribution 2015 - 2039 \n miroc SSP126")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_miroc_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
miroc_temps2 <- readRDS(here('data', 'miroc_forecast_temp2.rds'))
miroc_salts2 <- readRDS(here('data', 'miroc_forecast_salt2.rds'))

preds_fhslarvae2_miroc126 <- pred_loop(2040:2069, fhs_larvae, 130,
                                       5, 'ssp126', miroc_temps2, 
                                       miroc_salts2, larval_formula)

# Combine into one data frame
df_fhslarvae2_miroc126 <- list(preds_fhslarvae2_miroc126[[1]], preds_fhslarvae2_miroc126[[2]],
                               preds_fhslarvae2_miroc126[[3]], preds_fhslarvae2_miroc126[[4]],
                               preds_fhslarvae2_miroc126[[5]], preds_fhslarvae2_miroc126[[6]],
                               preds_fhslarvae2_miroc126[[7]], preds_fhslarvae2_miroc126[[8]],
                               preds_fhslarvae2_miroc126[[9]], preds_fhslarvae2_miroc126[[10]],
                               preds_fhslarvae2_miroc126[[11]], preds_fhslarvae2_miroc126[[12]],
                               preds_fhslarvae2_miroc126[[13]], preds_fhslarvae2_miroc126[[14]],
                               preds_fhslarvae2_miroc126[[15]], preds_fhslarvae2_miroc126[[16]],
                               preds_fhslarvae2_miroc126[[17]], preds_fhslarvae2_miroc126[[18]],
                               preds_fhslarvae2_miroc126[[19]], preds_fhslarvae2_miroc126[[20]],
                               preds_fhslarvae2_miroc126[[21]], preds_fhslarvae2_miroc126[[22]],
                               preds_fhslarvae2_miroc126[[23]], preds_fhslarvae2_miroc126[[24]],
                               preds_fhslarvae2_miroc126[[25]], preds_fhslarvae2_miroc126[[26]],
                               preds_fhslarvae2_miroc126[[27]], preds_fhslarvae2_miroc126[[28]],
                               preds_fhslarvae2_miroc126[[29]], preds_fhslarvae2_miroc126[[30]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhslarvae2_miroc126), fixed = T)
df_fhslarvae_avg2_miroc126 <- data.frame(lat = df_fhslarvae2_miroc126$lat, 
                                         lon = df_fhslarvae2_miroc126$lon, 
                                         dist = df_fhslarvae2_miroc126$dist,
                                         avg_pred = rowSums(df_fhslarvae2_miroc126[, x])/30)
saveRDS(df_fhslarvae_avg2_miroc126, file = here("data", "df_fhslarvae_avg2_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae_avg2_miroc126, "Forecasted Distribution 2040 - 2069 \n miroc SSP126")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_miroc_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
miroc_temps3 <- readRDS(here('data', 'miroc_forecast_temp3.rds'))
miroc_salts3 <- readRDS(here('data', 'miroc_forecast_salt3.rds'))

preds_fhslarvae3_miroc126 <- pred_loop(2070:2099, fhs_larvae, 130,
                                       5, 'ssp126', miroc_temps3, 
                                       miroc_salts3, larval_formula)

# Combine into one data frame
df_fhslarvae3_miroc126 <- list(preds_fhslarvae3_miroc126[[1]], preds_fhslarvae3_miroc126[[2]],
                               preds_fhslarvae3_miroc126[[3]], preds_fhslarvae3_miroc126[[4]],
                               preds_fhslarvae3_miroc126[[5]], preds_fhslarvae3_miroc126[[6]],
                               preds_fhslarvae3_miroc126[[7]], preds_fhslarvae3_miroc126[[8]],
                               preds_fhslarvae3_miroc126[[9]], preds_fhslarvae3_miroc126[[10]],
                               preds_fhslarvae3_miroc126[[11]], preds_fhslarvae3_miroc126[[12]],
                               preds_fhslarvae3_miroc126[[13]], preds_fhslarvae3_miroc126[[14]],
                               preds_fhslarvae3_miroc126[[15]], preds_fhslarvae3_miroc126[[16]],
                               preds_fhslarvae3_miroc126[[17]], preds_fhslarvae3_miroc126[[18]],
                               preds_fhslarvae3_miroc126[[19]], preds_fhslarvae3_miroc126[[20]],
                               preds_fhslarvae3_miroc126[[21]], preds_fhslarvae3_miroc126[[22]],
                               preds_fhslarvae3_miroc126[[23]], preds_fhslarvae3_miroc126[[24]],
                               preds_fhslarvae3_miroc126[[25]], preds_fhslarvae3_miroc126[[26]], 
                               preds_fhslarvae3_miroc126[[27]], preds_fhslarvae3_miroc126[[28]], 
                               preds_fhslarvae3_miroc126[[29]], preds_fhslarvae3_miroc126[[30]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhslarvae3_miroc126), fixed = T)
df_fhslarvae_avg3_miroc126 <- data.frame(lat = df_fhslarvae3_miroc126$lat, 
                                         lon = df_fhslarvae3_miroc126$lon, 
                                         dist = df_fhslarvae3_miroc126$dist,
                                         avg_pred = rowSums(df_fhslarvae3_miroc126[, x])/30)
saveRDS(df_fhslarvae_avg3_miroc126, file = here("data", "df_fhslarvae_avg3_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae_avg3_miroc126, "Forecasted Distribution 2070 - 2099 \n miroc SSP126")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_miroc_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


rm(miroc_temps1, miroc_temps2, miroc_temps3,
   miroc_salts1, miroc_salts2, miroc_salts3)

##### MIROC 585 -----------------------------------------------------------------------------------------------------------------
## 2015 - 2039
miroc_temps1 <- readRDS(here('data', 'miroc_forecast_temp1.rds'))
miroc_salts1 <- readRDS(here('data', 'miroc_forecast_salt1.rds'))

preds_fhslarvae1_miroc585 <- pred_loop(2015:2039, fhs_larvae, 130,
                                       5, 'ssp585', miroc_temps1, 
                                       miroc_salts1, larval_formula)

# Combine into one data frame
df_fhslarvae1_miroc585 <- list(preds_fhslarvae1_miroc585[[1]], preds_fhslarvae1_miroc585[[2]],
                               preds_fhslarvae1_miroc585[[3]], preds_fhslarvae1_miroc585[[4]],
                               preds_fhslarvae1_miroc585[[5]], preds_fhslarvae1_miroc585[[6]],
                               preds_fhslarvae1_miroc585[[7]], preds_fhslarvae1_miroc585[[8]],
                               preds_fhslarvae1_miroc585[[9]], preds_fhslarvae1_miroc585[[10]],
                               preds_fhslarvae1_miroc585[[11]], preds_fhslarvae1_miroc585[[12]],
                               preds_fhslarvae1_miroc585[[13]], preds_fhslarvae1_miroc585[[14]],
                               preds_fhslarvae1_miroc585[[15]], preds_fhslarvae1_miroc585[[16]],
                               preds_fhslarvae1_miroc585[[17]], preds_fhslarvae1_miroc585[[18]],
                               preds_fhslarvae1_miroc585[[19]], preds_fhslarvae1_miroc585[[20]],
                               preds_fhslarvae1_miroc585[[21]], preds_fhslarvae1_miroc585[[22]],
                               preds_fhslarvae1_miroc585[[23]], preds_fhslarvae1_miroc585[[24]],
                               preds_fhslarvae1_miroc585[[25]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhslarvae1_miroc585), fixed = T)
df_fhslarvae_avg1_miroc585 <- data.frame(lat = df_fhslarvae1_miroc585$lat, 
                                         lon = df_fhslarvae1_miroc585$lon, 
                                         dist = df_fhslarvae1_miroc585$dist,
                                         avg_pred = rowSums(df_fhslarvae1_miroc585[, x])/25)
saveRDS(df_fhslarvae_avg1_miroc585, file = here("data", "df_fhslarvae_avg1_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae_avg1_miroc585, "Forecasted Distribution 2015 - 2039 \n miroc SSP585")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_miroc_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2040 - 2069
miroc_temps2 <- readRDS(here('data', 'miroc_forecast_temp2.rds'))
miroc_salts2 <- readRDS(here('data', 'miroc_forecast_salt2.rds'))

preds_fhslarvae2_miroc585 <- pred_loop(2040:2069, fhs_larvae, 130,
                                       5, 'ssp585', miroc_temps2, 
                                       miroc_salts2, larval_formula)

# Combine into one data frame
df_fhslarvae2_miroc585 <- list(preds_fhslarvae2_miroc585[[1]], preds_fhslarvae2_miroc585[[2]],
                               preds_fhslarvae2_miroc585[[3]], preds_fhslarvae2_miroc585[[4]],
                               preds_fhslarvae2_miroc585[[5]], preds_fhslarvae2_miroc585[[6]],
                               preds_fhslarvae2_miroc585[[7]], preds_fhslarvae2_miroc585[[8]],
                               preds_fhslarvae2_miroc585[[9]], preds_fhslarvae2_miroc585[[10]],
                               preds_fhslarvae2_miroc585[[11]], preds_fhslarvae2_miroc585[[12]],
                               preds_fhslarvae2_miroc585[[13]], preds_fhslarvae2_miroc585[[14]],
                               preds_fhslarvae2_miroc585[[15]], preds_fhslarvae2_miroc585[[16]],
                               preds_fhslarvae2_miroc585[[17]], preds_fhslarvae2_miroc585[[18]],
                               preds_fhslarvae2_miroc585[[19]], preds_fhslarvae2_miroc585[[20]],
                               preds_fhslarvae2_miroc585[[21]], preds_fhslarvae2_miroc585[[22]],
                               preds_fhslarvae2_miroc585[[23]], preds_fhslarvae2_miroc585[[24]],
                               preds_fhslarvae2_miroc585[[25]], preds_fhslarvae2_miroc585[[26]],
                               preds_fhslarvae2_miroc585[[27]], preds_fhslarvae2_miroc585[[28]],
                               preds_fhslarvae2_miroc585[[29]], preds_fhslarvae2_miroc585[[30]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhslarvae2_miroc585), fixed = T)
df_fhslarvae_avg2_miroc585 <- data.frame(lat = df_fhslarvae2_miroc585$lat, 
                                         lon = df_fhslarvae2_miroc585$lon, 
                                         dist = df_fhslarvae2_miroc585$dist,
                                         avg_pred = rowSums(df_fhslarvae2_miroc585[, x])/30)
saveRDS(df_fhslarvae_avg2_miroc585, file = here("data", "df_fhslarvae_avg2_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae_avg2_miroc585, "Forecasted Distribution 2040 - 2069 \n miroc SSP585")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_miroc_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


## 2070 - 2099
miroc_temps3 <- readRDS(here('data', 'miroc_forecast_temp3.rds'))
miroc_salts3 <- readRDS(here('data', 'miroc_forecast_salt3.rds'))

preds_fhslarvae3_miroc585 <- pred_loop(2070:2099, fhs_larvae, 130,
                                       5, 'ssp585', miroc_temps3, 
                                       miroc_salts3, larval_formula)

# Combine into one data frame
df_fhslarvae3_miroc585 <- list(preds_fhslarvae3_miroc585[[1]], preds_fhslarvae3_miroc585[[2]],
                               preds_fhslarvae3_miroc585[[3]], preds_fhslarvae3_miroc585[[4]],
                               preds_fhslarvae3_miroc585[[5]], preds_fhslarvae3_miroc585[[6]],
                               preds_fhslarvae3_miroc585[[7]], preds_fhslarvae3_miroc585[[8]],
                               preds_fhslarvae3_miroc585[[9]], preds_fhslarvae3_miroc585[[10]],
                               preds_fhslarvae3_miroc585[[11]], preds_fhslarvae3_miroc585[[12]],
                               preds_fhslarvae3_miroc585[[13]], preds_fhslarvae3_miroc585[[14]],
                               preds_fhslarvae3_miroc585[[15]], preds_fhslarvae3_miroc585[[16]],
                               preds_fhslarvae3_miroc585[[17]], preds_fhslarvae3_miroc585[[18]],
                               preds_fhslarvae3_miroc585[[19]], preds_fhslarvae3_miroc585[[20]],
                               preds_fhslarvae3_miroc585[[21]], preds_fhslarvae3_miroc585[[22]],
                               preds_fhslarvae3_miroc585[[23]], preds_fhslarvae3_miroc585[[24]],
                               preds_fhslarvae3_miroc585[[25]], preds_fhslarvae3_miroc585[[26]], 
                               preds_fhslarvae3_miroc585[[27]], preds_fhslarvae3_miroc585[[28]], 
                               preds_fhslarvae3_miroc585[[29]], preds_fhslarvae3_miroc585[[30]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 


# Generate average prediction from all predictions
x <- grepl("pred", names(df_fhslarvae3_miroc585), fixed = T)
df_fhslarvae_avg3_miroc585 <- data.frame(lat = df_fhslarvae3_miroc585$lat, 
                                         lon = df_fhslarvae3_miroc585$lon, 
                                         dist = df_fhslarvae3_miroc585$dist,
                                         avg_pred = rowSums(df_fhslarvae3_miroc585[, x])/30)
saveRDS(df_fhslarvae_avg3_miroc585, file = here("data", "df_fhslarvae_avg3_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae_avg3_miroc585, "Forecasted Distribution 2070 - 2099 \n miroc SSP585")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_miroc_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


rm(miroc_temps1, miroc_temps2, miroc_temps3,
   miroc_salts1, miroc_salts2, miroc_salts3)


### Average Predictions ------------------------------------------------------------------------------------------------------------
# Change to have same scale for all
grid_predict_egg <- function(grid, title){
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
        zlim = c(0, 2522),
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
             zlim = c(0, 2522),
             legend.args = list("Avg. Predicted \n Occurrence",
                                side = 2, cex = 1))
}
grid_predict_larvae <- function(grid, title){
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
        zlim = c(0, 1695),
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
             zlim = c(0, 1695),
             legend.args = list("Avg. Predicted \n Occurrence",
                                side = 2, cex = 1))
}

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
grid_predict_egg(df_fhsegg_final1, "Forecasted Distribution 2015 - 2039")
dev.copy(jpeg,
         here('results/flathead_forecast/fhsegg_avgs',
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
grid_predict_larvae(df_fhslarvae_final1, "Forecasted Distribution 2015 - 2039")
dev.copy(jpeg,
         here('results/flathead_forecast/fhslarvae_avgs',
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
grid_predict_egg(df_fhsegg_final2, "Forecasted Distribution 2040 - 2069")
dev.copy(jpeg,
         here('results/flathead_forecast/fhsegg_avgs',
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
grid_predict_larvae(df_fhslarvae_final2, "Forecasted Distribution 2040 - 2069")
dev.copy(jpeg,
         here('results/flathead_forecast/fhslarvae_avgs',
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
grid_predict_egg(df_fhsegg_final3, "Forecasted Distribution 2070 - 2099")
dev.copy(jpeg,
         here('results/flathead_forecast/fhsegg_avgs',
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
grid_predict_larvae(df_fhslarvae_final3, "Forecasted Distribution 2070 - 2099")
dev.copy(jpeg,
         here('results/flathead_forecast/fhslarvae_avgs',
              'flathead_larvae_avg3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

#### Make GIFs -----------------------------------------------------------------------------------------------------------------------
library(magick)
base_dir <- getwd()
dir_out <- file.path(base_dir, 'results', 'flathead_forecast', 'fhsegg_avgs')

imgs <- list.files(dir_out, full.names = T)
img_list <- lapply(imgs, image_read)

img_joined <- image_join(img_list)

img_animated <- image_animate(img_joined, fps = 1)

img_animated

image_write(image = img_animated,
            path = here('results', 'flathead_forecast', "fhsegg_avgs.gif"))


dir_out2 <- file.path(base_dir, 'results', 'flathead_forecast', 'fhslarvae_avgs')

imgs2 <- list.files(dir_out2, full.names = T)
img_list2 <- lapply(imgs2, image_read)

img_joined2 <- image_join(img_list2)

img_animated2 <- image_animate(img_joined2, fps = 1)

img_animated2

image_write(image = img_animated2,
            path = here('results', 'flathead_forecast', "fhslarvae_avgs.gif"))
