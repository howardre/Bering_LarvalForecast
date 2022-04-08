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
source(here('code/functions', 'distance_function.R'))
source(here('code/functions', 'grid_predict.R'))
source(here('code/functions', 'get_preds.R'))
source(here('code/functions', 'pred_avgs.R'))

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
  gam(catch ~ s(year, bs = 're') +
        s(doy, k = 8) +
        s(lon, lat) +
        s(roms_temperature, k = 6) +
        s(roms_salinity, k = 6) +
        s(doy, by = mean_temp, k = 6), # phenology
      data = data,
      family = tw(link = 'log'),
      method = 'REML')
}

formula_geog <- function(data){
  gam(catch ~ s(year, bs = 're') +
        s(doy, k = 8) +
        s(lon, lat) +
        s(roms_temperature, k = 6) +
        s(roms_salinity, k = 6) +
        s(lat, lon, by = mean_temp), # geography
      data = data,
      family = tw(link = 'log'),
      method = 'REML')
}

# Load ROMS temperature means and forecast
roms_temps <- readRDS(here('data', 'roms_temps.rds'))

# Load fish data
pk_larvae <- load_data('pk_larvae.rds', pk_larvae, roms_temps)
fhs_egg <- load_data('fhs_egg.rds', fhs_egg, roms_temps)
fhs_larvae <- load_data('fhs_larvae.rds', fhs_larvae, roms_temps)
akp_egg <- load_data('akp_egg.rds', akp_egg, roms_temps)
yfs_larvae <- load_data('yfs_larvae.rds', yfs_larvae, roms_temps)
pcod_larvae <- load_data('pcod_larvae.rds', pcod_larvae, roms_temps)
nrs_larvae <- load_data('nrs_larvae.rds', nrs_larvae, roms_temps)

### Pollock Eggs --------------------------------------------------------------------------------------------------------------------------
pk_egg <- load_data('pk_egg.rds', pk_egg, roms_temps)
pkegg_formula <- formula_pheno(pk_egg)

#### Forecast and average into 3 time periods ---------------------------------------------------------------------------------------------
##### CESM 126 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
cesm_temps1 <- readRDS(here('data', 'cesm_forecast_temp1.rds'))
cesm_salts1 <- readRDS(here('data', 'cesm_forecast_salt1.rds'))

df_pkegg1_cesm126 <- pred_list1(2015:2039, pk_egg, 134,
                                5, 'ssp126', cesm_temps1, 
                                cesm_salts1, pkegg_formula)
avgs_pkegg1_cesm126 <- pred_avgs1(df_pkegg1_cesm126)
saveRDS(avgs_pkegg1_cesm126, file = here("data", "avgs_pkegg1_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(avgs_pkegg1_cesm126, 
             "Forecasted Distribution 2015 - 2039 \n CESM SSP126")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_egg_cesm_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
cesm_temps2 <- readRDS(here('data', 'cesm_forecast_temp2.rds'))
cesm_salts2 <- readRDS(here('data', 'cesm_forecast_salt2.rds'))

df_pkegg2_cesm126 <- pred_list2(2040:2069, pk_egg, 134,
                                5, 'ssp126', cesm_temps2, 
                                cesm_salts2, pkegg_formula)
avgs_pkegg2_cesm126 <- pred_avgs2(df_pkegg2_cesm126)
saveRDS(avgs_pkegg2_cesm126, file = here("data", "avgs_pkegg2_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(avgs_pkegg2_cesm126, 
             "Forecasted Distribution 2040 - 2069 \n CESM SSP126")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_egg_cesm_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
cesm_temps3 <- readRDS(here('data', 'cesm_forecast_temp3.rds'))
cesm_salts3 <- readRDS(here('data', 'cesm_forecast_salt3.rds'))

df_pkegg3_cesm126 <- pred_list2(2070:2099, pk_egg, 134,
                                5, 'ssp126', cesm_temps3, 
                                cesm_salts3, pkegg_formula)
avgs_pkegg3_cesm126 <- pred_avgs2(df_pkegg3_cesm126)
saveRDS(avgs_pkegg3_cesm126, file = here("data", "avgs_pkegg3_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(avgs_pkegg3_cesm126, 
             "Forecasted Distribution 2070 - 2099 \n CESM SSP126")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_egg_cesm_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### CESM 585 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
df_pkegg1_cesm585 <- pred_list1(2015:2039, pk_egg, 134,
                                5, 'ssp585', cesm_temps1, 
                                cesm_salts1, pkegg_formula)
avgs_pkegg1_cesm585 <- pred_avgs1(df_pkegg1_cesm585)
saveRDS(avgs_pkegg1_cesm585, file = here("data", "avgs_pkegg1_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(avgs_pkegg1_cesm585, 
             "Forecasted Distribution 2015 - 2039 \n CESM SSP585")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_egg_cesm_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
df_pkegg2_cesm585 <- pred_list2(2040:2069, pk_egg, 134,
                                5, 'ssp585', cesm_temps2, 
                                cesm_salts2, pkegg_formula)
avgs_pkegg2_cesm585 <- pred_avgs2(df_pkegg2_cesm585)
saveRDS(avgs_pkegg2_cesm585, file = here("data", "avgs_pkegg2_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(avgs_pkegg2_cesm585, 
             "Forecasted Distribution 2040 - 2069 \n CESM SSP585")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_egg_cesm_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
df_pkegg3_cesm585 <- pred_list2(2070:2099, pk_egg, 134,
                                5, 'ssp585', cesm_temps3, 
                                cesm_salts3, pkegg_formula)
avgs_pkegg3_cesm585 <- pred_avgs2(df_pkegg3_cesm585)
saveRDS(avgs_pkegg3_cesm585, file = here("data", "avgs_pkegg3_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(avgs_pkegg3_cesm585, 
             "Forecasted Distribution 2070 - 2099 \n CESM SSP585")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_egg_cesm_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

rm(cesm_temps1, cesm_temps2, cesm_temps3,
   cesm_salts1, cesm_salts2, cesm_salts3)

##### GFDL 126 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
gfdl_temps1 <- readRDS(here('data', 'gfdl_forecast_temp1.rds'))
gfdl_salts1 <- readRDS(here('data', 'gfdl_forecast_salt1.rds'))

df_pkegg1_gfdl126 <- pred_list1(2015:2039, pk_egg, 134,
                                5, 'ssp126', gfdl_temps1, 
                                gfdl_salts1, pkegg_formula)
avgs_pkegg1_gfdl126 <- pred_avgs1(df_pkegg1_gfdl126)
saveRDS(avgs_pkegg1_gfdl126, file = here("data", "avgs_pkegg1_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(avgs_pkegg1_gfdl126, 
             "Forecasted Distribution 2015 - 2039 \n GFDL SSP126")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_egg_gfdl_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
gfdl_temps2 <- readRDS(here('data', 'gfdl_forecast_temp2.rds'))
gfdl_salts2 <- readRDS(here('data', 'gfdl_forecast_salt2.rds'))

df_pkegg2_gfdl126 <- pred_list2(2040:2069, pk_egg, 134,
                                5, 'ssp126', gfdl_temps2, 
                                gfdl_salts2, pkegg_formula)
avgs_pkegg2_gfdl126 <- pred_avgs2(df_pkegg2_gfdl126)
saveRDS(avgs_pkegg2_gfdl126, file = here("data", "avgs_pkegg2_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(avgs_pkegg2_gfdl126, 
             "Forecasted Distribution 2040 - 2069 \n GFDL SSP126")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_egg_gfdl_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
gfdl_temps3 <- readRDS(here('data', 'gfdl_forecast_temp3.rds'))
gfdl_salts3 <- readRDS(here('data', 'gfdl_forecast_salt3.rds'))

df_pkegg3_gfdl126 <- pred_list2(2070:2099, pk_egg, 134,
                                5, 'ssp126', gfdl_temps3, 
                                gfdl_salts3, pkegg_formula)
avgs_pkegg3_gfdl126 <- pred_avgs2(df_pkegg3_gfdl126)
saveRDS(avgs_pkegg3_gfdl126, file = here("data", "avgs_pkegg3_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(avgs_pkegg3_gfdl126, 
             "Forecasted Distribution 2070 - 2099 \n GFDL SSP126")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_egg_gfdl_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### GFDL 585 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
df_pkegg1_gfdl585 <- pred_list1(2015:2039, pk_egg, 134,
                                5, 'ssp585', gfdl_temps1, 
                                gfdl_salts1, pkegg_formula)
avgs_pkegg1_gfdl585 <- pred_avgs1(df_pkegg1_gfdl585)
saveRDS(avgs_pkegg1_gfdl585, file = here("data", "avgs_pkegg1_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(avgs_pkegg1_gfdl585, 
             "Forecasted Distribution 2015 - 2039 \n GFDL SSP585")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_egg_gfdl_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
df_pkegg2_gfdl585 <- pred_list2(2040:2069, pk_egg, 134,
                                5, 'ssp585', gfdl_temps2, 
                                gfdl_salts2, pkegg_formula)
avgs_pkegg2_gfdl585 <- pred_avgs2(df_pkegg2_gfdl585)
saveRDS(avgs_pkegg2_gfdl585, file = here("data", "avgs_pkegg2_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(avgs_pkegg2_gfdl585, 
             "Forecasted Distribution 2040 - 2069 \n GFDL SSP585")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_egg_gfdl_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
df_pkegg3_gfdl585 <- pred_list2(2070:2099, pk_egg, 134,
                                5, 'ssp585', gfdl_temps3, 
                                gfdl_salts3, pkegg_formula)
avgs_pkegg3_gfdl585 <- pred_avgs2(df_pkegg3_gfdl585)
saveRDS(avgs_pkegg3_gfdl585, file = here("data", "avgs_pkegg3_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(avgs_pkegg3_gfdl585, 
             "Forecasted Distribution 2070 - 2099 \n GFDL SSP585")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_egg_gfdl_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

rm(gfdl_temps1, gfdl_temps2, gfdl_temps3,
   gfdl_salts1, gfdl_salts2, gfdl_salts3)

##### MIROC 126 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
miroc_temps1 <- readRDS(here('data', 'miroc_forecast_temp1.rds'))
miroc_salts1 <- readRDS(here('data', 'miroc_forecast_salt1.rds'))

df_pkegg1_miroc126 <- pred_list1(2015:2039, pk_egg, 134,
                                5, 'ssp126', miroc_temps1, 
                                miroc_salts1, pkegg_formula)
avgs_pkegg1_miroc126 <- pred_avgs1(df_pkegg1_miroc126)
saveRDS(avgs_pkegg1_miroc126, file = here("data", "avgs_pkegg1_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(avgs_pkegg1_miroc126, 
             "Forecasted Distribution 2015 - 2039 \n MIROC SSP126")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_egg_miroc_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
miroc_temps2 <- readRDS(here('data', 'miroc_forecast_temp2.rds'))
miroc_salts2 <- readRDS(here('data', 'miroc_forecast_salt2.rds'))

df_pkegg2_miroc126 <- pred_list2(2040:2069, pk_egg, 134,
                                5, 'ssp126', miroc_temps2, 
                                miroc_salts2, pkegg_formula)
avgs_pkegg2_miroc126 <- pred_avgs2(df_pkegg2_miroc126)
saveRDS(avgs_pkegg2_miroc126, file = here("data", "avgs_pkegg2_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(avgs_pkegg2_miroc126, 
             "Forecasted Distribution 2040 - 2069 \n MIROC SSP126")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_egg_miroc_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
miroc_temps3 <- readRDS(here('data', 'miroc_forecast_temp3.rds'))
miroc_salts3 <- readRDS(here('data', 'miroc_forecast_salt3.rds'))

df_pkegg3_miroc126 <- pred_list2(2070:2099, pk_egg, 134,
                                5, 'ssp126', miroc_temps3, 
                                miroc_salts3, pkegg_formula)
avgs_pkegg3_miroc126 <- pred_avgs2(df_pkegg3_miroc126)
saveRDS(avgs_pkegg3_miroc126, file = here("data", "avgs_pkegg3_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(avgs_pkegg3_miroc126, 
             "Forecasted Distribution 2070 - 2099 \n MIROC SSP126")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_egg_miroc_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### MIROC 585 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
df_pkegg1_miroc585 <- pred_list1(2015:2039, pk_egg, 134,
                                5, 'ssp585', miroc_temps1, 
                                miroc_salts1, pkegg_formula)
avgs_pkegg1_miroc585 <- pred_avgs1(df_pkegg1_miroc585)
saveRDS(avgs_pkegg1_miroc585, file = here("data", "avgs_pkegg1_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(avgs_pkegg1_miroc585, 
             "Forecasted Distribution 2015 - 2039 \n MIROC SSP585")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_egg_miroc_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
df_pkegg2_miroc585 <- pred_list2(2040:2069, pk_egg, 134,
                                5, 'ssp585', miroc_temps2, 
                                miroc_salts2, pkegg_formula)
avgs_pkegg2_miroc585 <- pred_avgs2(df_pkegg2_miroc585)
saveRDS(avgs_pkegg2_miroc585, file = here("data", "avgs_pkegg2_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(avgs_pkegg2_miroc585, 
             "Forecasted Distribution 2040 - 2069 \n MIROC SSP585")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_egg_miroc_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
df_pkegg3_miroc585 <- predict_cells(2070:2099, pk_egg, 134,
                                    5, 'ssp585', miroc_temps3,
                                    miroc_salts3, pkegg_formula)
avgs_pkegg3_miroc585 <- pred_avgs2(df_pkegg3_miroc585)
saveRDS(avgs_pkegg3_miroc585, file = here("data", "avgs_pkegg3_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pkegg3_miroc585, 
             "Forecasted Distribution 2070 - 2099 \n MIROC SSP585")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_egg_miroc_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

rm(miroc_temps1, miroc_temps2, miroc_temps3,
   miroc_salts1, miroc_salts2, miroc_salts3)

##### Multi-panel figure
tiff(here('results/pollock_forecast',
          'pollockegg_multipanel.tiff'),
     units = "in",
     width = 65,
     height = 90,
     res = 300)
par(mfrow = c(18, 3),
    mar = c(1, 1, 1, 1),
    family = "serif")
grid_predict(avgs_pkegg1_cesm126, 
             "Forecasted Distribution 2015 - 2039 \n CESM SSP126")
grid_predict(avgs_pkegg2_cesm126, 
             "Forecasted Distribution 2040 - 2069 \n CESM SSP126")
grid_predict(avgs_pkegg3_cesm126, 
             "Forecasted Distribution 2070 - 2099 \n CESM SSP126")
grid_predict(avgs_pkegg1_cesm585, 
             "Forecasted Distribution 2015 - 2039 \n CESM SSP585")
grid_predict(avgs_pkegg2_cesm585, 
             "Forecasted Distribution 2040 - 2069 \n CESM SSP585")
grid_predict(avgs_pkegg3_cesm585, 
             "Forecasted Distribution 2070 - 2099 \n CESM SSP585")
grid_predict(avgs_pkegg1_gfdl126, 
             "Forecasted Distribution 2015 - 2039 \n GFDL SSP126")
grid_predict(avgs_pkegg2_gfdl126, 
             "Forecasted Distribution 2040 - 2069 \n GFDL SSP126")
grid_predict(avgs_pkegg3_gfdl126, 
             "Forecasted Distribution 2070 - 2099 \n GFDL SSP126")
grid_predict(avgs_pkegg1_gfdl585, 
             "Forecasted Distribution 2015 - 2039 \n GFDL SSP585")
grid_predict(avgs_pkegg2_gfdl585, 
             "Forecasted Distribution 2040 - 2069 \n GFDL SSP585")
grid_predict(avgs_pkegg3_gfdl585, 
             "Forecasted Distribution 2070 - 2099 \n GFDL SSP585")
grid_predict(avgs_pkegg1_miroc126, 
             "Forecasted Distribution 2015 - 2039 \n MIROC SSP126")
grid_predict(avgs_pkegg2_miroc126, 
             "Forecasted Distribution 2040 - 2069 \n MIROC SSP126")
grid_predict(avgs_pkegg3_miroc126, 
             "Forecasted Distribution 2070 - 2099 \n MIROC SSP126")
grid_predict(avgs_pkegg1_miroc585, 
             "Forecasted Distribution 2015 - 2039 \n MIROC SSP585")
grid_predict(avgs_pkegg2_miroc585, 
             "Forecasted Distribution 2040 - 2069 \n MIROC SSP585")
grid_predict(avgs_pkegg3_miroc585, 
             "Forecasted Distribution 2070 - 2099 \n MIROC SSP585")



##### Averages ---------------------------------------------------------------------------------------------------------------------------
df_pkegg_merged1 <- list(avgs_pkegg1_cesm126, avgs_pkegg1_cesm585,
                         avgs_pkegg1_gfdl126, avgs_pkegg1_gfdl585,
                         avgs_pkegg1_miroc126, avgs_pkegg1_miroc585) %>%
  reduce(inner_join, by = c("lon", "lat", "dist"))

x <- grepl("pred", names(df_pkegg_merged1), fixed = T)
df_pkegg_final1 <- data.frame(lat = df_pkegg_merged1$lat,
                              lon = df_pkegg_merged1$lon,
                              avg_pred = rowSums(df_pkegg_merged1[, x])/6)
df_pkegg_final1$pred_scaled <- rescale(df_pkegg_final1$avg_pred)

windows(width = 6, height = 6, family = "serif")
grid_predict(df_pkegg_final1, "Forecasted Distribution 2015 - 2039")
dev.copy(jpeg,
         here('results/pollock_forecast/pkegg_avgs',
              'pollock_egg_avg1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()