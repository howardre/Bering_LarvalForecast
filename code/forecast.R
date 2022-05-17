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
        s(doy, k = 9) +
        s(lon, lat) +
        s(roms_temperature, k = 9) +
        s(roms_salinity, k = 9) +
        s(doy, by = mean_temp), # phenology
      data = data,
      family = tw(link = "log"),
      method = 'REML')
}

formula_geog <- function(data){
  gam(catch ~ s(year, bs = 're') +
        s(doy, k = 9) +
        s(lon, lat) +
        s(roms_temperature, k = 9) +
        s(roms_salinity, k = 9) +
        s(lat, lon, by = mean_temp), # geography
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

df_pkegg1_cesm126 <- predict_cells(2015:2039, pk_egg, 134,
                                   5, 'ssp126', cesm_temps1,
                                   cesm_salts1, pkegg_formula)
saveRDS(df_pkegg1_cesm126, file = here("data", "df_pkegg1_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pkegg1_cesm126, 
             "Forecasted Distribution 2015 - 2039 \n CESM SSP1-2.6")
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

df_pkegg2_cesm126 <- predict_cells(2040:2069, pk_egg, 134,
                                5, 'ssp126', cesm_temps2, 
                                cesm_salts2, pkegg_formula)
saveRDS(df_pkegg2_cesm126, file = here("data", "df_pkegg2_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pkegg2_cesm126, 
             "Forecasted Distribution 2040 - 2069 \n CESM SSP1-2.6")
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

df_pkegg3_cesm126 <- predict_cells(2070:2099, pk_egg, 134,
                                5, 'ssp126', cesm_temps3, 
                                cesm_salts3, pkegg_formula)
saveRDS(df_pkegg3_cesm126, file = here("data", "df_pkegg3_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pkegg3_cesm126, 
             "Forecasted Distribution 2070 - 2099 \n CESM SSP1-2.6")
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
df_pkegg1_cesm585 <- predict_cells(2015:2039, pk_egg, 134,
                                5, 'ssp585', cesm_temps1, 
                                cesm_salts1, pkegg_formula)
saveRDS(df_pkegg1_cesm585, file = here("data", "df_pkegg1_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pkegg1_cesm585, 
             "Forecasted Distribution 2015 - 2039 \n CESM SSP5-8.5")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_egg_cesm_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
df_pkegg2_cesm585 <- predict_cells(2040:2069, pk_egg, 134,
                                5, 'ssp585', cesm_temps2, 
                                cesm_salts2, pkegg_formula)
saveRDS(df_pkegg2_cesm585, file = here("data", "df_pkegg2_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pkegg2_cesm585, 
             "Forecasted Distribution 2040 - 2069 \n CESM SSP5-8.5")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_egg_cesm_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
df_pkegg3_cesm585 <- predict_cells(2070:2099, pk_egg, 134,
                                5, 'ssp585', cesm_temps3, 
                                cesm_salts3, pkegg_formula)
saveRDS(df_pkegg3_cesm585, file = here("data", "df_pkegg3_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pkegg3_cesm585, 
             "Forecasted Distribution 2070 - 2099 \n CESM SSP5-8.5")
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

df_pkegg1_gfdl126 <- predict_cells(2015:2039, pk_egg, 134,
                                5, 'ssp126', gfdl_temps1, 
                                gfdl_salts1, pkegg_formula)
saveRDS(df_pkegg1_gfdl126, file = here("data", "df_pkegg1_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pkegg1_gfdl126, 
             "Forecasted Distribution 2015 - 2039 \n GFDL SSP1-2.6")
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

df_pkegg2_gfdl126 <- predict_cells(2040:2069, pk_egg, 134,
                                5, 'ssp126', gfdl_temps2, 
                                gfdl_salts2, pkegg_formula)
saveRDS(df_pkegg2_gfdl126, file = here("data", "df_pkegg2_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pkegg2_gfdl126, 
             "Forecasted Distribution 2040 - 2069 \n GFDL SSP1-2.6")
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

df_pkegg3_gfdl126 <- predict_cells(2070:2099, pk_egg, 134,
                                5, 'ssp126', gfdl_temps3, 
                                gfdl_salts3, pkegg_formula)
saveRDS(df_pkegg3_gfdl126, file = here("data", "df_pkegg3_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pkegg3_gfdl126, 
             "Forecasted Distribution 2070 - 2099 \n GFDL SSP1-2.6")
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
df_pkegg1_gfdl585 <- predict_cells(2015:2039, pk_egg, 134,
                                5, 'ssp585', gfdl_temps1, 
                                gfdl_salts1, pkegg_formula)
saveRDS(df_pkegg1_gfdl585, file = here("data", "df_pkegg1_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pkegg1_gfdl585, 
             "Forecasted Distribution 2015 - 2039 \n GFDL SSP5-8.5")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_egg_gfdl_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
df_pkegg2_gfdl585 <- predict_cells(2040:2069, pk_egg, 134,
                                5, 'ssp585', gfdl_temps2, 
                                gfdl_salts2, pkegg_formula)
saveRDS(df_pkegg2_gfdl585, file = here("data", "df_pkegg2_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pkegg2_gfdl585, 
             "Forecasted Distribution 2040 - 2069 \n GFDL SSP5-8.5")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_egg_gfdl_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
df_pkegg3_gfdl585 <- predict_cells(2070:2099, pk_egg, 134,
                                5, 'ssp585', gfdl_temps3, 
                                gfdl_salts3, pkegg_formula)
saveRDS(df_pkegg3_gfdl585, file = here("data", "df_pkegg3_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pkegg3_gfdl585, 
             "Forecasted Distribution 2070 - 2099 \n GFDL SSP5-8.5")
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

df_pkegg1_miroc126 <- predict_cells(2015:2039, pk_egg, 134,
                                5, 'ssp126', miroc_temps1, 
                                miroc_salts1, pkegg_formula)
saveRDS(df_pkegg1_miroc126, file = here("data", "df_pkegg1_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pkegg1_miroc126, 
             "Forecasted Distribution 2015 - 2039 \n MIROC SSP1-2.6")
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

df_pkegg2_miroc126 <- predict_cells(2040:2069, pk_egg, 134,
                                5, 'ssp126', miroc_temps2, 
                                miroc_salts2, pkegg_formula)
saveRDS(df_pkegg2_miroc126, file = here("data", "df_pkegg2_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pkegg2_miroc126, 
             "Forecasted Distribution 2040 - 2069 \n MIROC SSP1-2.6")
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

df_pkegg3_miroc126 <- predict_cells(2070:2099, pk_egg, 134,
                                5, 'ssp126', miroc_temps3, 
                                miroc_salts3, pkegg_formula)
saveRDS(df_pkegg3_miroc126, file = here("data", "df_pkegg3_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pkegg3_miroc126, 
             "Forecasted Distribution 2070 - 2099 \n MIROC SSP1-2.6")
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
df_pkegg1_miroc585 <- predict_cells(2015:2039, pk_egg, 134,
                                5, 'ssp585', miroc_temps1, 
                                miroc_salts1, pkegg_formula)
saveRDS(df_pkegg1_miroc585, file = here("data", "df_pkegg1_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pkegg1_miroc585, 
             "Forecasted Distribution 2015 - 2039 \n MIROC SSP5-8.5")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_egg_miroc_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
df_pkegg2_miroc585 <- predict_cells(2040:2069, pk_egg, 134,
                                5, 'ssp585', miroc_temps2, 
                                miroc_salts2, pkegg_formula)
saveRDS(df_pkegg2_miroc585, file = here("data", "df_pkegg2_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pkegg2_miroc585, 
             "Forecasted Distribution 2040 - 2069 \n MIROC SSP5-8.5")
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
saveRDS(df_pkegg3_miroc585, file = here("data", "df_pkegg3_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pkegg3_miroc585, 
             "Forecasted Distribution 2070 - 2099 \n MIROC SSP5-8.5")
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

graphics.off()

##### Multi-panel figure -------------------------------------------------------------------------------------------------------------------------
tiff(here('results/pollock_forecast',
          'pollockegg_multipanel.tiff'),
     units = "in",
     width = 45,
     height = 90,
     res = 300)
par(mfrow = c(6, 3),
    mar = c(11, 12, 5, 1.5) + 0.1,
    oma = c(3, 25, 15, 1),
    mgp = c(10, 4, 0),
    family = "serif")
grid_multipanel(df_pkegg1_cesm126)
grid_multipanel(df_pkegg2_cesm126)
grid_multipanel(df_pkegg3_cesm126)
grid_multipanel(df_pkegg1_cesm585)
grid_multipanel(df_pkegg2_cesm585)
grid_multipanel(df_pkegg3_cesm585)
grid_multipanel(df_pkegg1_gfdl126)
grid_multipanel(df_pkegg2_gfdl126)
grid_multipanel(df_pkegg3_gfdl126)
grid_multipanel(df_pkegg1_gfdl585)
grid_multipanel(df_pkegg2_gfdl585)
grid_multipanel(df_pkegg3_gfdl585)
grid_multipanel(df_pkegg1_miroc126)
grid_multipanel(df_pkegg2_miroc126)
grid_multipanel(df_pkegg3_miroc126)
grid_multipanel(df_pkegg1_miroc585)
grid_multipanel(df_pkegg2_miroc585)
grid_multipanel(df_pkegg3_miroc585)
mtext("CESM SSP1-2.6", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.92)
mtext("CESM SSP5-8.5", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.76)
mtext("GFDL1-2.6", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.585)
mtext("GFDL5-8.5", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.42)
mtext("MIROC1-2.6", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.26)
mtext("MIROC5-8.5", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.095)
mtext("2015-2039", 
      side = 3, 
      line = 2, 
      outer = TRUE, 
      cex = 7,
      at = 0.17)
mtext("2040-2069", 
      side = 3, 
      line = 2, 
      outer = TRUE, 
      cex = 7,
      at = 0.51)
mtext("2070-2099", 
      side = 3, 
      line = 2, 
      outer = TRUE, 
      cex = 7,
      at = 0.84)
dev.off()

##### Averages ---------------------------------------------------------------------------------------------------------------------------
# 2015 - 2039
df_pkegg_merged1 <- list(df_pkegg1_cesm126, df_pkegg1_cesm585,
                         df_pkegg1_gfdl126, df_pkegg1_gfdl585,
                         df_pkegg1_miroc126, df_pkegg1_miroc585) 

avg_pkegg_merged1 <- predict_avgs(df_pkegg_merged1)

windows(width = 6, height = 6, family = "serif")
grid_predict(avg_pkegg_merged1, "Forecasted Distribution 2015 - 2039")
dev.copy(jpeg,
         here('results/pollock_forecast/pkegg_avgs',
              'pollock_egg_avg1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

# 2040-2069
df_pkegg_merged2 <- list(df_pkegg2_cesm126, df_pkegg2_cesm585,
                         df_pkegg2_gfdl126, df_pkegg2_gfdl585,
                         df_pkegg2_miroc126, df_pkegg2_miroc585) 

avg_pkegg_merged2 <- predict_avgs(df_pkegg_merged2)

windows(width = 6, height = 6, family = "serif")
grid_predict(avg_pkegg_merged2, "Forecasted Distribution 2040 - 2069")
dev.copy(jpeg,
         here('results/pollock_forecast/pkegg_avgs',
              'pollock_egg_avg2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

# 2070-2099
df_pkegg_merged3 <- list(df_pkegg3_cesm126, df_pkegg3_cesm585,
                         df_pkegg3_gfdl126, df_pkegg3_gfdl585,
                         df_pkegg3_miroc126, df_pkegg3_miroc585) 

avg_pkegg_merged3 <- predict_avgs(df_pkegg_merged3)

windows(width = 6, height = 6, family = "serif")
grid_predict(avg_pkegg_merged3, "Forecasted Distribution 2070 - 2099")
dev.copy(jpeg,
         here('results/pollock_forecast/pkegg_avgs',
              'pollock_egg_avg3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

##### GIFs -------------------------------------------------------------------------------------------------------------------------
pkegg_dir_out <- file.path(base_dir, 'results', 'pollock_forecast', 'pkegg_avgs')
pkegg_imgs <- list.files(pkegg_dir_out, full.names = T)
pkegg_img_list <- lapply(pkegg_imgs, image_read)
pkegg_img_joined <- image_join(pkegg_img_list)
pkegg_img_animated <- image_animate(pkegg_img_joined, fps = 1)
image_write(image = pkegg_img_animated,
            path = here('results', 'pollock_forecast', "pkegg_avgs.gif"))

##### Clear environment -------------------------------------------------------------------------------------------------------------------------
rm(df_pkegg1_cesm126, df_pkegg1_cesm585,
   df_pkegg1_gfdl126, df_pkegg1_gfdl585,
   df_pkegg1_miroc126, df_pkegg1_miroc585,
   df_pkegg2_cesm126, df_pkegg2_cesm585,
   df_pkegg2_gfdl126, df_pkegg2_gfdl585,
   df_pkegg2_miroc126, df_pkegg2_miroc585,
   df_pkegg3_cesm126, df_pkegg3_cesm585,
   df_pkegg3_gfdl126, df_pkegg3_gfdl585,
   df_pkegg3_miroc126, df_pkegg3_miroc585,
   avg_pkegg_merged1, avg_pkegg_merged2,
   avg_pkegg_merged3, df_pkegg_merged1,
   df_pkegg_merged2, df_pkegg_merged3,
   pk_egg, pkegg_formula, pkegg_img_animated,
   pkegg_dir_out, pkegg_img_joined, 
   pkegg_img_list, pkegg_imgs)


### Pollock Larvae --------------------------------------------------------------------------------------------------------------------------
pk_larvae <- load_data('pk_larvae.rds', pk_larvae, roms_temps)
pklarvae_formula <- formula_geog(pk_larvae)

#### Forecast and average into 3 time periods ---------------------------------------------------------------------------------------------
##### CESM 126 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
cesm_temps1 <- readRDS(here('data', 'cesm_forecast_temp1.rds'))
cesm_salts1 <- readRDS(here('data', 'cesm_forecast_salt1.rds'))

df_pklarvae1_cesm126 <- predict_cells(2015:2039, pk_larvae, 134,
                                   5, 'ssp126', cesm_temps1,
                                   cesm_salts1, pklarvae_formula)
saveRDS(df_pklarvae1_cesm126, file = here("data", "df_pklarvae1_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pklarvae1_cesm126, 
             "Forecasted Distribution 2015 - 2039 \n CESM SSP1-2.6")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_larvae_cesm_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
cesm_temps2 <- readRDS(here('data', 'cesm_forecast_temp2.rds'))
cesm_salts2 <- readRDS(here('data', 'cesm_forecast_salt2.rds'))

df_pklarvae2_cesm126 <- predict_cells(2040:2069, pk_larvae, 134,
                                   5, 'ssp126', cesm_temps2, 
                                   cesm_salts2, pklarvae_formula)
saveRDS(df_pklarvae2_cesm126, file = here("data", "df_pklarvae2_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pklarvae2_cesm126, 
             "Forecasted Distribution 2040 - 2069 \n CESM SSP1-2.6")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_larvae_cesm_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
cesm_temps3 <- readRDS(here('data', 'cesm_forecast_temp3.rds'))
cesm_salts3 <- readRDS(here('data', 'cesm_forecast_salt3.rds'))

df_pklarvae3_cesm126 <- predict_cells(2070:2099, pk_larvae, 134,
                                   5, 'ssp126', cesm_temps3, 
                                   cesm_salts3, pklarvae_formula)
saveRDS(df_pklarvae3_cesm126, file = here("data", "df_pklarvae3_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pklarvae3_cesm126, 
             "Forecasted Distribution 2070 - 2099 \n CESM SSP1-2.6")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_larvae_cesm_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### CESM 585 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
df_pklarvae1_cesm585 <- predict_cells(2015:2039, pk_larvae, 134,
                                   5, 'ssp585', cesm_temps1, 
                                   cesm_salts1, pklarvae_formula)
saveRDS(df_pklarvae1_cesm585, file = here("data", "df_pklarvae1_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pklarvae1_cesm585, 
             "Forecasted Distribution 2015 - 2039 \n CESM SSP5-8.5")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_larvae_cesm_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
df_pklarvae2_cesm585 <- predict_cells(2040:2069, pk_larvae, 134,
                                   5, 'ssp585', cesm_temps2, 
                                   cesm_salts2, pklarvae_formula)
saveRDS(df_pklarvae2_cesm585, file = here("data", "df_pklarvae2_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pklarvae2_cesm585, 
             "Forecasted Distribution 2040 - 2069 \n CESM SSP5-8.5")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_larvae_cesm_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
df_pklarvae3_cesm585 <- predict_cells(2070:2099, pk_larvae, 134,
                                   5, 'ssp585', cesm_temps3, 
                                   cesm_salts3, pklarvae_formula)
saveRDS(df_pklarvae3_cesm585, file = here("data", "df_pklarvae3_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pklarvae3_cesm585, 
             "Forecasted Distribution 2070 - 2099 \n CESM SSP5-8.5")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_larvae_cesm_ssp585_3.jpg'),
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

df_pklarvae1_gfdl126 <- predict_cells(2015:2039, pk_larvae, 134,
                                   5, 'ssp126', gfdl_temps1, 
                                   gfdl_salts1, pklarvae_formula)
saveRDS(df_pklarvae1_gfdl126, file = here("data", "df_pklarvae1_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pklarvae1_gfdl126, 
             "Forecasted Distribution 2015 - 2039 \n GFDL SSP1-2.6")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_larvae_gfdl_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
gfdl_temps2 <- readRDS(here('data', 'gfdl_forecast_temp2.rds'))
gfdl_salts2 <- readRDS(here('data', 'gfdl_forecast_salt2.rds'))

df_pklarvae2_gfdl126 <- predict_cells(2040:2069, pk_larvae, 134,
                                   5, 'ssp126', gfdl_temps2, 
                                   gfdl_salts2, pklarvae_formula)
saveRDS(df_pklarvae2_gfdl126, file = here("data", "df_pklarvae2_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pklarvae2_gfdl126, 
             "Forecasted Distribution 2040 - 2069 \n GFDL SSP1-2.6")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_larvae_gfdl_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
gfdl_temps3 <- readRDS(here('data', 'gfdl_forecast_temp3.rds'))
gfdl_salts3 <- readRDS(here('data', 'gfdl_forecast_salt3.rds'))

df_pklarvae3_gfdl126 <- predict_cells(2070:2099, pk_larvae, 134,
                                   5, 'ssp126', gfdl_temps3, 
                                   gfdl_salts3, pklarvae_formula)
saveRDS(df_pklarvae3_gfdl126, file = here("data", "df_pklarvae3_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pklarvae3_gfdl126, 
             "Forecasted Distribution 2070 - 2099 \n GFDL SSP1-2.6")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_larvae_gfdl_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### GFDL 585 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
df_pklarvae1_gfdl585 <- predict_cells(2015:2039, pk_larvae, 134,
                                   5, 'ssp585', gfdl_temps1, 
                                   gfdl_salts1, pklarvae_formula)
saveRDS(df_pklarvae1_gfdl585, file = here("data", "df_pklarvae1_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pklarvae1_gfdl585, 
             "Forecasted Distribution 2015 - 2039 \n GFDL SSP5-8.5")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_larvae_gfdl_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
df_pklarvae2_gfdl585 <- predict_cells(2040:2069, pk_larvae, 134,
                                   5, 'ssp585', gfdl_temps2, 
                                   gfdl_salts2, pklarvae_formula)
saveRDS(df_pklarvae2_gfdl585, file = here("data", "df_pklarvae2_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pklarvae2_gfdl585, 
             "Forecasted Distribution 2040 - 2069 \n GFDL SSP5-8.5")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_larvae_gfdl_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
df_pklarvae3_gfdl585 <- predict_cells(2070:2099, pk_larvae, 134,
                                   5, 'ssp585', gfdl_temps3, 
                                   gfdl_salts3, pklarvae_formula)
saveRDS(df_pklarvae3_gfdl585, file = here("data", "df_pklarvae3_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pklarvae3_gfdl585, 
             "Forecasted Distribution 2070 - 2099 \n GFDL SSP5-8.5")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_larvae_gfdl_ssp585_3.jpg'),
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

df_pklarvae1_miroc126 <- predict_cells(2015:2039, pk_larvae, 134,
                                    5, 'ssp126', miroc_temps1, 
                                    miroc_salts1, pklarvae_formula)
saveRDS(df_pklarvae1_miroc126, file = here("data", "df_pklarvae1_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pklarvae1_miroc126, 
             "Forecasted Distribution 2015 - 2039 \n MIROC SSP1-2.6")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_larvae_miroc_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
miroc_temps2 <- readRDS(here('data', 'miroc_forecast_temp2.rds'))
miroc_salts2 <- readRDS(here('data', 'miroc_forecast_salt2.rds'))

df_pklarvae2_miroc126 <- predict_cells(2040:2069, pk_larvae, 134,
                                    5, 'ssp126', miroc_temps2, 
                                    miroc_salts2, pklarvae_formula)
saveRDS(df_pklarvae2_miroc126, file = here("data", "df_pklarvae2_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pklarvae2_miroc126, 
             "Forecasted Distribution 2040 - 2069 \n MIROC SSP1-2.6")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_larvae_miroc_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
miroc_temps3 <- readRDS(here('data', 'miroc_forecast_temp3.rds'))
miroc_salts3 <- readRDS(here('data', 'miroc_forecast_salt3.rds'))

df_pklarvae3_miroc126 <- predict_cells(2070:2099, pk_larvae, 134,
                                    5, 'ssp126', miroc_temps3, 
                                    miroc_salts3, pklarvae_formula)
saveRDS(df_pklarvae3_miroc126, file = here("data", "df_pklarvae3_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pklarvae3_miroc126, 
             "Forecasted Distribution 2070 - 2099 \n MIROC SSP1-2.6")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_larvae_miroc_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### MIROC 585 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
df_pklarvae1_miroc585 <- predict_cells(2015:2039, pk_larvae, 134,
                                    5, 'ssp585', miroc_temps1, 
                                    miroc_salts1, pklarvae_formula)
saveRDS(df_pklarvae1_miroc585, file = here("data", "df_pklarvae1_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pklarvae1_miroc585, 
             "Forecasted Distribution 2015 - 2039 \n MIROC SSP5-8.5")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_larvae_miroc_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
df_pklarvae2_miroc585 <- predict_cells(2040:2069, pk_larvae, 134,
                                    5, 'ssp585', miroc_temps2, 
                                    miroc_salts2, pklarvae_formula)
saveRDS(df_pklarvae2_miroc585, file = here("data", "df_pklarvae2_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pklarvae2_miroc585, 
             "Forecasted Distribution 2040 - 2069 \n MIROC SSP5-8.5")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_larvae_miroc_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
df_pklarvae3_miroc585 <- predict_cells(2070:2099, pk_larvae, 134,
                                    5, 'ssp585', miroc_temps3,
                                    miroc_salts3, pklarvae_formula)
saveRDS(df_pklarvae3_miroc585, file = here("data", "df_pklarvae3_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pklarvae3_miroc585, 
             "Forecasted Distribution 2070 - 2099 \n MIROC SSP5-8.5")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_larvae_miroc_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

rm(miroc_temps1, miroc_temps2, miroc_temps3,
   miroc_salts1, miroc_salts2, miroc_salts3)

graphics.off()

##### Multi-panel figure -------------------------------------------------------------------------------------------------------------------------
tiff(here('results/pollock_forecast',
          'pollocklarvae_multipanel.tiff'),
     units = "in",
     width = 45,
     height = 90,
     res = 300)
par(mfrow = c(6, 3),
    mar = c(11, 12, 5, 1.5) + 0.1,
    oma = c(3, 25, 15, 1),
    mgp = c(10, 4, 0),
    family = "serif")
grid_multipanel(df_pklarvae1_cesm126)
grid_multipanel(df_pklarvae2_cesm126)
grid_multipanel(df_pklarvae3_cesm126)
grid_multipanel(df_pklarvae1_cesm585)
grid_multipanel(df_pklarvae2_cesm585)
grid_multipanel(df_pklarvae3_cesm585)
grid_multipanel(df_pklarvae1_gfdl126)
grid_multipanel(df_pklarvae2_gfdl126)
grid_multipanel(df_pklarvae3_gfdl126)
grid_multipanel(df_pklarvae1_gfdl585)
grid_multipanel(df_pklarvae2_gfdl585)
grid_multipanel(df_pklarvae3_gfdl585)
grid_multipanel(df_pklarvae1_miroc126)
grid_multipanel(df_pklarvae2_miroc126)
grid_multipanel(df_pklarvae3_miroc126)
grid_multipanel(df_pklarvae1_miroc585)
grid_multipanel(df_pklarvae2_miroc585)
grid_multipanel(df_pklarvae3_miroc585)
mtext("CESM SSP1-2.6", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.92)
mtext("CESM SSP5-8.5", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.76)
mtext("GFDL1-2.6", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.585)
mtext("GFDL5-8.5", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.42)
mtext("MIROC1-2.6", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.26)
mtext("MIROC5-8.5", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.095)
mtext("2015-2039", 
      side = 3, 
      line = 2, 
      outer = TRUE, 
      cex = 7,
      at = 0.17)
mtext("2040-2069", 
      side = 3, 
      line = 2, 
      outer = TRUE, 
      cex = 7,
      at = 0.51)
mtext("2070-2099", 
      side = 3, 
      line = 2, 
      outer = TRUE, 
      cex = 7,
      at = 0.84)
dev.off()

##### Averages ---------------------------------------------------------------------------------------------------------------------------
# 2015 - 2039
df_pklarvae_merged1 <- list(df_pklarvae1_cesm126, df_pklarvae1_cesm585,
                         df_pklarvae1_gfdl126, df_pklarvae1_gfdl585,
                         df_pklarvae1_miroc126, df_pklarvae1_miroc585) 

avg_pklarvae_merged1 <- predict_avgs(df_pklarvae_merged1)

windows(width = 6, height = 6, family = "serif")
grid_predict(avg_pklarvae_merged1, "Forecasted Distribution 2015 - 2039")
dev.copy(jpeg,
         here('results/pollock_forecast/pklarvae_avgs',
              'pollock_larvae_avg1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

# 2040-2069
df_pklarvae_merged2 <- list(df_pklarvae2_cesm126, df_pklarvae2_cesm585,
                         df_pklarvae2_gfdl126, df_pklarvae2_gfdl585,
                         df_pklarvae2_miroc126, df_pklarvae2_miroc585) 

avg_pklarvae_merged2 <- predict_avgs(df_pklarvae_merged2)

windows(width = 6, height = 6, family = "serif")
grid_predict(avg_pklarvae_merged2, "Forecasted Distribution 2040 - 2069")
dev.copy(jpeg,
         here('results/pollock_forecast/pklarvae_avgs',
              'pollock_larvae_avg2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

# 2070-2099
df_pklarvae_merged3 <- list(df_pklarvae3_cesm126, df_pklarvae3_cesm585,
                         df_pklarvae3_gfdl126, df_pklarvae3_gfdl585,
                         df_pklarvae3_miroc126, df_pklarvae3_miroc585) 

avg_pklarvae_merged3 <- predict_avgs(df_pklarvae_merged3)

windows(width = 6, height = 6, family = "serif")
grid_predict(avg_pklarvae_merged3, "Forecasted Distribution 2070 - 2099")
dev.copy(jpeg,
         here('results/pollock_forecast/pklarvae_avgs',
              'pollock_larvae_avg3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

##### GIFs -------------------------------------------------------------------------------------------------------------------------
pklarvae_dir_out <- file.path(base_dir, 'results', 'pollock_forecast', 'pklarvae_avgs')
pklarvae_imgs <- list.files(pklarvae_dir_out, full.names = T)
pklarvae_img_list <- lapply(pklarvae_imgs, image_read)
pklarvae_img_joined <- image_join(pklarvae_img_list)
pklarvae_img_animated <- image_animate(pklarvae_img_joined, fps = 1)
image_write(image = pklarvae_img_animated,
            path = here('results', 'pollock_forecast', "pklarvae_avgs.gif"))

##### Clear environment -------------------------------------------------------------------------------------------------------------------------
rm(df_pklarvae1_cesm126, df_pklarvae1_cesm585,
   df_pklarvae1_gfdl126, df_pklarvae1_gfdl585,
   df_pklarvae1_miroc126, df_pklarvae1_miroc585,
   df_pklarvae2_cesm126, df_pklarvae2_cesm585,
   df_pklarvae2_gfdl126, df_pklarvae2_gfdl585,
   df_pklarvae2_miroc126, df_pklarvae2_miroc585,
   df_pklarvae3_cesm126, df_pklarvae3_cesm585,
   df_pklarvae3_gfdl126, df_pklarvae3_gfdl585,
   df_pklarvae3_miroc126, df_pklarvae3_miroc585,
   avg_pklarvae_merged1, avg_pklarvae_merged2,
   avg_pklarvae_merged3, df_pklarvae_merged1,
   df_pklarvae_merged2, df_pklarvae_merged3,
   pk_larvae, pklarvae_formula, pklarvae_img_animated,
   pklarvae_dir_out, pklarvae_img_joined, 
   pklarvae_img_list, pklarvae_imgs)

### Flathead Eggs --------------------------------------------------------------------------------------------------------------------------
fhs_egg <- load_data('fhs_egg.rds', fhs_egg, roms_temps)
fhsegg_formula <- formula_geog(fhs_egg)

#### Forecast and average into 3 time periods ---------------------------------------------------------------------------------------------
##### CESM 126 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
cesm_temps1 <- readRDS(here('data', 'cesm_forecast_temp1.rds'))
cesm_salts1 <- readRDS(here('data', 'cesm_forecast_salt1.rds'))

df_fhsegg1_cesm126 <- predict_cells(2015:2039, fhs_egg, 154,
                                   6, 'ssp126', cesm_temps1,
                                   cesm_salts1, fhsegg_formula)
saveRDS(df_fhsegg1_cesm126, file = here("data", "df_fhsegg1_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg1_cesm126, 
             "Forecasted Distribution 2015 - 2039 \n CESM SSP1-2.6")
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

df_fhsegg2_cesm126 <- predict_cells(2040:2069, fhs_egg, 154,
                                   6, 'ssp126', cesm_temps2, 
                                   cesm_salts2, fhsegg_formula)
saveRDS(df_fhsegg2_cesm126, file = here("data", "df_fhsegg2_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg2_cesm126, 
             "Forecasted Distribution 2040 - 2069 \n CESM SSP1-2.6")
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

df_fhsegg3_cesm126 <- predict_cells(2070:2099, fhs_egg, 154,
                                   6, 'ssp126', cesm_temps3, 
                                   cesm_salts3, fhsegg_formula)
saveRDS(df_fhsegg3_cesm126, file = here("data", "df_fhsegg3_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg3_cesm126, 
             "Forecasted Distribution 2070 - 2099 \n CESM SSP1-2.6")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_cesm_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### CESM 585 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
df_fhsegg1_cesm585 <- predict_cells(2015:2039, fhs_egg, 154,
                                   6, 'ssp585', cesm_temps1, 
                                   cesm_salts1, fhsegg_formula)
saveRDS(df_fhsegg1_cesm585, file = here("data", "df_fhsegg1_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg1_cesm585, 
             "Forecasted Distribution 2015 - 2039 \n CESM SSP5-8.5")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_cesm_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
df_fhsegg2_cesm585 <- predict_cells(2040:2069, fhs_egg, 154,
                                   6, 'ssp585', cesm_temps2, 
                                   cesm_salts2, fhsegg_formula)
saveRDS(df_fhsegg2_cesm585, file = here("data", "df_fhsegg2_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg2_cesm585, 
             "Forecasted Distribution 2040 - 2069 \n CESM SSP5-8.5")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_cesm_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
df_fhsegg3_cesm585 <- predict_cells(2070:2099, fhs_egg, 154,
                                   6, 'ssp585', cesm_temps3, 
                                   cesm_salts3, fhsegg_formula)
saveRDS(df_fhsegg3_cesm585, file = here("data", "df_fhsegg3_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg3_cesm585, 
             "Forecasted Distribution 2070 - 2099 \n CESM SSP5-8.5")
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
## 2015 - 2039
gfdl_temps1 <- readRDS(here('data', 'gfdl_forecast_temp1.rds'))
gfdl_salts1 <- readRDS(here('data', 'gfdl_forecast_salt1.rds'))

df_fhsegg1_gfdl126 <- predict_cells(2015:2039, fhs_egg, 154,
                                   6, 'ssp126', gfdl_temps1, 
                                   gfdl_salts1, fhsegg_formula)
saveRDS(df_fhsegg1_gfdl126, file = here("data", "df_fhsegg1_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg1_gfdl126, 
             "Forecasted Distribution 2015 - 2039 \n GFDL SSP1-2.6")
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

df_fhsegg2_gfdl126 <- predict_cells(2040:2069, fhs_egg, 154,
                                   6, 'ssp126', gfdl_temps2, 
                                   gfdl_salts2, fhsegg_formula)
saveRDS(df_fhsegg2_gfdl126, file = here("data", "df_fhsegg2_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg2_gfdl126, 
             "Forecasted Distribution 2040 - 2069 \n GFDL SSP1-2.6")
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

df_fhsegg3_gfdl126 <- predict_cells(2070:2099, fhs_egg, 154,
                                   6, 'ssp126', gfdl_temps3, 
                                   gfdl_salts3, fhsegg_formula)
saveRDS(df_fhsegg3_gfdl126, file = here("data", "df_fhsegg3_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg3_gfdl126, 
             "Forecasted Distribution 2070 - 2099 \n GFDL SSP1-2.6")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_gfdl_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### GFDL 585 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
df_fhsegg1_gfdl585 <- predict_cells(2015:2039, fhs_egg, 154,
                                   6, 'ssp585', gfdl_temps1, 
                                   gfdl_salts1, fhsegg_formula)
saveRDS(df_fhsegg1_gfdl585, file = here("data", "df_fhsegg1_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg1_gfdl585, 
             "Forecasted Distribution 2015 - 2039 \n GFDL SSP5-8.5")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_gfdl_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
df_fhsegg2_gfdl585 <- predict_cells(2040:2069, fhs_egg, 154,
                                   6, 'ssp585', gfdl_temps2, 
                                   gfdl_salts2, fhsegg_formula)
saveRDS(df_fhsegg2_gfdl585, file = here("data", "df_fhsegg2_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg2_gfdl585, 
             "Forecasted Distribution 2040 - 2069 \n GFDL SSP5-8.5")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_gfdl_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
df_fhsegg3_gfdl585 <- predict_cells(2070:2099, fhs_egg, 154,
                                   6, 'ssp585', gfdl_temps3, 
                                   gfdl_salts3, fhsegg_formula)
saveRDS(df_fhsegg3_gfdl585, file = here("data", "df_fhsegg3_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg3_gfdl585, 
             "Forecasted Distribution 2070 - 2099 \n GFDL SSP5-8.5")
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
## 2015 - 2039
miroc_temps1 <- readRDS(here('data', 'miroc_forecast_temp1.rds'))
miroc_salts1 <- readRDS(here('data', 'miroc_forecast_salt1.rds'))

df_fhsegg1_miroc126 <- predict_cells(2015:2039, fhs_egg, 154,
                                    6, 'ssp126', miroc_temps1, 
                                    miroc_salts1, fhsegg_formula)
saveRDS(df_fhsegg1_miroc126, file = here("data", "df_fhsegg1_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg1_miroc126, 
             "Forecasted Distribution 2015 - 2039 \n MIROC SSP1-2.6")
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

df_fhsegg2_miroc126 <- predict_cells(2040:2069, fhs_egg, 154,
                                    6, 'ssp126', miroc_temps2, 
                                    miroc_salts2, fhsegg_formula)
saveRDS(df_fhsegg2_miroc126, file = here("data", "df_fhsegg2_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg2_miroc126, 
             "Forecasted Distribution 2040 - 2069 \n MIROC SSP1-2.6")
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

df_fhsegg3_miroc126 <- predict_cells(2070:2099, fhs_egg, 154,
                                    6, 'ssp126', miroc_temps3, 
                                    miroc_salts3, fhsegg_formula)
saveRDS(df_fhsegg3_miroc126, file = here("data", "df_fhsegg3_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg3_miroc126, 
             "Forecasted Distribution 2070 - 2099 \n MIROC SSP1-2.6")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_miroc_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### MIROC 585 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
df_fhsegg1_miroc585 <- predict_cells(2015:2039, fhs_egg, 154,
                                    6, 'ssp585', miroc_temps1, 
                                    miroc_salts1, fhsegg_formula)
saveRDS(df_fhsegg1_miroc585, file = here("data", "df_fhsegg1_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg1_miroc585, 
             "Forecasted Distribution 2015 - 2039 \n MIROC SSP5-8.5")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_miroc_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
df_fhsegg2_miroc585 <- predict_cells(2040:2069, fhs_egg, 154,
                                    6, 'ssp585', miroc_temps2, 
                                    miroc_salts2, fhsegg_formula)
saveRDS(df_fhsegg2_miroc585, file = here("data", "df_fhsegg2_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg2_miroc585, 
             "Forecasted Distribution 2040 - 2069 \n MIROC SSP5-8.5")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_miroc_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
df_fhsegg3_miroc585 <- predict_cells(2070:2099, fhs_egg, 154,
                                    6, 'ssp585', miroc_temps3,
                                    miroc_salts3, fhsegg_formula)
saveRDS(df_fhsegg3_miroc585, file = here("data", "df_fhsegg3_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhsegg3_miroc585, 
             "Forecasted Distribution 2070 - 2099 \n MIROC SSP5-8.5")
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

graphics.off()

##### Multi-panel figure -------------------------------------------------------------------------------------------------------------------------
tiff(here('results/flathead_forecast',
          'flatheadegg_multipanel.tiff'),
     units = "in",
     width = 45,
     height = 90,
     res = 300)
par(mfrow = c(6, 3),
    mar = c(11, 12, 5, 1.5) + 0.1,
    oma = c(3, 25, 15, 1),
    mgp = c(10, 4, 0),
    family = "serif")
grid_multipanel(df_fhsegg1_cesm126)
grid_multipanel(df_fhsegg2_cesm126)
grid_multipanel(df_fhsegg3_cesm126)
grid_multipanel(df_fhsegg1_cesm585)
grid_multipanel(df_fhsegg2_cesm585)
grid_multipanel(df_fhsegg3_cesm585)
grid_multipanel(df_fhsegg1_gfdl126)
grid_multipanel(df_fhsegg2_gfdl126)
grid_multipanel(df_fhsegg3_gfdl126)
grid_multipanel(df_fhsegg1_gfdl585)
grid_multipanel(df_fhsegg2_gfdl585)
grid_multipanel(df_fhsegg3_gfdl585)
grid_multipanel(df_fhsegg1_miroc126)
grid_multipanel(df_fhsegg2_miroc126)
grid_multipanel(df_fhsegg3_miroc126)
grid_multipanel(df_fhsegg1_miroc585)
grid_multipanel(df_fhsegg2_miroc585)
grid_multipanel(df_fhsegg3_miroc585)
mtext("CESM SSP1-2.6", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.92)
mtext("CESM SSP5-8.5", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.76)
mtext("GFDL1-2.6", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.585)
mtext("GFDL5-8.5", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.42)
mtext("MIROC1-2.6", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.26)
mtext("MIROC5-8.5", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.095)
mtext("2015-2039", 
      side = 3, 
      line = 2, 
      outer = TRUE, 
      cex = 7,
      at = 0.17)
mtext("2040-2069", 
      side = 3, 
      line = 2, 
      outer = TRUE, 
      cex = 7,
      at = 0.51)
mtext("2070-2099", 
      side = 3, 
      line = 2, 
      outer = TRUE, 
      cex = 7,
      at = 0.84)
dev.off()

##### Averages ---------------------------------------------------------------------------------------------------------------------------
# 2015 - 2039
df_fhsegg_merged1 <- list(df_fhsegg1_cesm126, df_fhsegg1_cesm585,
                         df_fhsegg1_gfdl126, df_fhsegg1_gfdl585,
                         df_fhsegg1_miroc126, df_fhsegg1_miroc585) 

avg_fhsegg_merged1 <- predict_avgs(df_fhsegg_merged1)

windows(width = 6, height = 6, family = "serif")
grid_predict(avg_fhsegg_merged1, "Forecasted Distribution 2015 - 2039")
dev.copy(jpeg,
         here('results/flathead_forecast/fhsegg_avgs',
              'flathead_egg_avg1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

# 2040-2069
df_fhsegg_merged2 <- list(df_fhsegg2_cesm126, df_fhsegg2_cesm585,
                         df_fhsegg2_gfdl126, df_fhsegg2_gfdl585,
                         df_fhsegg2_miroc126, df_fhsegg2_miroc585) 

avg_fhsegg_merged2 <- predict_avgs(df_fhsegg_merged2)

windows(width = 6, height = 6, family = "serif")
grid_predict(avg_fhsegg_merged2, "Forecasted Distribution 2040 - 2069")
dev.copy(jpeg,
         here('results/flathead_forecast/fhsegg_avgs',
              'flathead_egg_avg2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

# 2070-2099
df_fhsegg_merged3 <- list(df_fhsegg3_cesm126, df_fhsegg3_cesm585,
                         df_fhsegg3_gfdl126, df_fhsegg3_gfdl585,
                         df_fhsegg3_miroc126, df_fhsegg3_miroc585) 

avg_fhsegg_merged3 <- predict_avgs(df_fhsegg_merged3)

windows(width = 6, height = 6, family = "serif")
grid_predict(avg_fhsegg_merged3, "Forecasted Distribution 2070 - 2099")
dev.copy(jpeg,
         here('results/flathead_forecast/fhsegg_avgs',
              'flathead_egg_avg3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

##### GIFs -------------------------------------------------------------------------------------------------------------------------
fhsegg_dir_out <- file.path(base_dir, 'results', 'flathead_forecast', 'fhsegg_avgs')
fhsegg_imgs <- list.files(fhsegg_dir_out, full.names = T)
fhsegg_img_list <- lapply(fhsegg_imgs, image_read)
fhsegg_img_joined <- image_join(fhsegg_img_list)
fhsegg_img_animated <- image_animate(fhsegg_img_joined, fps = 1)
image_write(image = fhsegg_img_animated,
            path = here('results', 'flathead_forecast', "fhsegg_avgs.gif"))

##### Clear environment -------------------------------------------------------------------------------------------------------------------------
rm(df_fhsegg1_cesm126, df_fhsegg1_cesm585,
   df_fhsegg1_gfdl126, df_fhsegg1_gfdl585,
   df_fhsegg1_miroc126, df_fhsegg1_miroc585,
   df_fhsegg2_cesm126, df_fhsegg2_cesm585,
   df_fhsegg2_gfdl126, df_fhsegg2_gfdl585,
   df_fhsegg2_miroc126, df_fhsegg2_miroc585,
   df_fhsegg3_cesm126, df_fhsegg3_cesm585,
   df_fhsegg3_gfdl126, df_fhsegg3_gfdl585,
   df_fhsegg3_miroc126, df_fhsegg3_miroc585,
   avg_fhsegg_merged1, avg_fhsegg_merged2,
   avg_fhsegg_merged3, df_fhsegg_merged1,
   df_fhsegg_merged2, df_fhsegg_merged3,
   fhs_egg, fhsegg_formula, fhsegg_img_animated,
   fhsegg_dir_out, fhsegg_img_joined, 
   fhsegg_img_list, fhsegg_imgs)


### Flathead Larvae --------------------------------------------------------------------------------------------------------------------------
fhs_larvae <- load_data('fhs_larvae.rds', fhs_larvae, roms_temps)
fhslarvae_formula <- formula_geog(fhs_larvae)

#### Forecast and average into 3 time periods ---------------------------------------------------------------------------------------------
##### CESM 126 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
cesm_temps1 <- readRDS(here('data', 'cesm_forecast_temp1.rds'))
cesm_salts1 <- readRDS(here('data', 'cesm_forecast_salt1.rds'))

df_fhslarvae1_cesm126 <- predict_cells(2015:2039, fhs_larvae, 138,
                                      5, 'ssp126', cesm_temps1,
                                      cesm_salts1, fhslarvae_formula)
saveRDS(df_fhslarvae1_cesm126, file = here("data", "df_fhslarvae1_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae1_cesm126, 
             "Forecasted Distribution 2015 - 2039 \n CESM SSP1-2.6")
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

df_fhslarvae2_cesm126 <- predict_cells(2040:2069, fhs_larvae, 138,
                                      5, 'ssp126', cesm_temps2, 
                                      cesm_salts2, fhslarvae_formula)
saveRDS(df_fhslarvae2_cesm126, file = here("data", "df_fhslarvae2_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae2_cesm126, 
             "Forecasted Distribution 2040 - 2069 \n CESM SSP1-2.6")
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

df_fhslarvae3_cesm126 <- predict_cells(2070:2099, fhs_larvae, 138,
                                      5, 'ssp126', cesm_temps3, 
                                      cesm_salts3, fhslarvae_formula)
saveRDS(df_fhslarvae3_cesm126, file = here("data", "df_fhslarvae3_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae3_cesm126, 
             "Forecasted Distribution 2070 - 2099 \n CESM SSP1-2.6")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_cesm_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### CESM 585 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
df_fhslarvae1_cesm585 <- predict_cells(2015:2039, fhs_larvae, 138,
                                      5, 'ssp585', cesm_temps1, 
                                      cesm_salts1, fhslarvae_formula)
saveRDS(df_fhslarvae1_cesm585, file = here("data", "df_fhslarvae1_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae1_cesm585, 
             "Forecasted Distribution 2015 - 2039 \n CESM SSP5-8.5")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_cesm_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
df_fhslarvae2_cesm585 <- predict_cells(2040:2069, fhs_larvae, 138,
                                      5, 'ssp585', cesm_temps2, 
                                      cesm_salts2, fhslarvae_formula)
saveRDS(df_fhslarvae2_cesm585, file = here("data", "df_fhslarvae2_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae2_cesm585, 
             "Forecasted Distribution 2040 - 2069 \n CESM SSP5-8.5")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_cesm_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
df_fhslarvae3_cesm585 <- predict_cells(2070:2099, fhs_larvae, 138,
                                      5, 'ssp585', cesm_temps3, 
                                      cesm_salts3, fhslarvae_formula)
saveRDS(df_fhslarvae3_cesm585, file = here("data", "df_fhslarvae3_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae3_cesm585, 
             "Forecasted Distribution 2070 - 2099 \n CESM SSP5-8.5")
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
## 2015 - 2039
gfdl_temps1 <- readRDS(here('data', 'gfdl_forecast_temp1.rds'))
gfdl_salts1 <- readRDS(here('data', 'gfdl_forecast_salt1.rds'))

df_fhslarvae1_gfdl126 <- predict_cells(2015:2039, fhs_larvae, 138,
                                      5, 'ssp126', gfdl_temps1, 
                                      gfdl_salts1, fhslarvae_formula)
saveRDS(df_fhslarvae1_gfdl126, file = here("data", "df_fhslarvae1_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae1_gfdl126, 
             "Forecasted Distribution 2015 - 2039 \n GFDL SSP1-2.6")
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

df_fhslarvae2_gfdl126 <- predict_cells(2040:2069, fhs_larvae, 138,
                                      5, 'ssp126', gfdl_temps2, 
                                      gfdl_salts2, fhslarvae_formula)
saveRDS(df_fhslarvae2_gfdl126, file = here("data", "df_fhslarvae2_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae2_gfdl126, 
             "Forecasted Distribution 2040 - 2069 \n GFDL SSP1-2.6")
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

df_fhslarvae3_gfdl126 <- predict_cells(2070:2099, fhs_larvae, 138,
                                      5, 'ssp126', gfdl_temps3, 
                                      gfdl_salts3, fhslarvae_formula)
saveRDS(df_fhslarvae3_gfdl126, file = here("data", "df_fhslarvae3_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae3_gfdl126, 
             "Forecasted Distribution 2070 - 2099 \n GFDL SSP1-2.6")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_gfdl_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### GFDL 585 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
df_fhslarvae1_gfdl585 <- predict_cells(2015:2039, fhs_larvae, 138,
                                      5, 'ssp585', gfdl_temps1, 
                                      gfdl_salts1, fhslarvae_formula)
saveRDS(df_fhslarvae1_gfdl585, file = here("data", "df_fhslarvae1_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae1_gfdl585, 
             "Forecasted Distribution 2015 - 2039 \n GFDL SSP5-8.5")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_gfdl_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
df_fhslarvae2_gfdl585 <- predict_cells(2040:2069, fhs_larvae, 138,
                                      5, 'ssp585', gfdl_temps2, 
                                      gfdl_salts2, fhslarvae_formula)
saveRDS(df_fhslarvae2_gfdl585, file = here("data", "df_fhslarvae2_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae2_gfdl585, 
             "Forecasted Distribution 2040 - 2069 \n GFDL SSP5-8.5")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_gfdl_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
df_fhslarvae3_gfdl585 <- predict_cells(2070:2099, fhs_larvae, 138,
                                      5, 'ssp585', gfdl_temps3, 
                                      gfdl_salts3, fhslarvae_formula)
saveRDS(df_fhslarvae3_gfdl585, file = here("data", "df_fhslarvae3_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae3_gfdl585, 
             "Forecasted Distribution 2070 - 2099 \n GFDL SSP5-8.5")
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
## 2015 - 2039
miroc_temps1 <- readRDS(here('data', 'miroc_forecast_temp1.rds'))
miroc_salts1 <- readRDS(here('data', 'miroc_forecast_salt1.rds'))

df_fhslarvae1_miroc126 <- predict_cells(2015:2039, fhs_larvae, 138,
                                       5, 'ssp126', miroc_temps1, 
                                       miroc_salts1, fhslarvae_formula)
saveRDS(df_fhslarvae1_miroc126, file = here("data", "df_fhslarvae1_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae1_miroc126, 
             "Forecasted Distribution 2015 - 2039 \n MIROC SSP1-2.6")
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

df_fhslarvae2_miroc126 <- predict_cells(2040:2069, fhs_larvae, 138,
                                       5, 'ssp126', miroc_temps2, 
                                       miroc_salts2, fhslarvae_formula)
saveRDS(df_fhslarvae2_miroc126, file = here("data", "df_fhslarvae2_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae2_miroc126, 
             "Forecasted Distribution 2040 - 2069 \n MIROC SSP1-2.6")
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

df_fhslarvae3_miroc126 <- predict_cells(2070:2099, fhs_larvae, 138,
                                       5, 'ssp126', miroc_temps3, 
                                       miroc_salts3, fhslarvae_formula)
saveRDS(df_fhslarvae3_miroc126, file = here("data", "df_fhslarvae3_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae3_miroc126, 
             "Forecasted Distribution 2070 - 2099 \n MIROC SSP1-2.6")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_miroc_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### MIROC 585 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
df_fhslarvae1_miroc585 <- predict_cells(2015:2039, fhs_larvae, 138,
                                       5, 'ssp585', miroc_temps1, 
                                       miroc_salts1, fhslarvae_formula)
saveRDS(df_fhslarvae1_miroc585, file = here("data", "df_fhslarvae1_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae1_miroc585, 
             "Forecasted Distribution 2015 - 2039 \n MIROC SSP5-8.5")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_miroc_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
df_fhslarvae2_miroc585 <- predict_cells(2040:2069, fhs_larvae, 138,
                                       5, 'ssp585', miroc_temps2, 
                                       miroc_salts2, fhslarvae_formula)
saveRDS(df_fhslarvae2_miroc585, file = here("data", "df_fhslarvae2_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae2_miroc585, 
             "Forecasted Distribution 2040 - 2069 \n MIROC SSP5-8.5")
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_miroc_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
df_fhslarvae3_miroc585 <- predict_cells(2070:2099, fhs_larvae, 138,
                                       5, 'ssp585', miroc_temps3,
                                       miroc_salts3, fhslarvae_formula)
saveRDS(df_fhslarvae3_miroc585, file = here("data", "df_fhslarvae3_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_fhslarvae3_miroc585, 
             "Forecasted Distribution 2070 - 2099 \n MIROC SSP5-8.5")
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

graphics.off()

##### Multi-panel figure -------------------------------------------------------------------------------------------------------------------------
tiff(here('results/flathead_forecast',
          'flatheadlarvae_multipanel.tiff'),
     units = "in",
     width = 45,
     height = 90,
     res = 300)
par(mfrow = c(6, 3),
    mar = c(11, 12, 5, 1.5) + 0.1,
    oma = c(3, 25, 15, 1),
    mgp = c(10, 4, 0),
    family = "serif")
grid_multipanel(df_fhslarvae1_cesm126)
grid_multipanel(df_fhslarvae2_cesm126)
grid_multipanel(df_fhslarvae3_cesm126)
grid_multipanel(df_fhslarvae1_cesm585)
grid_multipanel(df_fhslarvae2_cesm585)
grid_multipanel(df_fhslarvae3_cesm585)
grid_multipanel(df_fhslarvae1_gfdl126)
grid_multipanel(df_fhslarvae2_gfdl126)
grid_multipanel(df_fhslarvae3_gfdl126)
grid_multipanel(df_fhslarvae1_gfdl585)
grid_multipanel(df_fhslarvae2_gfdl585)
grid_multipanel(df_fhslarvae3_gfdl585)
grid_multipanel(df_fhslarvae1_miroc126)
grid_multipanel(df_fhslarvae2_miroc126)
grid_multipanel(df_fhslarvae3_miroc126)
grid_multipanel(df_fhslarvae1_miroc585)
grid_multipanel(df_fhslarvae2_miroc585)
grid_multipanel(df_fhslarvae3_miroc585)
mtext("CESM SSP1-2.6", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.92)
mtext("CESM SSP5-8.5", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.76)
mtext("GFDL1-2.6", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.585)
mtext("GFDL5-8.5", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.42)
mtext("MIROC1-2.6", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.26)
mtext("MIROC5-8.5", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.095)
mtext("2015-2039", 
      side = 3, 
      line = 2, 
      outer = TRUE, 
      cex = 7,
      at = 0.17)
mtext("2040-2069", 
      side = 3, 
      line = 2, 
      outer = TRUE, 
      cex = 7,
      at = 0.51)
mtext("2070-2099", 
      side = 3, 
      line = 2, 
      outer = TRUE, 
      cex = 7,
      at = 0.84)
dev.off()

##### Averages ---------------------------------------------------------------------------------------------------------------------------
# 2015 - 2039
df_fhslarvae_merged1 <- list(df_fhslarvae1_cesm126, df_fhslarvae1_cesm585,
                            df_fhslarvae1_gfdl126, df_fhslarvae1_gfdl585,
                            df_fhslarvae1_miroc126, df_fhslarvae1_miroc585) 

avg_fhslarvae_merged1 <- predict_avgs(df_fhslarvae_merged1)

windows(width = 6, height = 6, family = "serif")
grid_predict(avg_fhslarvae_merged1, "Forecasted Distribution 2015 - 2039")
dev.copy(jpeg,
         here('results/flathead_forecast/fhslarvae_avgs',
              'flathead_larvae_avg1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

# 2040-2069
df_fhslarvae_merged2 <- list(df_fhslarvae2_cesm126, df_fhslarvae2_cesm585,
                            df_fhslarvae2_gfdl126, df_fhslarvae2_gfdl585,
                            df_fhslarvae2_miroc126, df_fhslarvae2_miroc585) 

avg_fhslarvae_merged2 <- predict_avgs(df_fhslarvae_merged2)

windows(width = 6, height = 6, family = "serif")
grid_predict(avg_fhslarvae_merged2, "Forecasted Distribution 2040 - 2069")
dev.copy(jpeg,
         here('results/flathead_forecast/fhslarvae_avgs',
              'flathead_larvae_avg2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

# 2070-2099
df_fhslarvae_merged3 <- list(df_fhslarvae3_cesm126, df_fhslarvae3_cesm585,
                            df_fhslarvae3_gfdl126, df_fhslarvae3_gfdl585,
                            df_fhslarvae3_miroc126, df_fhslarvae3_miroc585) 

avg_fhslarvae_merged3 <- predict_avgs(df_fhslarvae_merged3)

windows(width = 6, height = 6, family = "serif")
grid_predict(avg_fhslarvae_merged3, "Forecasted Distribution 2070 - 2099")
dev.copy(jpeg,
         here('results/flathead_forecast/fhslarvae_avgs',
              'flathead_larvae_avg3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

##### GIFs -------------------------------------------------------------------------------------------------------------------------
fhslarvae_dir_out <- file.path(base_dir, 'results', 'flathead_forecast', 'fhslarvae_avgs')
fhslarvae_imgs <- list.files(fhslarvae_dir_out, full.names = T)
fhslarvae_img_list <- lapply(fhslarvae_imgs, image_read)
fhslarvae_img_joined <- image_join(fhslarvae_img_list)
fhslarvae_img_animated <- image_animate(fhslarvae_img_joined, fps = 1)
image_write(image = fhslarvae_img_animated,
            path = here('results', 'flathead_forecast', "fhslarvae_avgs.gif"))

##### Clear environment -------------------------------------------------------------------------------------------------------------------------
rm(df_fhslarvae1_cesm126, df_fhslarvae1_cesm585,
   df_fhslarvae1_gfdl126, df_fhslarvae1_gfdl585,
   df_fhslarvae1_miroc126, df_fhslarvae1_miroc585,
   df_fhslarvae2_cesm126, df_fhslarvae2_cesm585,
   df_fhslarvae2_gfdl126, df_fhslarvae2_gfdl585,
   df_fhslarvae2_miroc126, df_fhslarvae2_miroc585,
   df_fhslarvae3_cesm126, df_fhslarvae3_cesm585,
   df_fhslarvae3_gfdl126, df_fhslarvae3_gfdl585,
   df_fhslarvae3_miroc126, df_fhslarvae3_miroc585,
   avg_fhslarvae_merged1, avg_fhslarvae_merged2,
   avg_fhslarvae_merged3, df_fhslarvae_merged1,
   df_fhslarvae_merged2, df_fhslarvae_merged3,
   fhs_larvae, fhslarvae_formula, fhslarvae_img_animated,
   fhslarvae_dir_out, fhslarvae_img_joined, 
   fhslarvae_img_list, fhslarvae_imgs)

### Alaska Plaice Eggs --------------------------------------------------------------------------------------------------------------------------
akp_egg <- load_data('akp_egg.rds', akp_egg, roms_temps)
akpegg_formula <- formula_geog(akp_egg)

#### Forecast and average into 3 time periods ---------------------------------------------------------------------------------------------
##### CESM 126 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
cesm_temps1 <- readRDS(here('data', 'cesm_forecast_temp1.rds'))
cesm_salts1 <- readRDS(here('data', 'cesm_forecast_salt1.rds'))

df_akpegg1_cesm126 <- predict_cells(2015:2039, akp_egg, 134,
                                   5, 'ssp126', cesm_temps1,
                                   cesm_salts1, akpegg_formula)
saveRDS(df_akpegg1_cesm126, file = here("data", "df_akpegg1_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akpegg1_cesm126, 
             "Forecasted Distribution 2015 - 2039 \n CESM SSP1-2.6")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_cesm_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
cesm_temps2 <- readRDS(here('data', 'cesm_forecast_temp2.rds'))
cesm_salts2 <- readRDS(here('data', 'cesm_forecast_salt2.rds'))

df_akpegg2_cesm126 <- predict_cells(2040:2069, akp_egg, 134,
                                   5, 'ssp126', cesm_temps2, 
                                   cesm_salts2, akpegg_formula)
saveRDS(df_akpegg2_cesm126, file = here("data", "df_akpegg2_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akpegg2_cesm126, 
             "Forecasted Distribution 2040 - 2069 \n CESM SSP1-2.6")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_cesm_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
cesm_temps3 <- readRDS(here('data', 'cesm_forecast_temp3.rds'))
cesm_salts3 <- readRDS(here('data', 'cesm_forecast_salt3.rds'))

df_akpegg3_cesm126 <- predict_cells(2070:2099, akp_egg, 134,
                                   5, 'ssp126', cesm_temps3, 
                                   cesm_salts3, akpegg_formula)
saveRDS(df_akpegg3_cesm126, file = here("data", "df_akpegg3_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akpegg3_cesm126, 
             "Forecasted Distribution 2070 - 2099 \n CESM SSP1-2.6")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_cesm_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### CESM 585 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
df_akpegg1_cesm585 <- predict_cells(2015:2039, akp_egg, 134,
                                   5, 'ssp585', cesm_temps1, 
                                   cesm_salts1, akpegg_formula)
saveRDS(df_akpegg1_cesm585, file = here("data", "df_akpegg1_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akpegg1_cesm585, 
             "Forecasted Distribution 2015 - 2039 \n CESM SSP5-8.5")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_cesm_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
df_akpegg2_cesm585 <- predict_cells(2040:2069, akp_egg, 134,
                                   5, 'ssp585', cesm_temps2, 
                                   cesm_salts2, akpegg_formula)
saveRDS(df_akpegg2_cesm585, file = here("data", "df_akpegg2_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akpegg2_cesm585, 
             "Forecasted Distribution 2040 - 2069 \n CESM SSP5-8.5")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_cesm_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
df_akpegg3_cesm585 <- predict_cells(2070:2099, akp_egg, 134,
                                   5, 'ssp585', cesm_temps3, 
                                   cesm_salts3, akpegg_formula)
saveRDS(df_akpegg3_cesm585, file = here("data", "df_akpegg3_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akpegg3_cesm585, 
             "Forecasted Distribution 2070 - 2099 \n CESM SSP5-8.5")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_cesm_ssp585_3.jpg'),
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

df_akpegg1_gfdl126 <- predict_cells(2015:2039, akp_egg, 134,
                                   5, 'ssp126', gfdl_temps1, 
                                   gfdl_salts1, akpegg_formula)
saveRDS(df_akpegg1_gfdl126, file = here("data", "df_akpegg1_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akpegg1_gfdl126, 
             "Forecasted Distribution 2015 - 2039 \n GFDL SSP1-2.6")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_gfdl_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
gfdl_temps2 <- readRDS(here('data', 'gfdl_forecast_temp2.rds'))
gfdl_salts2 <- readRDS(here('data', 'gfdl_forecast_salt2.rds'))

df_akpegg2_gfdl126 <- predict_cells(2040:2069, akp_egg, 134,
                                   5, 'ssp126', gfdl_temps2, 
                                   gfdl_salts2, akpegg_formula)
saveRDS(df_akpegg2_gfdl126, file = here("data", "df_akpegg2_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akpegg2_gfdl126, 
             "Forecasted Distribution 2040 - 2069 \n GFDL SSP1-2.6")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_gfdl_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
gfdl_temps3 <- readRDS(here('data', 'gfdl_forecast_temp3.rds'))
gfdl_salts3 <- readRDS(here('data', 'gfdl_forecast_salt3.rds'))

df_akpegg3_gfdl126 <- predict_cells(2070:2099, akp_egg, 134,
                                   5, 'ssp126', gfdl_temps3, 
                                   gfdl_salts3, akpegg_formula)
saveRDS(df_akpegg3_gfdl126, file = here("data", "df_akpegg3_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akpegg3_gfdl126, 
             "Forecasted Distribution 2070 - 2099 \n GFDL SSP1-2.6")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_gfdl_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### GFDL 585 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
df_akpegg1_gfdl585 <- predict_cells(2015:2039, akp_egg, 134,
                                   5, 'ssp585', gfdl_temps1, 
                                   gfdl_salts1, akpegg_formula)
saveRDS(df_akpegg1_gfdl585, file = here("data", "df_akpegg1_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akpegg1_gfdl585, 
             "Forecasted Distribution 2015 - 2039 \n GFDL SSP5-8.5")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_gfdl_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
df_akpegg2_gfdl585 <- predict_cells(2040:2069, akp_egg, 134,
                                   5, 'ssp585', gfdl_temps2, 
                                   gfdl_salts2, akpegg_formula)
saveRDS(df_akpegg2_gfdl585, file = here("data", "df_akpegg2_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akpegg2_gfdl585, 
             "Forecasted Distribution 2040 - 2069 \n GFDL SSP5-8.5")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_gfdl_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
df_akpegg3_gfdl585 <- predict_cells(2070:2099, akp_egg, 134,
                                   5, 'ssp585', gfdl_temps3, 
                                   gfdl_salts3, akpegg_formula)
saveRDS(df_akpegg3_gfdl585, file = here("data", "df_akpegg3_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akpegg3_gfdl585, 
             "Forecasted Distribution 2070 - 2099 \n GFDL SSP5-8.5")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_gfdl_ssp585_3.jpg'),
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

df_akpegg1_miroc126 <- predict_cells(2015:2039, akp_egg, 134,
                                    5, 'ssp126', miroc_temps1, 
                                    miroc_salts1, akpegg_formula)
saveRDS(df_akpegg1_miroc126, file = here("data", "df_akpegg1_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akpegg1_miroc126, 
             "Forecasted Distribution 2015 - 2039 \n MIROC SSP1-2.6")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_miroc_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
miroc_temps2 <- readRDS(here('data', 'miroc_forecast_temp2.rds'))
miroc_salts2 <- readRDS(here('data', 'miroc_forecast_salt2.rds'))

df_akpegg2_miroc126 <- predict_cells(2040:2069, akp_egg, 134,
                                    5, 'ssp126', miroc_temps2, 
                                    miroc_salts2, akpegg_formula)
saveRDS(df_akpegg2_miroc126, file = here("data", "df_akpegg2_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akpegg2_miroc126, 
             "Forecasted Distribution 2040 - 2069 \n MIROC SSP1-2.6")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_miroc_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
miroc_temps3 <- readRDS(here('data', 'miroc_forecast_temp3.rds'))
miroc_salts3 <- readRDS(here('data', 'miroc_forecast_salt3.rds'))

df_akpegg3_miroc126 <- predict_cells(2070:2099, akp_egg, 134,
                                    5, 'ssp126', miroc_temps3, 
                                    miroc_salts3, akpegg_formula)
saveRDS(df_akpegg3_miroc126, file = here("data", "df_akpegg3_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akpegg3_miroc126, 
             "Forecasted Distribution 2070 - 2099 \n MIROC SSP1-2.6")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_miroc_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### MIROC 585 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
df_akpegg1_miroc585 <- predict_cells(2015:2039, akp_egg, 134,
                                    5, 'ssp585', miroc_temps1, 
                                    miroc_salts1, akpegg_formula)
saveRDS(df_akpegg1_miroc585, file = here("data", "df_akpegg1_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akpegg1_miroc585, 
             "Forecasted Distribution 2015 - 2039 \n MIROC SSP5-8.5")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_miroc_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
df_akpegg2_miroc585 <- predict_cells(2040:2069, akp_egg, 134,
                                    5, 'ssp585', miroc_temps2, 
                                    miroc_salts2, akpegg_formula)
saveRDS(df_akpegg2_miroc585, file = here("data", "df_akpegg2_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akpegg2_miroc585, 
             "Forecasted Distribution 2040 - 2069 \n MIROC SSP5-8.5")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_miroc_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
df_akpegg3_miroc585 <- predict_cells(2070:2099, akp_egg, 134,
                                    5, 'ssp585', miroc_temps3,
                                    miroc_salts3, akpegg_formula)
saveRDS(df_akpegg3_miroc585, file = here("data", "df_akpegg3_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akpegg3_miroc585, 
             "Forecasted Distribution 2070 - 2099 \n MIROC SSP5-8.5")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_miroc_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

rm(miroc_temps1, miroc_temps2, miroc_temps3,
   miroc_salts1, miroc_salts2, miroc_salts3)

graphics.off()

##### Multi-panel figure -------------------------------------------------------------------------------------------------------------------------
tiff(here('results/plaice_forecast',
          'plaiceegg_multipanel.tiff'),
     units = "in",
     width = 45,
     height = 90,
     res = 300)
par(mfrow = c(6, 3),
    mar = c(11, 12, 5, 1.5) + 0.1,
    oma = c(3, 25, 15, 1),
    mgp = c(10, 4, 0),
    family = "serif")
grid_multipanel(df_akpegg1_cesm126)
grid_multipanel(df_akpegg2_cesm126)
grid_multipanel(df_akpegg3_cesm126)
grid_multipanel(df_akpegg1_cesm585)
grid_multipanel(df_akpegg2_cesm585)
grid_multipanel(df_akpegg3_cesm585)
grid_multipanel(df_akpegg1_gfdl126)
grid_multipanel(df_akpegg2_gfdl126)
grid_multipanel(df_akpegg3_gfdl126)
grid_multipanel(df_akpegg1_gfdl585)
grid_multipanel(df_akpegg2_gfdl585)
grid_multipanel(df_akpegg3_gfdl585)
grid_multipanel(df_akpegg1_miroc126)
grid_multipanel(df_akpegg2_miroc126)
grid_multipanel(df_akpegg3_miroc126)
grid_multipanel(df_akpegg1_miroc585)
grid_multipanel(df_akpegg2_miroc585)
grid_multipanel(df_akpegg3_miroc585)
mtext("CESM SSP1-2.6", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.92)
mtext("CESM SSP5-8.5", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.76)
mtext("GFDL1-2.6", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.585)
mtext("GFDL5-8.5", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.42)
mtext("MIROC1-2.6", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.26)
mtext("MIROC5-8.5", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.095)
mtext("2015-2039", 
      side = 3, 
      line = 2, 
      outer = TRUE, 
      cex = 7,
      at = 0.17)
mtext("2040-2069", 
      side = 3, 
      line = 2, 
      outer = TRUE, 
      cex = 7,
      at = 0.51)
mtext("2070-2099", 
      side = 3, 
      line = 2, 
      outer = TRUE, 
      cex = 7,
      at = 0.84)
dev.off()

##### Averages ---------------------------------------------------------------------------------------------------------------------------
# 2015 - 2039
df_akpegg_merged1 <- list(df_akpegg1_cesm126, df_akpegg1_cesm585,
                         df_akpegg1_gfdl126, df_akpegg1_gfdl585,
                         df_akpegg1_miroc126, df_akpegg1_miroc585) 

avg_akpegg_merged1 <- predict_avgs(df_akpegg_merged1)

windows(width = 6, height = 6, family = "serif")
grid_predict(avg_akpegg_merged1, "Forecasted Distribution 2015 - 2039")
dev.copy(jpeg,
         here('results/plaice_forecast/akpegg_avgs',
              'plaice_egg_avg1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

# 2040-2069
df_akpegg_merged2 <- list(df_akpegg2_cesm126, df_akpegg2_cesm585,
                         df_akpegg2_gfdl126, df_akpegg2_gfdl585,
                         df_akpegg2_miroc126, df_akpegg2_miroc585) 

avg_akpegg_merged2 <- predict_avgs(df_akpegg_merged2)

windows(width = 6, height = 6, family = "serif")
grid_predict(avg_akpegg_merged2, "Forecasted Distribution 2040 - 2069")
dev.copy(jpeg,
         here('results/plaice_forecast/akpegg_avgs',
              'plaice_egg_avg2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

# 2070-2099
df_akpegg_merged3 <- list(df_akpegg3_cesm126, df_akpegg3_cesm585,
                         df_akpegg3_gfdl126, df_akpegg3_gfdl585,
                         df_akpegg3_miroc126, df_akpegg3_miroc585) 

avg_akpegg_merged3 <- predict_avgs(df_akpegg_merged3)

windows(width = 6, height = 6, family = "serif")
grid_predict(avg_akpegg_merged3, "Forecasted Distribution 2070 - 2099")
dev.copy(jpeg,
         here('results/plaice_forecast/akpegg_avgs',
              'plaice_egg_avg3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

##### GIFs -------------------------------------------------------------------------------------------------------------------------
akpegg_dir_out <- file.path(base_dir, 'results', 'plaice_forecast', 'akpegg_avgs')
akpegg_imgs <- list.files(akpegg_dir_out, full.names = T)
akpegg_img_list <- lapply(akpegg_imgs, image_read)
akpegg_img_joined <- image_join(akpegg_img_list)
akpegg_img_animated <- image_animate(akpegg_img_joined, fps = 1)
image_write(image = akpegg_img_animated,
            path = here('results', 'plaice_forecast', "akpegg_avgs.gif"))

##### Clear environment -------------------------------------------------------------------------------------------------------------------------
rm(df_akpegg1_cesm126, df_akpegg1_cesm585,
   df_akpegg1_gfdl126, df_akpegg1_gfdl585,
   df_akpegg1_miroc126, df_akpegg1_miroc585,
   df_akpegg2_cesm126, df_akpegg2_cesm585,
   df_akpegg2_gfdl126, df_akpegg2_gfdl585,
   df_akpegg2_miroc126, df_akpegg2_miroc585,
   df_akpegg3_cesm126, df_akpegg3_cesm585,
   df_akpegg3_gfdl126, df_akpegg3_gfdl585,
   df_akpegg3_miroc126, df_akpegg3_miroc585,
   avg_akpegg_merged1, avg_akpegg_merged2,
   avg_akpegg_merged3, df_akpegg_merged1,
   df_akpegg_merged2, df_akpegg_merged3,
   akp_egg, akpegg_formula, akpegg_img_animated,
   akpegg_dir_out, akpegg_img_joined, 
   akpegg_img_list, akpegg_imgs)


### Alaska Plaice Larvae --------------------------------------------------------------------------------------------------------------------------
akp_larvae <- load_data('akp_larvae.rds', akp_larvae, roms_temps)
akplarvae_formula <- formula_geog(akp_larvae)

#### Forecast and average into 3 time periods ---------------------------------------------------------------------------------------------
##### CESM 126 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
cesm_temps1 <- readRDS(here('data', 'cesm_forecast_temp1.rds'))
cesm_salts1 <- readRDS(here('data', 'cesm_forecast_salt1.rds'))

df_akplarvae1_cesm126 <- predict_cells(2015:2039, akp_larvae, 140,
                                      5, 'ssp126', cesm_temps1,
                                      cesm_salts1, akplarvae_formula)
saveRDS(df_akplarvae1_cesm126, file = here("data", "df_akplarvae1_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae1_cesm126, 
             "Forecasted Distribution 2015 - 2039 \n CESM SSP1-2.6")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_cesm_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
cesm_temps2 <- readRDS(here('data', 'cesm_forecast_temp2.rds'))
cesm_salts2 <- readRDS(here('data', 'cesm_forecast_salt2.rds'))

df_akplarvae2_cesm126 <- predict_cells(2040:2069, akp_larvae, 140,
                                      5, 'ssp126', cesm_temps2, 
                                      cesm_salts2, akplarvae_formula)
saveRDS(df_akplarvae2_cesm126, file = here("data", "df_akplarvae2_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae2_cesm126, 
             "Forecasted Distribution 2040 - 2069 \n CESM SSP1-2.6")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_cesm_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
cesm_temps3 <- readRDS(here('data', 'cesm_forecast_temp3.rds'))
cesm_salts3 <- readRDS(here('data', 'cesm_forecast_salt3.rds'))

df_akplarvae3_cesm126 <- predict_cells(2070:2099, akp_larvae, 140,
                                      5, 'ssp126', cesm_temps3, 
                                      cesm_salts3, akplarvae_formula)
saveRDS(df_akplarvae3_cesm126, file = here("data", "df_akplarvae3_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae3_cesm126, 
             "Forecasted Distribution 2070 - 2099 \n CESM SSP1-2.6")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_cesm_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### CESM 585 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
df_akplarvae1_cesm585 <- predict_cells(2015:2039, akp_larvae, 140,
                                      5, 'ssp585', cesm_temps1, 
                                      cesm_salts1, akplarvae_formula)
saveRDS(df_akplarvae1_cesm585, file = here("data", "df_akplarvae1_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae1_cesm585, 
             "Forecasted Distribution 2015 - 2039 \n CESM SSP5-8.5")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_cesm_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
df_akplarvae2_cesm585 <- predict_cells(2040:2069, akp_larvae, 140,
                                      5, 'ssp585', cesm_temps2, 
                                      cesm_salts2, akplarvae_formula)
saveRDS(df_akplarvae2_cesm585, file = here("data", "df_akplarvae2_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae2_cesm585, 
             "Forecasted Distribution 2040 - 2069 \n CESM SSP5-8.5")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_cesm_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
df_akplarvae3_cesm585 <- predict_cells(2070:2099, akp_larvae, 140,
                                      5, 'ssp585', cesm_temps3, 
                                      cesm_salts3, akplarvae_formula)
saveRDS(df_akplarvae3_cesm585, file = here("data", "df_akplarvae3_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae3_cesm585, 
             "Forecasted Distribution 2070 - 2099 \n CESM SSP5-8.5")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_cesm_ssp585_3.jpg'),
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

df_akplarvae1_gfdl126 <- predict_cells(2015:2039, akp_larvae, 140,
                                      5, 'ssp126', gfdl_temps1, 
                                      gfdl_salts1, akplarvae_formula)
saveRDS(df_akplarvae1_gfdl126, file = here("data", "df_akplarvae1_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae1_gfdl126, 
             "Forecasted Distribution 2015 - 2039 \n GFDL SSP1-2.6")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_gfdl_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
gfdl_temps2 <- readRDS(here('data', 'gfdl_forecast_temp2.rds'))
gfdl_salts2 <- readRDS(here('data', 'gfdl_forecast_salt2.rds'))

df_akplarvae2_gfdl126 <- predict_cells(2040:2069, akp_larvae, 140,
                                      5, 'ssp126', gfdl_temps2, 
                                      gfdl_salts2, akplarvae_formula)
saveRDS(df_akplarvae2_gfdl126, file = here("data", "df_akplarvae2_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae2_gfdl126, 
             "Forecasted Distribution 2040 - 2069 \n GFDL SSP1-2.6")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_gfdl_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
gfdl_temps3 <- readRDS(here('data', 'gfdl_forecast_temp3.rds'))
gfdl_salts3 <- readRDS(here('data', 'gfdl_forecast_salt3.rds'))

df_akplarvae3_gfdl126 <- predict_cells(2070:2099, akp_larvae, 140,
                                      5, 'ssp126', gfdl_temps3, 
                                      gfdl_salts3, akplarvae_formula)
saveRDS(df_akplarvae3_gfdl126, file = here("data", "df_akplarvae3_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae3_gfdl126, 
             "Forecasted Distribution 2070 - 2099 \n GFDL SSP1-2.6")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_gfdl_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### GFDL 585 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
df_akplarvae1_gfdl585 <- predict_cells(2015:2039, akp_larvae, 140,
                                      5, 'ssp585', gfdl_temps1, 
                                      gfdl_salts1, akplarvae_formula)
saveRDS(df_akplarvae1_gfdl585, file = here("data", "df_akplarvae1_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae1_gfdl585, 
             "Forecasted Distribution 2015 - 2039 \n GFDL SSP5-8.5")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_gfdl_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
df_akplarvae2_gfdl585 <- predict_cells(2040:2069, akp_larvae, 140,
                                      5, 'ssp585', gfdl_temps2, 
                                      gfdl_salts2, akplarvae_formula)
saveRDS(df_akplarvae2_gfdl585, file = here("data", "df_akplarvae2_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae2_gfdl585, 
             "Forecasted Distribution 2040 - 2069 \n GFDL SSP5-8.5")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_gfdl_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
df_akplarvae3_gfdl585 <- predict_cells(2070:2099, akp_larvae, 140,
                                      5, 'ssp585', gfdl_temps3, 
                                      gfdl_salts3, akplarvae_formula)
saveRDS(df_akplarvae3_gfdl585, file = here("data", "df_akplarvae3_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae3_gfdl585, 
             "Forecasted Distribution 2070 - 2099 \n GFDL SSP5-8.5")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_gfdl_ssp585_3.jpg'),
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

df_akplarvae1_miroc126 <- predict_cells(2015:2039, akp_larvae, 140,
                                       5, 'ssp126', miroc_temps1, 
                                       miroc_salts1, akplarvae_formula)
saveRDS(df_akplarvae1_miroc126, file = here("data", "df_akplarvae1_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae1_miroc126, 
             "Forecasted Distribution 2015 - 2039 \n MIROC SSP1-2.6")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_miroc_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
miroc_temps2 <- readRDS(here('data', 'miroc_forecast_temp2.rds'))
miroc_salts2 <- readRDS(here('data', 'miroc_forecast_salt2.rds'))

df_akplarvae2_miroc126 <- predict_cells(2040:2069, akp_larvae, 140,
                                       5, 'ssp126', miroc_temps2, 
                                       miroc_salts2, akplarvae_formula)
saveRDS(df_akplarvae2_miroc126, file = here("data", "df_akplarvae2_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae2_miroc126, 
             "Forecasted Distribution 2040 - 2069 \n MIROC SSP1-2.6")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_miroc_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
miroc_temps3 <- readRDS(here('data', 'miroc_forecast_temp3.rds'))
miroc_salts3 <- readRDS(here('data', 'miroc_forecast_salt3.rds'))

df_akplarvae3_miroc126 <- predict_cells(2070:2099, akp_larvae, 140,
                                       5, 'ssp126', miroc_temps3, 
                                       miroc_salts3, akplarvae_formula)
saveRDS(df_akplarvae3_miroc126, file = here("data", "df_akplarvae3_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae3_miroc126, 
             "Forecasted Distribution 2070 - 2099 \n MIROC SSP1-2.6")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_miroc_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### MIROC 585 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
df_akplarvae1_miroc585 <- predict_cells(2015:2039, akp_larvae, 140,
                                       5, 'ssp585', miroc_temps1, 
                                       miroc_salts1, akplarvae_formula)
saveRDS(df_akplarvae1_miroc585, file = here("data", "df_akplarvae1_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae1_miroc585, 
             "Forecasted Distribution 2015 - 2039 \n MIROC SSP5-8.5")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_miroc_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
df_akplarvae2_miroc585 <- predict_cells(2040:2069, akp_larvae, 140,
                                       5, 'ssp585', miroc_temps2, 
                                       miroc_salts2, akplarvae_formula)
saveRDS(df_akplarvae2_miroc585, file = here("data", "df_akplarvae2_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae2_miroc585, 
             "Forecasted Distribution 2040 - 2069 \n MIROC SSP5-8.5")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_miroc_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
df_akplarvae3_miroc585 <- predict_cells(2070:2099, akp_larvae, 140,
                                       5, 'ssp585', miroc_temps3,
                                       miroc_salts3, akplarvae_formula)
saveRDS(df_akplarvae3_miroc585, file = here("data", "df_akplarvae3_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_akplarvae3_miroc585, 
             "Forecasted Distribution 2070 - 2099 \n MIROC SSP5-8.5")
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_miroc_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

rm(miroc_temps1, miroc_temps2, miroc_temps3,
   miroc_salts1, miroc_salts2, miroc_salts3)

graphics.off()

##### Multi-panel figure -------------------------------------------------------------------------------------------------------------------------
tiff(here('results/plaice_forecast',
          'plaicelarvae_multipanel.tiff'),
     units = "in",
     width = 45,
     height = 90,
     res = 300)
par(mfrow = c(6, 3),
    mar = c(11, 12, 5, 1.5) + 0.1,
    oma = c(3, 25, 15, 1),
    mgp = c(10, 4, 0),
    family = "serif")
grid_multipanel(df_akplarvae1_cesm126)
grid_multipanel(df_akplarvae2_cesm126)
grid_multipanel(df_akplarvae3_cesm126)
grid_multipanel(df_akplarvae1_cesm585)
grid_multipanel(df_akplarvae2_cesm585)
grid_multipanel(df_akplarvae3_cesm585)
grid_multipanel(df_akplarvae1_gfdl126)
grid_multipanel(df_akplarvae2_gfdl126)
grid_multipanel(df_akplarvae3_gfdl126)
grid_multipanel(df_akplarvae1_gfdl585)
grid_multipanel(df_akplarvae2_gfdl585)
grid_multipanel(df_akplarvae3_gfdl585)
grid_multipanel(df_akplarvae1_miroc126)
grid_multipanel(df_akplarvae2_miroc126)
grid_multipanel(df_akplarvae3_miroc126)
grid_multipanel(df_akplarvae1_miroc585)
grid_multipanel(df_akplarvae2_miroc585)
grid_multipanel(df_akplarvae3_miroc585)
mtext("CESM SSP1-2.6", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.92)
mtext("CESM SSP5-8.5", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.76)
mtext("GFDL1-2.6", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.585)
mtext("GFDL5-8.5", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.42)
mtext("MIROC1-2.6", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.26)
mtext("MIROC5-8.5", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.095)
mtext("2015-2039", 
      side = 3, 
      line = 2, 
      outer = TRUE, 
      cex = 7,
      at = 0.17)
mtext("2040-2069", 
      side = 3, 
      line = 2, 
      outer = TRUE, 
      cex = 7,
      at = 0.51)
mtext("2070-2099", 
      side = 3, 
      line = 2, 
      outer = TRUE, 
      cex = 7,
      at = 0.84)
dev.off()

##### Averages ---------------------------------------------------------------------------------------------------------------------------
# 2015 - 2039
df_akplarvae_merged1 <- list(df_akplarvae1_cesm126, df_akplarvae1_cesm585,
                            df_akplarvae1_gfdl126, df_akplarvae1_gfdl585,
                            df_akplarvae1_miroc126, df_akplarvae1_miroc585) 

avg_akplarvae_merged1 <- predict_avgs(df_akplarvae_merged1)

windows(width = 6, height = 6, family = "serif")
grid_predict(avg_akplarvae_merged1, "Forecasted Distribution 2015 - 2039")
dev.copy(jpeg,
         here('results/plaice_forecast/akplarvae_avgs',
              'plaice_larvae_avg1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

# 2040-2069
df_akplarvae_merged2 <- list(df_akplarvae2_cesm126, df_akplarvae2_cesm585,
                            df_akplarvae2_gfdl126, df_akplarvae2_gfdl585,
                            df_akplarvae2_miroc126, df_akplarvae2_miroc585) 

avg_akplarvae_merged2 <- predict_avgs(df_akplarvae_merged2)

windows(width = 6, height = 6, family = "serif")
grid_predict(avg_akplarvae_merged2, "Forecasted Distribution 2040 - 2069")
dev.copy(jpeg,
         here('results/plaice_forecast/akplarvae_avgs',
              'plaice_larvae_avg2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

# 2070-2099
df_akplarvae_merged3 <- list(df_akplarvae3_cesm126, df_akplarvae3_cesm585,
                            df_akplarvae3_gfdl126, df_akplarvae3_gfdl585,
                            df_akplarvae3_miroc126, df_akplarvae3_miroc585) 

avg_akplarvae_merged3 <- predict_avgs(df_akplarvae_merged3)

windows(width = 6, height = 6, family = "serif")
grid_predict(avg_akplarvae_merged3, "Forecasted Distribution 2070 - 2099")
dev.copy(jpeg,
         here('results/plaice_forecast/akplarvae_avgs',
              'plaice_larvae_avg3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

##### GIFs -------------------------------------------------------------------------------------------------------------------------
akplarvae_dir_out <- file.path(base_dir, 'results', 'plaice_forecast', 'akplarvae_avgs')
akplarvae_imgs <- list.files(akplarvae_dir_out, full.names = T)
akplarvae_img_list <- lapply(akplarvae_imgs, image_read)
akplarvae_img_joined <- image_join(akplarvae_img_list)
akplarvae_img_animated <- image_animate(akplarvae_img_joined, fps = 1)
image_write(image = akplarvae_img_animated,
            path = here('results', 'plaice_forecast', "akplarvae_avgs.gif"))

##### Clear environment -------------------------------------------------------------------------------------------------------------------------
rm(df_akplarvae1_cesm126, df_akplarvae1_cesm585,
   df_akplarvae1_gfdl126, df_akplarvae1_gfdl585,
   df_akplarvae1_miroc126, df_akplarvae1_miroc585,
   df_akplarvae2_cesm126, df_akplarvae2_cesm585,
   df_akplarvae2_gfdl126, df_akplarvae2_gfdl585,
   df_akplarvae2_miroc126, df_akplarvae2_miroc585,
   df_akplarvae3_cesm126, df_akplarvae3_cesm585,
   df_akplarvae3_gfdl126, df_akplarvae3_gfdl585,
   df_akplarvae3_miroc126, df_akplarvae3_miroc585,
   avg_akplarvae_merged1, avg_akplarvae_merged2,
   avg_akplarvae_merged3, df_akplarvae_merged1,
   df_akplarvae_merged2, df_akplarvae_merged3,
   akp_larvae, akplarvae_formula, akplarvae_img_animated,
   akplarvae_dir_out, akplarvae_img_joined, 
   akplarvae_img_list, akplarvae_imgs)


### Yellowfin Sole Larvae --------------------------------------------------------------------------------------------------------------------------
yfs_larvae <- load_data('yfs_larvae.rds', yfs_larvae, roms_temps)
yfslarvae_formula <- formula_geog(yfs_larvae)

#### Forecast and average into 3 time periods ---------------------------------------------------------------------------------------------
##### CESM 126 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
cesm_temps1 <- readRDS(here('data', 'cesm_forecast_temp1.rds'))
cesm_salts1 <- readRDS(here('data', 'cesm_forecast_salt1.rds'))

df_yfslarvae1_cesm126 <- predict_cells(2015:2039, yfs_larvae, 254,
                                       9, 'ssp126', cesm_temps1,
                                       cesm_salts1, yfslarvae_formula)
saveRDS(df_yfslarvae1_cesm126, file = here("data", "df_yfslarvae1_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae1_cesm126, 
             "Forecasted Distribution 2015 - 2039 \n CESM SSP1-2.6")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_cesm_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
cesm_temps2 <- readRDS(here('data', 'cesm_forecast_temp2.rds'))
cesm_salts2 <- readRDS(here('data', 'cesm_forecast_salt2.rds'))

df_yfslarvae2_cesm126 <- predict_cells(2040:2069, yfs_larvae, 254,
                                       9, 'ssp126', cesm_temps2, 
                                       cesm_salts2, yfslarvae_formula)
saveRDS(df_yfslarvae2_cesm126, file = here("data", "df_yfslarvae2_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae2_cesm126, 
             "Forecasted Distribution 2040 - 2069 \n CESM SSP1-2.6")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_cesm_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
cesm_temps3 <- readRDS(here('data', 'cesm_forecast_temp3.rds'))
cesm_salts3 <- readRDS(here('data', 'cesm_forecast_salt3.rds'))

df_yfslarvae3_cesm126 <- predict_cells(2070:2099, yfs_larvae, 254,
                                       9, 'ssp126', cesm_temps3, 
                                       cesm_salts3, yfslarvae_formula)
saveRDS(df_yfslarvae3_cesm126, file = here("data", "df_yfslarvae3_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae3_cesm126, 
             "Forecasted Distribution 2070 - 2099 \n CESM SSP1-2.6")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_cesm_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### CESM 585 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
df_yfslarvae1_cesm585 <- predict_cells(2015:2039, yfs_larvae, 254,
                                       9, 'ssp585', cesm_temps1, 
                                       cesm_salts1, yfslarvae_formula)
saveRDS(df_yfslarvae1_cesm585, file = here("data", "df_yfslarvae1_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae1_cesm585, 
             "Forecasted Distribution 2015 - 2039 \n CESM SSP5-8.5")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_cesm_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
df_yfslarvae2_cesm585 <- predict_cells(2040:2069, yfs_larvae, 254,
                                       9, 'ssp585', cesm_temps2, 
                                       cesm_salts2, yfslarvae_formula)
saveRDS(df_yfslarvae2_cesm585, file = here("data", "df_yfslarvae2_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae2_cesm585, 
             "Forecasted Distribution 2040 - 2069 \n CESM SSP5-8.5")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_cesm_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
df_yfslarvae3_cesm585 <- predict_cells(2070:2099, yfs_larvae, 254,
                                       9, 'ssp585', cesm_temps3, 
                                       cesm_salts3, yfslarvae_formula)
saveRDS(df_yfslarvae3_cesm585, file = here("data", "df_yfslarvae3_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae3_cesm585, 
             "Forecasted Distribution 2070 - 2099 \n CESM SSP5-8.5")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_cesm_ssp585_3.jpg'),
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

df_yfslarvae1_gfdl126 <- predict_cells(2015:2039, yfs_larvae, 254,
                                       9, 'ssp126', gfdl_temps1, 
                                       gfdl_salts1, yfslarvae_formula)
saveRDS(df_yfslarvae1_gfdl126, file = here("data", "df_yfslarvae1_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae1_gfdl126, 
             "Forecasted Distribution 2015 - 2039 \n GFDL SSP1-2.6")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_gfdl_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
gfdl_temps2 <- readRDS(here('data', 'gfdl_forecast_temp2.rds'))
gfdl_salts2 <- readRDS(here('data', 'gfdl_forecast_salt2.rds'))

df_yfslarvae2_gfdl126 <- predict_cells(2040:2069, yfs_larvae, 254,
                                       9, 'ssp126', gfdl_temps2, 
                                       gfdl_salts2, yfslarvae_formula)
saveRDS(df_yfslarvae2_gfdl126, file = here("data", "df_yfslarvae2_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae2_gfdl126, 
             "Forecasted Distribution 2040 - 2069 \n GFDL SSP1-2.6")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_gfdl_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
gfdl_temps3 <- readRDS(here('data', 'gfdl_forecast_temp3.rds'))
gfdl_salts3 <- readRDS(here('data', 'gfdl_forecast_salt3.rds'))

df_yfslarvae3_gfdl126 <- predict_cells(2070:2099, yfs_larvae, 254,
                                       9, 'ssp126', gfdl_temps3, 
                                       gfdl_salts3, yfslarvae_formula)
saveRDS(df_yfslarvae3_gfdl126, file = here("data", "df_yfslarvae3_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae3_gfdl126, 
             "Forecasted Distribution 2070 - 2099 \n GFDL SSP1-2.6")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_gfdl_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### GFDL 585 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
df_yfslarvae1_gfdl585 <- predict_cells(2015:2039, yfs_larvae, 254,
                                       9, 'ssp585', gfdl_temps1, 
                                       gfdl_salts1, yfslarvae_formula)
saveRDS(df_yfslarvae1_gfdl585, file = here("data", "df_yfslarvae1_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae1_gfdl585, 
             "Forecasted Distribution 2015 - 2039 \n GFDL SSP5-8.5")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_gfdl_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
df_yfslarvae2_gfdl585 <- predict_cells(2040:2069, yfs_larvae, 254,
                                       9, 'ssp585', gfdl_temps2, 
                                       gfdl_salts2, yfslarvae_formula)
saveRDS(df_yfslarvae2_gfdl585, file = here("data", "df_yfslarvae2_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae2_gfdl585, 
             "Forecasted Distribution 2040 - 2069 \n GFDL SSP5-8.5")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_gfdl_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
df_yfslarvae3_gfdl585 <- predict_cells(2070:2099, yfs_larvae, 254,
                                       9, 'ssp585', gfdl_temps3, 
                                       gfdl_salts3, yfslarvae_formula)
saveRDS(df_yfslarvae3_gfdl585, file = here("data", "df_yfslarvae3_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae3_gfdl585, 
             "Forecasted Distribution 2070 - 2099 \n GFDL SSP5-8.5")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_gfdl_ssp585_3.jpg'),
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

df_yfslarvae1_miroc126 <- predict_cells(2015:2039, yfs_larvae, 254,
                                        9, 'ssp126', miroc_temps1, 
                                        miroc_salts1, yfslarvae_formula)
saveRDS(df_yfslarvae1_miroc126, file = here("data", "df_yfslarvae1_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae1_miroc126, 
             "Forecasted Distribution 2015 - 2039 \n MIROC SSP1-2.6")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_miroc_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
miroc_temps2 <- readRDS(here('data', 'miroc_forecast_temp2.rds'))
miroc_salts2 <- readRDS(here('data', 'miroc_forecast_salt2.rds'))

df_yfslarvae2_miroc126 <- predict_cells(2040:2069, yfs_larvae, 254,
                                        9, 'ssp126', miroc_temps2, 
                                        miroc_salts2, yfslarvae_formula)
saveRDS(df_yfslarvae2_miroc126, file = here("data", "df_yfslarvae2_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae2_miroc126, 
             "Forecasted Distribution 2040 - 2069 \n MIROC SSP1-2.6")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_miroc_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
miroc_temps3 <- readRDS(here('data', 'miroc_forecast_temp3.rds'))
miroc_salts3 <- readRDS(here('data', 'miroc_forecast_salt3.rds'))

df_yfslarvae3_miroc126 <- predict_cells(2070:2099, yfs_larvae, 254,
                                        9, 'ssp126', miroc_temps3, 
                                        miroc_salts3, yfslarvae_formula)
saveRDS(df_yfslarvae3_miroc126, file = here("data", "df_yfslarvae3_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae3_miroc126, 
             "Forecasted Distribution 2070 - 2099 \n MIROC SSP1-2.6")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_miroc_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### MIROC 585 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
df_yfslarvae1_miroc585 <- predict_cells(2015:2039, yfs_larvae, 254,
                                        9, 'ssp585', miroc_temps1, 
                                        miroc_salts1, yfslarvae_formula)
saveRDS(df_yfslarvae1_miroc585, file = here("data", "df_yfslarvae1_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae1_miroc585, 
             "Forecasted Distribution 2015 - 2039 \n MIROC SSP5-8.5")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_miroc_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
df_yfslarvae2_miroc585 <- predict_cells(2040:2069, yfs_larvae, 254,
                                        9, 'ssp585', miroc_temps2, 
                                        miroc_salts2, yfslarvae_formula)
saveRDS(df_yfslarvae2_miroc585, file = here("data", "df_yfslarvae2_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae2_miroc585, 
             "Forecasted Distribution 2040 - 2069 \n MIROC SSP5-8.5")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_miroc_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
df_yfslarvae3_miroc585 <- predict_cells(2070:2099, yfs_larvae, 254,
                                        9, 'ssp585', miroc_temps3,
                                        miroc_salts3, yfslarvae_formula)
saveRDS(df_yfslarvae3_miroc585, file = here("data", "df_yfslarvae3_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_yfslarvae3_miroc585, 
             "Forecasted Distribution 2070 - 2099 \n MIROC SSP5-8.5")
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_miroc_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

rm(miroc_temps1, miroc_temps2, miroc_temps3,
   miroc_salts1, miroc_salts2, miroc_salts3)

graphics.off()

##### Multi-panel figure -------------------------------------------------------------------------------------------------------------------------
tiff(here('results/yellowfin_forecast',
          'yellowfinlarvae_multipanel.tiff'),
     units = "in",
     width = 45,
     height = 90,
     res = 300)
par(mfrow = c(6, 3),
    mar = c(11, 12, 5, 1.5) + 0.1,
    oma = c(3, 25, 15, 1),
    mgp = c(10, 4, 0),
    family = "serif")
grid_multipanel(df_yfslarvae1_cesm126)
grid_multipanel(df_yfslarvae2_cesm126)
grid_multipanel(df_yfslarvae3_cesm126)
grid_multipanel(df_yfslarvae1_cesm585)
grid_multipanel(df_yfslarvae2_cesm585)
grid_multipanel(df_yfslarvae3_cesm585)
grid_multipanel(df_yfslarvae1_gfdl126)
grid_multipanel(df_yfslarvae2_gfdl126)
grid_multipanel(df_yfslarvae3_gfdl126)
grid_multipanel(df_yfslarvae1_gfdl585)
grid_multipanel(df_yfslarvae2_gfdl585)
grid_multipanel(df_yfslarvae3_gfdl585)
grid_multipanel(df_yfslarvae1_miroc126)
grid_multipanel(df_yfslarvae2_miroc126)
grid_multipanel(df_yfslarvae3_miroc126)
grid_multipanel(df_yfslarvae1_miroc585)
grid_multipanel(df_yfslarvae2_miroc585)
grid_multipanel(df_yfslarvae3_miroc585)
mtext("CESM SSP1-2.6", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.92)
mtext("CESM SSP5-8.5", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.76)
mtext("GFDL1-2.6", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.585)
mtext("GFDL5-8.5", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.42)
mtext("MIROC1-2.6", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.26)
mtext("MIROC5-8.5", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.095)
mtext("2015-2039", 
      side = 3, 
      line = 2, 
      outer = TRUE, 
      cex = 7,
      at = 0.17)
mtext("2040-2069", 
      side = 3, 
      line = 2, 
      outer = TRUE, 
      cex = 7,
      at = 0.51)
mtext("2070-2099", 
      side = 3, 
      line = 2, 
      outer = TRUE, 
      cex = 7,
      at = 0.84)
dev.off()

##### Averages ---------------------------------------------------------------------------------------------------------------------------
# 2015 - 2039
df_yfslarvae_merged1 <- list(df_yfslarvae1_cesm126, df_yfslarvae1_cesm585,
                             df_yfslarvae1_gfdl126, df_yfslarvae1_gfdl585,
                             df_yfslarvae1_miroc126, df_yfslarvae1_miroc585) 

avg_yfslarvae_merged1 <- predict_avgs(df_yfslarvae_merged1)

windows(width = 6, height = 6, family = "serif")
grid_predict(avg_yfslarvae_merged1, "Forecasted Distribution 2015 - 2039")
dev.copy(jpeg,
         here('results/yellowfin_forecast/yfslarvae_avgs',
              'yellowfin_larvae_avg1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

# 2040-2069
df_yfslarvae_merged2 <- list(df_yfslarvae2_cesm126, df_yfslarvae2_cesm585,
                             df_yfslarvae2_gfdl126, df_yfslarvae2_gfdl585,
                             df_yfslarvae2_miroc126, df_yfslarvae2_miroc585) 

avg_yfslarvae_merged2 <- predict_avgs(df_yfslarvae_merged2)

windows(width = 6, height = 6, family = "serif")
grid_predict(avg_yfslarvae_merged2, "Forecasted Distribution 2040 - 2069")
dev.copy(jpeg,
         here('results/yellowfin_forecast/yfslarvae_avgs',
              'yellowfin_larvae_avg2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

# 2070-2099
df_yfslarvae_merged3 <- list(df_yfslarvae3_cesm126, df_yfslarvae3_cesm585,
                             df_yfslarvae3_gfdl126, df_yfslarvae3_gfdl585,
                             df_yfslarvae3_miroc126, df_yfslarvae3_miroc585) 

avg_yfslarvae_merged3 <- predict_avgs(df_yfslarvae_merged3)

windows(width = 6, height = 6, family = "serif")
grid_predict(avg_yfslarvae_merged3, "Forecasted Distribution 2070 - 2099")
dev.copy(jpeg,
         here('results/yellowfin_forecast/yfslarvae_avgs',
              'yellowfin_larvae_avg3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

##### GIFs -------------------------------------------------------------------------------------------------------------------------
yfslarvae_dir_out <- file.path(base_dir, 'results', 'yellowfin_forecast', 'yfslarvae_avgs')
yfslarvae_imgs <- list.files(yfslarvae_dir_out, full.names = T)
yfslarvae_img_list <- lapply(yfslarvae_imgs, image_read)
yfslarvae_img_joined <- image_join(yfslarvae_img_list)
yfslarvae_img_animated <- image_animate(yfslarvae_img_joined, fps = 1)
image_write(image = yfslarvae_img_animated,
            path = here('results', 'yellowfin_forecast', "yfslarvae_avgs.gif"))

##### Clear environment -------------------------------------------------------------------------------------------------------------------------
rm(df_yfslarvae1_cesm126, df_yfslarvae1_cesm585,
   df_yfslarvae1_gfdl126, df_yfslarvae1_gfdl585,
   df_yfslarvae1_miroc126, df_yfslarvae1_miroc585,
   df_yfslarvae2_cesm126, df_yfslarvae2_cesm585,
   df_yfslarvae2_gfdl126, df_yfslarvae2_gfdl585,
   df_yfslarvae2_miroc126, df_yfslarvae2_miroc585,
   df_yfslarvae3_cesm126, df_yfslarvae3_cesm585,
   df_yfslarvae3_gfdl126, df_yfslarvae3_gfdl585,
   df_yfslarvae3_miroc126, df_yfslarvae3_miroc585,
   avg_yfslarvae_merged1, avg_yfslarvae_merged2,
   avg_yfslarvae_merged3, df_yfslarvae_merged1,
   df_yfslarvae_merged2, df_yfslarvae_merged3,
   yfs_larvae, yfslarvae_formula, yfslarvae_img_animated,
   yfslarvae_dir_out, yfslarvae_img_joined, 
   yfslarvae_img_list, yfslarvae_imgs)


### Northern Rock Sole Larvae --------------------------------------------------------------------------------------------------------------------------
nrs_larvae <- load_data('nrs_larvae.rds', nrs_larvae, roms_temps)
nrslarvae_formula <- formula_geog(nrs_larvae)

#### Forecast and average into 3 time periods ---------------------------------------------------------------------------------------------
##### CESM 126 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
cesm_temps1 <- readRDS(here('data', 'cesm_forecast_temp1.rds'))
cesm_salts1 <- readRDS(here('data', 'cesm_forecast_salt1.rds'))

df_nrslarvae1_cesm126 <- predict_cells(2015:2039, nrs_larvae, 134,
                                       5, 'ssp126', cesm_temps1,
                                       cesm_salts1, nrslarvae_formula)
saveRDS(df_nrslarvae1_cesm126, file = here("data", "df_nrslarvae1_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae1_cesm126, 
             "Forecasted Distribution 2015 - 2039 \n CESM SSP1-2.6")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_cesm_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
cesm_temps2 <- readRDS(here('data', 'cesm_forecast_temp2.rds'))
cesm_salts2 <- readRDS(here('data', 'cesm_forecast_salt2.rds'))

df_nrslarvae2_cesm126 <- predict_cells(2040:2069, nrs_larvae, 134,
                                       5, 'ssp126', cesm_temps2, 
                                       cesm_salts2, nrslarvae_formula)
saveRDS(df_nrslarvae2_cesm126, file = here("data", "df_nrslarvae2_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae2_cesm126, 
             "Forecasted Distribution 2040 - 2069 \n CESM SSP1-2.6")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_cesm_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
cesm_temps3 <- readRDS(here('data', 'cesm_forecast_temp3.rds'))
cesm_salts3 <- readRDS(here('data', 'cesm_forecast_salt3.rds'))

df_nrslarvae3_cesm126 <- predict_cells(2070:2099, nrs_larvae, 134,
                                       5, 'ssp126', cesm_temps3, 
                                       cesm_salts3, nrslarvae_formula)
saveRDS(df_nrslarvae3_cesm126, file = here("data", "df_nrslarvae3_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae3_cesm126, 
             "Forecasted Distribution 2070 - 2099 \n CESM SSP1-2.6")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_cesm_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### CESM 585 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
df_nrslarvae1_cesm585 <- predict_cells(2015:2039, nrs_larvae, 134,
                                       5, 'ssp585', cesm_temps1, 
                                       cesm_salts1, nrslarvae_formula)
saveRDS(df_nrslarvae1_cesm585, file = here("data", "df_nrslarvae1_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae1_cesm585, 
             "Forecasted Distribution 2015 - 2039 \n CESM SSP5-8.5")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_cesm_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
df_nrslarvae2_cesm585 <- predict_cells(2040:2069, nrs_larvae, 134,
                                       5, 'ssp585', cesm_temps2, 
                                       cesm_salts2, nrslarvae_formula)
saveRDS(df_nrslarvae2_cesm585, file = here("data", "df_nrslarvae2_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae2_cesm585, 
             "Forecasted Distribution 2040 - 2069 \n CESM SSP5-8.5")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_cesm_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
df_nrslarvae3_cesm585 <- predict_cells(2070:2099, nrs_larvae, 134,
                                       5, 'ssp585', cesm_temps3, 
                                       cesm_salts3, nrslarvae_formula)
saveRDS(df_nrslarvae3_cesm585, file = here("data", "df_nrslarvae3_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae3_cesm585, 
             "Forecasted Distribution 2070 - 2099 \n CESM SSP5-8.5")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_cesm_ssp585_3.jpg'),
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

df_nrslarvae1_gfdl126 <- predict_cells(2015:2039, nrs_larvae, 134,
                                       5, 'ssp126', gfdl_temps1, 
                                       gfdl_salts1, nrslarvae_formula)
saveRDS(df_nrslarvae1_gfdl126, file = here("data", "df_nrslarvae1_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae1_gfdl126, 
             "Forecasted Distribution 2015 - 2039 \n GFDL SSP1-2.6")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_gfdl_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
gfdl_temps2 <- readRDS(here('data', 'gfdl_forecast_temp2.rds'))
gfdl_salts2 <- readRDS(here('data', 'gfdl_forecast_salt2.rds'))

df_nrslarvae2_gfdl126 <- predict_cells(2040:2069, nrs_larvae, 134,
                                       5, 'ssp126', gfdl_temps2, 
                                       gfdl_salts2, nrslarvae_formula)
saveRDS(df_nrslarvae2_gfdl126, file = here("data", "df_nrslarvae2_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae2_gfdl126, 
             "Forecasted Distribution 2040 - 2069 \n GFDL SSP1-2.6")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_gfdl_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
gfdl_temps3 <- readRDS(here('data', 'gfdl_forecast_temp3.rds'))
gfdl_salts3 <- readRDS(here('data', 'gfdl_forecast_salt3.rds'))

df_nrslarvae3_gfdl126 <- predict_cells(2070:2099, nrs_larvae, 134,
                                       5, 'ssp126', gfdl_temps3, 
                                       gfdl_salts3, nrslarvae_formula)
saveRDS(df_nrslarvae3_gfdl126, file = here("data", "df_nrslarvae3_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae3_gfdl126, 
             "Forecasted Distribution 2070 - 2099 \n GFDL SSP1-2.6")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_gfdl_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### GFDL 585 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
df_nrslarvae1_gfdl585 <- predict_cells(2015:2039, nrs_larvae, 134,
                                       5, 'ssp585', gfdl_temps1, 
                                       gfdl_salts1, nrslarvae_formula)
saveRDS(df_nrslarvae1_gfdl585, file = here("data", "df_nrslarvae1_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae1_gfdl585, 
             "Forecasted Distribution 2015 - 2039 \n GFDL SSP5-8.5")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_gfdl_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
df_nrslarvae2_gfdl585 <- predict_cells(2040:2069, nrs_larvae, 134,
                                       5, 'ssp585', gfdl_temps2, 
                                       gfdl_salts2, nrslarvae_formula)
saveRDS(df_nrslarvae2_gfdl585, file = here("data", "df_nrslarvae2_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae2_gfdl585, 
             "Forecasted Distribution 2040 - 2069 \n GFDL SSP5-8.5")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_gfdl_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
df_nrslarvae3_gfdl585 <- predict_cells(2070:2099, nrs_larvae, 134,
                                       5, 'ssp585', gfdl_temps3, 
                                       gfdl_salts3, nrslarvae_formula)
saveRDS(df_nrslarvae3_gfdl585, file = here("data", "df_nrslarvae3_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae3_gfdl585, 
             "Forecasted Distribution 2070 - 2099 \n GFDL SSP5-8.5")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_gfdl_ssp585_3.jpg'),
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

df_nrslarvae1_miroc126 <- predict_cells(2015:2039, nrs_larvae, 134,
                                        5, 'ssp126', miroc_temps1, 
                                        miroc_salts1, nrslarvae_formula)
saveRDS(df_nrslarvae1_miroc126, file = here("data", "df_nrslarvae1_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae1_miroc126, 
             "Forecasted Distribution 2015 - 2039 \n MIROC SSP1-2.6")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_miroc_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
miroc_temps2 <- readRDS(here('data', 'miroc_forecast_temp2.rds'))
miroc_salts2 <- readRDS(here('data', 'miroc_forecast_salt2.rds'))

df_nrslarvae2_miroc126 <- predict_cells(2040:2069, nrs_larvae, 134,
                                        5, 'ssp126', miroc_temps2, 
                                        miroc_salts2, nrslarvae_formula)
saveRDS(df_nrslarvae2_miroc126, file = here("data", "df_nrslarvae2_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae2_miroc126, 
             "Forecasted Distribution 2040 - 2069 \n MIROC SSP1-2.6")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_miroc_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
miroc_temps3 <- readRDS(here('data', 'miroc_forecast_temp3.rds'))
miroc_salts3 <- readRDS(here('data', 'miroc_forecast_salt3.rds'))

df_nrslarvae3_miroc126 <- predict_cells(2070:2099, nrs_larvae, 134,
                                        5, 'ssp126', miroc_temps3, 
                                        miroc_salts3, nrslarvae_formula)
saveRDS(df_nrslarvae3_miroc126, file = here("data", "df_nrslarvae3_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae3_miroc126, 
             "Forecasted Distribution 2070 - 2099 \n MIROC SSP1-2.6")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_miroc_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### MIROC 585 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
df_nrslarvae1_miroc585 <- predict_cells(2015:2039, nrs_larvae, 134,
                                        5, 'ssp585', miroc_temps1, 
                                        miroc_salts1, nrslarvae_formula)
saveRDS(df_nrslarvae1_miroc585, file = here("data", "df_nrslarvae1_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae1_miroc585, 
             "Forecasted Distribution 2015 - 2039 \n MIROC SSP5-8.5")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_miroc_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
df_nrslarvae2_miroc585 <- predict_cells(2040:2069, nrs_larvae, 134,
                                        5, 'ssp585', miroc_temps2, 
                                        miroc_salts2, nrslarvae_formula)
saveRDS(df_nrslarvae2_miroc585, file = here("data", "df_nrslarvae2_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae2_miroc585, 
             "Forecasted Distribution 2040 - 2069 \n MIROC SSP5-8.5")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_miroc_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
df_nrslarvae3_miroc585 <- predict_cells(2070:2099, nrs_larvae, 134,
                                        5, 'ssp585', miroc_temps3,
                                        miroc_salts3, nrslarvae_formula)
saveRDS(df_nrslarvae3_miroc585, file = here("data", "df_nrslarvae3_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_nrslarvae3_miroc585, 
             "Forecasted Distribution 2070 - 2099 \n MIROC SSP5-8.5")
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_miroc_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

rm(miroc_temps1, miroc_temps2, miroc_temps3,
   miroc_salts1, miroc_salts2, miroc_salts3)

graphics.off()

##### Multi-panel figure -------------------------------------------------------------------------------------------------------------------------
tiff(here('results/rocksole_forecast',
          'rocksolelarvae_multipanel.tiff'),
     units = "in",
     width = 45,
     height = 90,
     res = 300)
par(mfrow = c(6, 3),
    mar = c(11, 12, 5, 1.5) + 0.1,
    oma = c(3, 25, 15, 1),
    mgp = c(10, 4, 0),
    family = "serif")
grid_multipanel(df_nrslarvae1_cesm126)
grid_multipanel(df_nrslarvae2_cesm126)
grid_multipanel(df_nrslarvae3_cesm126)
grid_multipanel(df_nrslarvae1_cesm585)
grid_multipanel(df_nrslarvae2_cesm585)
grid_multipanel(df_nrslarvae3_cesm585)
grid_multipanel(df_nrslarvae1_gfdl126)
grid_multipanel(df_nrslarvae2_gfdl126)
grid_multipanel(df_nrslarvae3_gfdl126)
grid_multipanel(df_nrslarvae1_gfdl585)
grid_multipanel(df_nrslarvae2_gfdl585)
grid_multipanel(df_nrslarvae3_gfdl585)
grid_multipanel(df_nrslarvae1_miroc126)
grid_multipanel(df_nrslarvae2_miroc126)
grid_multipanel(df_nrslarvae3_miroc126)
grid_multipanel(df_nrslarvae1_miroc585)
grid_multipanel(df_nrslarvae2_miroc585)
grid_multipanel(df_nrslarvae3_miroc585)
mtext("CESM SSP1-2.6", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.92)
mtext("CESM SSP5-8.5", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.76)
mtext("GFDL1-2.6", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.585)
mtext("GFDL5-8.5", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.42)
mtext("MIROC1-2.6", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.26)
mtext("MIROC5-8.5", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.095)
mtext("2015-2039", 
      side = 3, 
      line = 2, 
      outer = TRUE, 
      cex = 7,
      at = 0.17)
mtext("2040-2069", 
      side = 3, 
      line = 2, 
      outer = TRUE, 
      cex = 7,
      at = 0.51)
mtext("2070-2099", 
      side = 3, 
      line = 2, 
      outer = TRUE, 
      cex = 7,
      at = 0.84)
dev.off()

##### Averages ---------------------------------------------------------------------------------------------------------------------------
# 2015 - 2039
df_nrslarvae_merged1 <- list(df_nrslarvae1_cesm126, df_nrslarvae1_cesm585,
                             df_nrslarvae1_gfdl126, df_nrslarvae1_gfdl585,
                             df_nrslarvae1_miroc126, df_nrslarvae1_miroc585) 

avg_nrslarvae_merged1 <- predict_avgs(df_nrslarvae_merged1)

windows(width = 6, height = 6, family = "serif")
grid_predict(avg_nrslarvae_merged1, "Forecasted Distribution 2015 - 2039")
dev.copy(jpeg,
         here('results/rocksole_forecast/nrslarvae_avgs',
              'rocksole_larvae_avg1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

# 2040-2069
df_nrslarvae_merged2 <- list(df_nrslarvae2_cesm126, df_nrslarvae2_cesm585,
                             df_nrslarvae2_gfdl126, df_nrslarvae2_gfdl585,
                             df_nrslarvae2_miroc126, df_nrslarvae2_miroc585) 

avg_nrslarvae_merged2 <- predict_avgs(df_nrslarvae_merged2)

windows(width = 6, height = 6, family = "serif")
grid_predict(avg_nrslarvae_merged2, "Forecasted Distribution 2040 - 2069")
dev.copy(jpeg,
         here('results/rocksole_forecast/nrslarvae_avgs',
              'rocksole_larvae_avg2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

# 2070-2099
df_nrslarvae_merged3 <- list(df_nrslarvae3_cesm126, df_nrslarvae3_cesm585,
                             df_nrslarvae3_gfdl126, df_nrslarvae3_gfdl585,
                             df_nrslarvae3_miroc126, df_nrslarvae3_miroc585) 

avg_nrslarvae_merged3 <- predict_avgs(df_nrslarvae_merged3)

windows(width = 6, height = 6, family = "serif")
grid_predict(avg_nrslarvae_merged3, "Forecasted Distribution 2070 - 2099")
dev.copy(jpeg,
         here('results/rocksole_forecast/nrslarvae_avgs',
              'rocksole_larvae_avg3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

##### GIFs -------------------------------------------------------------------------------------------------------------------------
nrslarvae_dir_out <- file.path(base_dir, 'results', 'rocksole_forecast', 'nrslarvae_avgs')
nrslarvae_imgs <- list.files(nrslarvae_dir_out, full.names = T)
nrslarvae_img_list <- lapply(nrslarvae_imgs, image_read)
nrslarvae_img_joined <- image_join(nrslarvae_img_list)
nrslarvae_img_animated <- image_animate(nrslarvae_img_joined, fps = 1)
image_write(image = nrslarvae_img_animated,
            path = here('results', 'rocksole_forecast', "nrslarvae_avgs.gif"))

##### Clear environment -------------------------------------------------------------------------------------------------------------------------
rm(df_nrslarvae1_cesm126, df_nrslarvae1_cesm585,
   df_nrslarvae1_gfdl126, df_nrslarvae1_gfdl585,
   df_nrslarvae1_miroc126, df_nrslarvae1_miroc585,
   df_nrslarvae2_cesm126, df_nrslarvae2_cesm585,
   df_nrslarvae2_gfdl126, df_nrslarvae2_gfdl585,
   df_nrslarvae2_miroc126, df_nrslarvae2_miroc585,
   df_nrslarvae3_cesm126, df_nrslarvae3_cesm585,
   df_nrslarvae3_gfdl126, df_nrslarvae3_gfdl585,
   df_nrslarvae3_miroc126, df_nrslarvae3_miroc585,
   avg_nrslarvae_merged1, avg_nrslarvae_merged2,
   avg_nrslarvae_merged3, df_nrslarvae_merged1,
   df_nrslarvae_merged2, df_nrslarvae_merged3,
   nrs_larvae, nrslarvae_formula, nrslarvae_img_animated,
   nrslarvae_dir_out, nrslarvae_img_joined, 
   nrslarvae_img_list, nrslarvae_imgs)

### Pacific Cod Larvae --------------------------------------------------------------------------------------------------------------------------
pcod_larvae <- load_data('pcod_larvae.rds', pcod_larvae, roms_temps)
pcodlarvae_formula <- formula_pheno(pcod_larvae)

#### Forecast and average into 3 time periods ---------------------------------------------------------------------------------------------
##### CESM 126 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
cesm_temps1 <- readRDS(here('data', 'cesm_forecast_temp1.rds'))
cesm_salts1 <- readRDS(here('data', 'cesm_forecast_salt1.rds'))

df_pcodlarvae1_cesm126 <- predict_cells(2015:2039, pcod_larvae, 134,
                                       5, 'ssp126', cesm_temps1,
                                       cesm_salts1, pcodlarvae_formula)
saveRDS(df_pcodlarvae1_cesm126, file = here("data", "df_pcodlarvae1_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pcodlarvae1_cesm126, 
             "Forecasted Distribution 2015 - 2039 \n CESM SSP1-2.6")
dev.copy(jpeg,
         here('results/cod_forecast',
              'cod_larvae_cesm_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
cesm_temps2 <- readRDS(here('data', 'cesm_forecast_temp2.rds'))
cesm_salts2 <- readRDS(here('data', 'cesm_forecast_salt2.rds'))

df_pcodlarvae2_cesm126 <- predict_cells(2040:2069, pcod_larvae, 134,
                                       5, 'ssp126', cesm_temps2, 
                                       cesm_salts2, pcodlarvae_formula)
saveRDS(df_pcodlarvae2_cesm126, file = here("data", "df_pcodlarvae2_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pcodlarvae2_cesm126, 
             "Forecasted Distribution 2040 - 2069 \n CESM SSP1-2.6")
dev.copy(jpeg,
         here('results/cod_forecast',
              'cod_larvae_cesm_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
cesm_temps3 <- readRDS(here('data', 'cesm_forecast_temp3.rds'))
cesm_salts3 <- readRDS(here('data', 'cesm_forecast_salt3.rds'))

df_pcodlarvae3_cesm126 <- predict_cells(2070:2099, pcod_larvae, 134,
                                       5, 'ssp126', cesm_temps3, 
                                       cesm_salts3, pcodlarvae_formula)
saveRDS(df_pcodlarvae3_cesm126, file = here("data", "df_pcodlarvae3_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pcodlarvae3_cesm126, 
             "Forecasted Distribution 2070 - 2099 \n CESM SSP1-2.6")
dev.copy(jpeg,
         here('results/cod_forecast',
              'cod_larvae_cesm_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### CESM 585 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
df_pcodlarvae1_cesm585 <- predict_cells(2015:2039, pcod_larvae, 134,
                                       5, 'ssp585', cesm_temps1, 
                                       cesm_salts1, pcodlarvae_formula)
saveRDS(df_pcodlarvae1_cesm585, file = here("data", "df_pcodlarvae1_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pcodlarvae1_cesm585, 
             "Forecasted Distribution 2015 - 2039 \n CESM SSP5-8.5")
dev.copy(jpeg,
         here('results/cod_forecast',
              'cod_larvae_cesm_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
df_pcodlarvae2_cesm585 <- predict_cells(2040:2069, pcod_larvae, 134,
                                       5, 'ssp585', cesm_temps2, 
                                       cesm_salts2, pcodlarvae_formula)
saveRDS(df_pcodlarvae2_cesm585, file = here("data", "df_pcodlarvae2_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pcodlarvae2_cesm585, 
             "Forecasted Distribution 2040 - 2069 \n CESM SSP5-8.5")
dev.copy(jpeg,
         here('results/cod_forecast',
              'cod_larvae_cesm_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
df_pcodlarvae3_cesm585 <- predict_cells(2070:2099, pcod_larvae, 134,
                                       5, 'ssp585', cesm_temps3, 
                                       cesm_salts3, pcodlarvae_formula)
saveRDS(df_pcodlarvae3_cesm585, file = here("data", "df_pcodlarvae3_cesm585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pcodlarvae3_cesm585, 
             "Forecasted Distribution 2070 - 2099 \n CESM SSP5-8.5")
dev.copy(jpeg,
         here('results/cod_forecast',
              'cod_larvae_cesm_ssp585_3.jpg'),
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

df_pcodlarvae1_gfdl126 <- predict_cells(2015:2039, pcod_larvae, 134,
                                       5, 'ssp126', gfdl_temps1, 
                                       gfdl_salts1, pcodlarvae_formula)
saveRDS(df_pcodlarvae1_gfdl126, file = here("data", "df_pcodlarvae1_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pcodlarvae1_gfdl126, 
             "Forecasted Distribution 2015 - 2039 \n GFDL SSP1-2.6")
dev.copy(jpeg,
         here('results/cod_forecast',
              'cod_larvae_gfdl_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
gfdl_temps2 <- readRDS(here('data', 'gfdl_forecast_temp2.rds'))
gfdl_salts2 <- readRDS(here('data', 'gfdl_forecast_salt2.rds'))

df_pcodlarvae2_gfdl126 <- predict_cells(2040:2069, pcod_larvae, 134,
                                       5, 'ssp126', gfdl_temps2, 
                                       gfdl_salts2, pcodlarvae_formula)
saveRDS(df_pcodlarvae2_gfdl126, file = here("data", "df_pcodlarvae2_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pcodlarvae2_gfdl126, 
             "Forecasted Distribution 2040 - 2069 \n GFDL SSP1-2.6")
dev.copy(jpeg,
         here('results/cod_forecast',
              'cod_larvae_gfdl_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
gfdl_temps3 <- readRDS(here('data', 'gfdl_forecast_temp3.rds'))
gfdl_salts3 <- readRDS(here('data', 'gfdl_forecast_salt3.rds'))

df_pcodlarvae3_gfdl126 <- predict_cells(2070:2099, pcod_larvae, 134,
                                       5, 'ssp126', gfdl_temps3, 
                                       gfdl_salts3, pcodlarvae_formula)
saveRDS(df_pcodlarvae3_gfdl126, file = here("data", "df_pcodlarvae3_gfdl126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pcodlarvae3_gfdl126, 
             "Forecasted Distribution 2070 - 2099 \n GFDL SSP1-2.6")
dev.copy(jpeg,
         here('results/cod_forecast',
              'cod_larvae_gfdl_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### GFDL 585 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
df_pcodlarvae1_gfdl585 <- predict_cells(2015:2039, pcod_larvae, 134,
                                       5, 'ssp585', gfdl_temps1, 
                                       gfdl_salts1, pcodlarvae_formula)
saveRDS(df_pcodlarvae1_gfdl585, file = here("data", "df_pcodlarvae1_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pcodlarvae1_gfdl585, 
             "Forecasted Distribution 2015 - 2039 \n GFDL SSP5-8.5")
dev.copy(jpeg,
         here('results/cod_forecast',
              'cod_larvae_gfdl_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
df_pcodlarvae2_gfdl585 <- predict_cells(2040:2069, pcod_larvae, 134,
                                       5, 'ssp585', gfdl_temps2, 
                                       gfdl_salts2, pcodlarvae_formula)
saveRDS(df_pcodlarvae2_gfdl585, file = here("data", "df_pcodlarvae2_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pcodlarvae2_gfdl585, 
             "Forecasted Distribution 2040 - 2069 \n GFDL SSP5-8.5")
dev.copy(jpeg,
         here('results/cod_forecast',
              'cod_larvae_gfdl_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
df_pcodlarvae3_gfdl585 <- predict_cells(2070:2099, pcod_larvae, 134,
                                       5, 'ssp585', gfdl_temps3, 
                                       gfdl_salts3, pcodlarvae_formula)
saveRDS(df_pcodlarvae3_gfdl585, file = here("data", "df_pcodlarvae3_gfdl585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pcodlarvae3_gfdl585, 
             "Forecasted Distribution 2070 - 2099 \n GFDL SSP5-8.5")
dev.copy(jpeg,
         here('results/cod_forecast',
              'cod_larvae_gfdl_ssp585_3.jpg'),
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

df_pcodlarvae1_miroc126 <- predict_cells(2015:2039, pcod_larvae, 134,
                                        5, 'ssp126', miroc_temps1, 
                                        miroc_salts1, pcodlarvae_formula)
saveRDS(df_pcodlarvae1_miroc126, file = here("data", "df_pcodlarvae1_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pcodlarvae1_miroc126, 
             "Forecasted Distribution 2015 - 2039 \n MIROC SSP1-2.6")
dev.copy(jpeg,
         here('results/cod_forecast',
              'cod_larvae_miroc_ssp126_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
miroc_temps2 <- readRDS(here('data', 'miroc_forecast_temp2.rds'))
miroc_salts2 <- readRDS(here('data', 'miroc_forecast_salt2.rds'))

df_pcodlarvae2_miroc126 <- predict_cells(2040:2069, pcod_larvae, 134,
                                        5, 'ssp126', miroc_temps2, 
                                        miroc_salts2, pcodlarvae_formula)
saveRDS(df_pcodlarvae2_miroc126, file = here("data", "df_pcodlarvae2_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pcodlarvae2_miroc126, 
             "Forecasted Distribution 2040 - 2069 \n MIROC SSP1-2.6")
dev.copy(jpeg,
         here('results/cod_forecast',
              'cod_larvae_miroc_ssp126_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
miroc_temps3 <- readRDS(here('data', 'miroc_forecast_temp3.rds'))
miroc_salts3 <- readRDS(here('data', 'miroc_forecast_salt3.rds'))

df_pcodlarvae3_miroc126 <- predict_cells(2070:2099, pcod_larvae, 134,
                                        5, 'ssp126', miroc_temps3, 
                                        miroc_salts3, pcodlarvae_formula)
saveRDS(df_pcodlarvae3_miroc126, file = here("data", "df_pcodlarvae3_miroc126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pcodlarvae3_miroc126, 
             "Forecasted Distribution 2070 - 2099 \n MIROC SSP1-2.6")
dev.copy(jpeg,
         here('results/cod_forecast',
              'cod_larvae_miroc_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


##### MIROC 585 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
df_pcodlarvae1_miroc585 <- predict_cells(2015:2039, pcod_larvae, 134,
                                        5, 'ssp585', miroc_temps1, 
                                        miroc_salts1, pcodlarvae_formula)
saveRDS(df_pcodlarvae1_miroc585, file = here("data", "df_pcodlarvae1_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pcodlarvae1_miroc585, 
             "Forecasted Distribution 2015 - 2039 \n MIROC SSP5-8.5")
dev.copy(jpeg,
         here('results/cod_forecast',
              'cod_larvae_miroc_ssp585_1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2040 - 2069
df_pcodlarvae2_miroc585 <- predict_cells(2040:2069, pcod_larvae, 134,
                                        5, 'ssp585', miroc_temps2, 
                                        miroc_salts2, pcodlarvae_formula)
saveRDS(df_pcodlarvae2_miroc585, file = here("data", "df_pcodlarvae2_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pcodlarvae2_miroc585, 
             "Forecasted Distribution 2040 - 2069 \n MIROC SSP5-8.5")
dev.copy(jpeg,
         here('results/cod_forecast',
              'cod_larvae_miroc_ssp585_2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## 2070 - 2099
df_pcodlarvae3_miroc585 <- predict_cells(2070:2099, pcod_larvae, 134,
                                        5, 'ssp585', miroc_temps3,
                                        miroc_salts3, pcodlarvae_formula)
saveRDS(df_pcodlarvae3_miroc585, file = here("data", "df_pcodlarvae3_miroc585.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pcodlarvae3_miroc585, 
             "Forecasted Distribution 2070 - 2099 \n MIROC SSP5-8.5")
dev.copy(jpeg,
         here('results/cod_forecast',
              'cod_larvae_miroc_ssp585_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

rm(miroc_temps1, miroc_temps2, miroc_temps3,
   miroc_salts1, miroc_salts2, miroc_salts3)

graphics.off()

##### Multi-panel figure -------------------------------------------------------------------------------------------------------------------------
tiff(here('results/cod_forecast',
          'codlarvae_multipanel.tiff'),
     units = "in",
     width = 45,
     height = 90,
     res = 300)
par(mfrow = c(6, 3),
    mar = c(11, 12, 5, 1.5) + 0.1,
    oma = c(3, 25, 15, 1),
    mgp = c(10, 4, 0),
    family = "serif")
grid_multipanel(df_pcodlarvae1_cesm126)
grid_multipanel(df_pcodlarvae2_cesm126)
grid_multipanel(df_pcodlarvae3_cesm126)
grid_multipanel(df_pcodlarvae1_cesm585)
grid_multipanel(df_pcodlarvae2_cesm585)
grid_multipanel(df_pcodlarvae3_cesm585)
grid_multipanel(df_pcodlarvae1_gfdl126)
grid_multipanel(df_pcodlarvae2_gfdl126)
grid_multipanel(df_pcodlarvae3_gfdl126)
grid_multipanel(df_pcodlarvae1_gfdl585)
grid_multipanel(df_pcodlarvae2_gfdl585)
grid_multipanel(df_pcodlarvae3_gfdl585)
grid_multipanel(df_pcodlarvae1_miroc126)
grid_multipanel(df_pcodlarvae2_miroc126)
grid_multipanel(df_pcodlarvae3_miroc126)
grid_multipanel(df_pcodlarvae1_miroc585)
grid_multipanel(df_pcodlarvae2_miroc585)
grid_multipanel(df_pcodlarvae3_miroc585)
mtext("CESM SSP1-2.6", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.92)
mtext("CESM SSP5-8.5", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.76)
mtext("GFDL1-2.6", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.585)
mtext("GFDL5-8.5", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.42)
mtext("MIROC1-2.6", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.26)
mtext("MIROC5-8.5", 
      side = 2, 
      line = 12, 
      outer = TRUE, 
      cex = 7,
      at = 0.095)
mtext("2015-2039", 
      side = 3, 
      line = 2, 
      outer = TRUE, 
      cex = 7,
      at = 0.17)
mtext("2040-2069", 
      side = 3, 
      line = 2, 
      outer = TRUE, 
      cex = 7,
      at = 0.51)
mtext("2070-2099", 
      side = 3, 
      line = 2, 
      outer = TRUE, 
      cex = 7,
      at = 0.84)
dev.off()

##### Averages ---------------------------------------------------------------------------------------------------------------------------
# 2015 - 2039
df_pcodlarvae_merged1 <- list(df_pcodlarvae1_cesm126, df_pcodlarvae1_cesm585,
                             df_pcodlarvae1_gfdl126, df_pcodlarvae1_gfdl585,
                             df_pcodlarvae1_miroc126, df_pcodlarvae1_miroc585) 

avg_pcodlarvae_merged1 <- predict_avgs(df_pcodlarvae_merged1)

windows(width = 6, height = 6, family = "serif")
grid_predict(avg_pcodlarvae_merged1, "Forecasted Distribution 2015 - 2039")
dev.copy(jpeg,
         here('results/cod_forecast/pcodlarvae_avgs',
              'cod_larvae_avg1.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

# 2040-2069
df_pcodlarvae_merged2 <- list(df_pcodlarvae2_cesm126, df_pcodlarvae2_cesm585,
                             df_pcodlarvae2_gfdl126, df_pcodlarvae2_gfdl585,
                             df_pcodlarvae2_miroc126, df_pcodlarvae2_miroc585) 

avg_pcodlarvae_merged2 <- predict_avgs(df_pcodlarvae_merged2)

windows(width = 6, height = 6, family = "serif")
grid_predict(avg_pcodlarvae_merged2, "Forecasted Distribution 2040 - 2069")
dev.copy(jpeg,
         here('results/cod_forecast/pcodlarvae_avgs',
              'cod_larvae_avg2.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

# 2070-2099
df_pcodlarvae_merged3 <- list(df_pcodlarvae3_cesm126, df_pcodlarvae3_cesm585,
                             df_pcodlarvae3_gfdl126, df_pcodlarvae3_gfdl585,
                             df_pcodlarvae3_miroc126, df_pcodlarvae3_miroc585) 

avg_pcodlarvae_merged3 <- predict_avgs(df_pcodlarvae_merged3)

windows(width = 6, height = 6, family = "serif")
grid_predict(avg_pcodlarvae_merged3, "Forecasted Distribution 2070 - 2099")
dev.copy(jpeg,
         here('results/cod_forecast/pcodlarvae_avgs',
              'cod_larvae_avg3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

##### GIFs -------------------------------------------------------------------------------------------------------------------------
pcodlarvae_dir_out <- file.path(base_dir, 'results', 'cod_forecast', 'pcodlarvae_avgs')
pcodlarvae_imgs <- list.files(pcodlarvae_dir_out, full.names = T)
pcodlarvae_img_list <- lapply(pcodlarvae_imgs, image_read)
pcodlarvae_img_joined <- image_join(pcodlarvae_img_list)
pcodlarvae_img_animated <- image_animate(pcodlarvae_img_joined, fps = 1)
image_write(image = pcodlarvae_img_animated,
            path = here('results', 'cod_forecast', "pcodlarvae_avgs.gif"))

##### Clear environment -------------------------------------------------------------------------------------------------------------------------
rm(df_pcodlarvae1_cesm126, df_pcodlarvae1_cesm585,
   df_pcodlarvae1_gfdl126, df_pcodlarvae1_gfdl585,
   df_pcodlarvae1_miroc126, df_pcodlarvae1_miroc585,
   df_pcodlarvae2_cesm126, df_pcodlarvae2_cesm585,
   df_pcodlarvae2_gfdl126, df_pcodlarvae2_gfdl585,
   df_pcodlarvae2_miroc126, df_pcodlarvae2_miroc585,
   df_pcodlarvae3_cesm126, df_pcodlarvae3_cesm585,
   df_pcodlarvae3_gfdl126, df_pcodlarvae3_gfdl585,
   df_pcodlarvae3_miroc126, df_pcodlarvae3_miroc585,
   avg_pcodlarvae_merged1, avg_pcodlarvae_merged2,
   avg_pcodlarvae_merged3, df_pcodlarvae_merged1,
   df_pcodlarvae_merged2, df_pcodlarvae_merged3,
   pcod_larvae, pcodlarvae_formula, pcodlarvae_img_animated,
   pcodlarvae_dir_out, pcodlarvae_img_joined, 
   pcodlarvae_img_list, pcodlarvae_imgs)