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
      family = tw(link = 'log'))
}

formula_geog <- function(data){
  gam(catch ~ s(year, bs = 're') +
        s(doy, k = 8) +
        s(lon, lat) +
        s(roms_temperature, k = 6) +
        s(roms_salinity, k = 6) +
        s(lat, lon, by = mean_temp), # geography
      data = data,
      family = tw(link = 'log'))
}

# Load ROMS temperature means and forecast
roms_temps <- readRDS(here('data', 'roms_temps.rds'))

# Load fish data
pk_egg <- load_data('pk_egg.rds', pk_egg, roms_temps)
pk_larvae <- load_data('pk_larvae.rds', pk_larvae, roms_temps)
fhs_egg <- load_data('fhs_egg.rds', fhs_egg, roms_temps)
fhs_larvae <- load_data('fhs_larvae.rds', fhs_larvae, roms_temps)
akp_egg <- load_data('akp_egg.rds', akp_egg, roms_temps)
yfs_larvae <- load_data('yfs_larvae.rds', yfs_larvae, roms_temps)
pcod_larvae <- load_data('pcod_larvae.rds', pcod_larvae, roms_temps)
nrs_larvae <- load_data('nrs_larvae.rds', nrs_larvae, roms_temps)

### Pollock Eggs --------------------------------------------------------------------------------------------------------------------------
#### Forecast and average into 3 time periods ---------------------------------------------------------------------------------------------
##### CESM 126 ----------------------------------------------------------------------------------------------------------------------------
## 2015 - 2039
cesm_temps1 <- readRDS(here('data', 'cesm_forecast_temp1.rds'))
cesm_salts1 <- readRDS(here('data', 'cesm_forecast_salt1.rds'))
pkegg_formula <- formula_pheno(pk_egg)

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

preds_pkegg2_cesm126 <- pred_loop(2040:2069, pk_egg, 130,
                                  5, 'ssp126', cesm_temps2, 
                                  cesm_salts2, pkegg_formula)

# Combine into one data frame
df_pkegg2_cesm126 <- list(preds_pkegg2_cesm126[[1]], preds_pkegg2_cesm126[[2]],
                          preds_pkegg2_cesm126[[3]], preds_pkegg2_cesm126[[4]],
                          preds_pkegg2_cesm126[[5]], preds_pkegg2_cesm126[[6]],
                          preds_pkegg2_cesm126[[7]], preds_pkegg2_cesm126[[8]],
                          preds_pkegg2_cesm126[[9]], preds_pkegg2_cesm126[[10]],
                          preds_pkegg2_cesm126[[11]], preds_pkegg2_cesm126[[12]],
                          preds_pkegg2_cesm126[[13]], preds_pkegg2_cesm126[[14]],
                          preds_pkegg2_cesm126[[15]], preds_pkegg2_cesm126[[16]],
                          preds_pkegg2_cesm126[[17]], preds_pkegg2_cesm126[[18]],
                          preds_pkegg2_cesm126[[19]], preds_pkegg2_cesm126[[20]],
                          preds_pkegg2_cesm126[[21]], preds_pkegg2_cesm126[[22]],
                          preds_pkegg2_cesm126[[23]], preds_pkegg2_cesm126[[24]],
                          preds_pkegg2_cesm126[[25]], preds_pkegg2_cesm126[[26]],
                          preds_pkegg2_cesm126[[27]], preds_pkegg2_cesm126[[28]],
                          preds_pkegg2_cesm126[[29]], preds_pkegg2_cesm126[[30]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 

# Calculate mean predicted abundance per year
preds_pkegg2_cesm126_avgs <- sapply(preds_pkegg2_cesm126, function(x) colMeans(select(x, pred), na.rm = T))
df_pkegg2_cesm126_avgs <- data.frame(year = c(2040:2069), 
                                     avg_pred = preds_pkegg2_cesm126_avgs)

ggplot(df_pkegg2_cesm126_avgs) +
  geom_line(aes(x = year,
                y = avg_pred))


# Generate average prediction from all predictions
x <- grepl("pred", names(df_pkegg2_cesm126), fixed = T)
df_pkegg_avg2_cesm126 <- data.frame(lat = df_pkegg2_cesm126$lat, 
                                    lon = df_pkegg2_cesm126$lon, 
                                    dist = df_pkegg2_cesm126$dist,
                                    avg_pred = rowSums(df_pkegg2_cesm126[, x])/30)
saveRDS(df_pkegg_avg2_cesm126, file = here("data", "df_pkegg_avg2_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pkegg_avg2_cesm126, "Forecasted Distribution 2040 - 2069 \n CESM SSP126")
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

preds_pkegg3_cesm126 <- pred_loop(2070:2099, pk_egg, 130,
                                  5, 'ssp126', cesm_temps3, 
                                  cesm_salts3, egg_formula)

# Combine into one data frame
df_pkegg3_cesm126 <- list(preds_pkegg3_cesm126[[1]], preds_pkegg3_cesm126[[2]],
                          preds_pkegg3_cesm126[[3]], preds_pkegg3_cesm126[[4]],
                          preds_pkegg3_cesm126[[5]], preds_pkegg3_cesm126[[6]],
                          preds_pkegg3_cesm126[[7]], preds_pkegg3_cesm126[[8]],
                          preds_pkegg3_cesm126[[9]], preds_pkegg3_cesm126[[10]],
                          preds_pkegg3_cesm126[[11]], preds_pkegg3_cesm126[[12]],
                          preds_pkegg3_cesm126[[13]], preds_pkegg3_cesm126[[14]],
                          preds_pkegg3_cesm126[[15]], preds_pkegg3_cesm126[[16]],
                          preds_pkegg3_cesm126[[17]], preds_pkegg3_cesm126[[18]],
                          preds_pkegg3_cesm126[[19]], preds_pkegg3_cesm126[[20]],
                          preds_pkegg3_cesm126[[21]], preds_pkegg3_cesm126[[22]],
                          preds_pkegg3_cesm126[[23]], preds_pkegg3_cesm126[[24]],
                          preds_pkegg3_cesm126[[25]], preds_pkegg3_cesm126[[26]], 
                          preds_pkegg3_cesm126[[27]], preds_pkegg3_cesm126[[28]], 
                          preds_pkegg3_cesm126[[29]], preds_pkegg3_cesm126[[30]]) %>%
  reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 

# Calculate mean predicted abundance per year
preds_pkegg3_cesm126_avgs <- sapply(preds_pkegg3_cesm126, function(x) colMeans(select(x, pred), na.rm = T))
df_pkegg3_cesm126_avgs <- data.frame(year = c(2070:2099), 
                                     avg_pred = preds_pkegg3_cesm126_avgs)

ggplot(df_pkegg3_cesm126_avgs) +
  geom_line(aes(x = year,
                y = avg_pred))


# Generate average prediction from all predictions
x <- grepl("pred", names(df_pkegg3_cesm126), fixed = T)
df_pkegg_avg3_cesm126 <- data.frame(lat = df_pkegg3_cesm126$lat, 
                                    lon = df_pkegg3_cesm126$lon, 
                                    dist = df_pkegg3_cesm126$dist,
                                    avg_pred = rowSums(df_pkegg3_cesm126[, x])/30)
saveRDS(df_pkegg_avg3_cesm126, file = here("data", "df_pkegg_avg3_cesm126.rds"))

# Plot
windows(width = 6, height = 6, family = "serif")
grid_predict(df_pkegg_avg3_cesm126, "Forecasted Distribution 2070 - 2099 \n CESM SSP126")
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_egg_cesm_ssp126_3.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


rm(cesm_temps1, cesm_temps2, cesm_temps3,
   cesm_salts1, cesm_salts2, cesm_salts3)