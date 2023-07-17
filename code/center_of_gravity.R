# Title: Calculate center of gravity
# Authors: Rebecca Howard
# Date: 07/17/2023

## Libraries ----
library(dplyr)
library(here)

## Functions ----
# Calculate center of gravity per time period
COG <- function(data){
  temp1 <- data %>% dplyr::mutate(pred = pred_scaled * lat) %>%
    dplyr::summarize(value = sum(pred, na.rm = TRUE)) %>%
    dplyr::select(value)
  
  temp2 <- data %>%
    dplyr::summarize(value = sum(pred_scaled, na.rm = TRUE)) %>%
    dplyr::select(value)
  
  COG_y <- temp1 / temp2
  
  temp3 <- data %>% dplyr::mutate(pred = pred_scaled * lon) %>%
    dplyr::summarize(value1 = sum(pred, na.rm = TRUE)) %>%
    dplyr::select(value1)
  
  temp4 <- data %>%
    dplyr::summarize(value1 = sum(pred_scaled, na.rm = TRUE)) %>%
    dplyr::select(value1)
  
  COG_x <- temp3 / temp4
  
  COG <- as.data.frame(cbind(COG_y, COG_x))
}

COG_loop <- function(data, the_list){
  for (i in 1:length(data)){
    the_list[[i]] <- COG(data[[i]])
  }
  print(the_list)
}

# Make loop to calculate per time frame
COG_time <- function(hindcast, list1, list2, list3){
  time1_COG <- list()
  time2_COG <- list()
  time3_COG <- list()
  
  hindcast_COG <- COG(hindcast)
  time1_COG <- COG_loop(list1, time1_COG)
  time2_COG <- COG_loop(list2, time2_COG)
  time3_COG <- COG_loop(list3, time3_COG)

  avg1_lat <- rowMeans(do.call(cbind, lapply(time1_COG, "[", "value")))
  avg1_lon <- rowMeans(do.call(cbind, lapply(time1_COG, "[", "value1")))
  
  avg2_lat <- rowMeans(do.call(cbind, lapply(time2_COG, "[", "value")))
  avg2_lon <- rowMeans(do.call(cbind, lapply(time2_COG, "[", "value1")))
  
  avg3_lat <- rowMeans(do.call(cbind, lapply(time3_COG, "[", "value")))
  avg3_lon <- rowMeans(do.call(cbind, lapply(time3_COG, "[", "value1")))
  
  lats <- rbind(hindcast_COG$value, avg1_lat, avg2_lat, avg3_lat)
  lons <- rbind(hindcast_COG$value1, avg1_lon, avg2_lon, avg3_lon)
  
  avg_df <- cbind(lats, lons)
  colnames(avg_df) <- c("lat", "lon")
  rownames(avg_df) <- c("hindcast", "time1", "time2", "time3")
  
  return(avg_df)
}

## Pollock ----
### Eggs ----
# Load data
pkegg_hindcast <- readRDS(here('data', 'pk_egg_hindcast'))
df_pkegg1_cesm126 <- readRDS(here('data', 'df_pkegg1_cesm126.rds'))
df_pkegg2_cesm126 <- readRDS(here('data', 'df_pkegg2_cesm126.rds'))
df_pkegg3_cesm126 <- readRDS(here('data', 'df_pkegg3_cesm126.rds'))
df_pkegg1_cesm585 <- readRDS(here('data', 'df_pkegg1_cesm585.rds'))
df_pkegg2_cesm585 <- readRDS(here('data', 'df_pkegg2_cesm585.rds'))
df_pkegg3_cesm585 <- readRDS(here('data', 'df_pkegg3_cesm585.rds'))
df_pkegg1_gfdl126 <- readRDS(here('data', 'df_pkegg1_gfdl126.rds'))
df_pkegg2_gfdl126 <- readRDS(here('data', 'df_pkegg2_gfdl126.rds'))
df_pkegg3_gfdl126 <- readRDS(here('data', 'df_pkegg3_gfdl126.rds'))
df_pkegg1_gfdl585 <- readRDS(here('data', 'df_pkegg1_gfdl585.rds'))
df_pkegg2_gfdl585 <- readRDS(here('data', 'df_pkegg2_gfdl585.rds'))
df_pkegg3_gfdl585 <- readRDS(here('data', 'df_pkegg3_gfdl585.rds'))
df_pkegg1_miroc126 <- readRDS(here('data', 'df_pkegg1_miroc126.rds'))
df_pkegg2_miroc126 <- readRDS(here('data', 'df_pkegg2_miroc126.rds'))
df_pkegg3_miroc126 <- readRDS(here('data', 'df_pkegg3_miroc126.rds'))
df_pkegg1_miroc585 <- readRDS(here('data', 'df_pkegg1_miroc585.rds'))
df_pkegg2_miroc585 <- readRDS(here('data', 'df_pkegg2_miroc585.rds'))
df_pkegg3_miroc585 <- readRDS(here('data', 'df_pkegg3_miroc585.rds'))

pkegg1_list <- list(df_pkegg1_cesm126, df_pkegg1_cesm585,
                   df_pkegg1_gfdl126, df_pkegg1_gfdl585,
                   df_pkegg1_miroc126, df_pkegg1_miroc585)

pkegg2_list <- list(df_pkegg2_cesm126, df_pkegg2_cesm585,
                   df_pkegg2_gfdl126, df_pkegg2_gfdl585,
                   df_pkegg2_miroc126, df_pkegg2_miroc585)

pkegg3_list <- list(df_pkegg3_cesm126, df_pkegg3_cesm585,
                   df_pkegg3_gfdl126, df_pkegg3_gfdl585,
                   df_pkegg3_miroc126, df_pkegg3_miroc585)

pkegg_COG <- COG_time(pkegg_hindcast, pkegg1_list, pkegg2_list, pkegg3_list)
saveRDS(pkegg_COG, here('data', 'pkegg_COG'))

rm(df_pkegg1_cesm126, df_pkegg1_cesm585, df_pkegg1_gfdl126, df_pkegg1_gfdl585,
   df_pkegg1_miroc126, df_pkegg1_miroc585, df_pkegg2_cesm126, df_pkegg2_cesm585, 
   df_pkegg2_gfdl126, df_pkegg2_gfdl585, df_pkegg2_miroc126, df_pkegg2_miroc585,
   df_pkegg3_cesm126, df_pkegg3_cesm585, df_pkegg3_gfdl126, df_pkegg3_gfdl585,
   df_pkegg3_miroc126, df_pkegg3_miroc585, pkegg1_list, pkegg2_list, pkegg3_list)

### Larvae ----
# Load data
df_pklarvae1_cesm126 <- readRDS(here('data', 'df_pklarvae1_cesm126.rds'))
df_pklarvae2_cesm126 <- readRDS(here('data', 'df_pklarvae2_cesm126.rds'))
df_pklarvae3_cesm126 <- readRDS(here('data', 'df_pklarvae3_cesm126.rds'))
df_pklarvae1_cesm585 <- readRDS(here('data', 'df_pklarvae1_cesm585.rds'))
df_pklarvae2_cesm585 <- readRDS(here('data', 'df_pklarvae2_cesm585.rds'))
df_pklarvae3_cesm585 <- readRDS(here('data', 'df_pklarvae3_cesm585.rds'))
df_pklarvae1_gfdl126 <- readRDS(here('data', 'df_pklarvae1_gfdl126.rds'))
df_pklarvae2_gfdl126 <- readRDS(here('data', 'df_pklarvae2_gfdl126.rds'))
df_pklarvae3_gfdl126 <- readRDS(here('data', 'df_pklarvae3_gfdl126.rds'))
df_pklarvae1_gfdl585 <- readRDS(here('data', 'df_pklarvae1_gfdl585.rds'))
df_pklarvae2_gfdl585 <- readRDS(here('data', 'df_pklarvae2_gfdl585.rds'))
df_pklarvae3_gfdl585 <- readRDS(here('data', 'df_pklarvae3_gfdl585.rds'))
df_pklarvae1_miroc126 <- readRDS(here('data', 'df_pklarvae1_miroc126.rds'))
df_pklarvae2_miroc126 <- readRDS(here('data', 'df_pklarvae2_miroc126.rds'))
df_pklarvae3_miroc126 <- readRDS(here('data', 'df_pklarvae3_miroc126.rds'))
df_pklarvae1_miroc585 <- readRDS(here('data', 'df_pklarvae1_miroc585.rds'))
df_pklarvae2_miroc585 <- readRDS(here('data', 'df_pklarvae2_miroc585.rds'))
df_pklarvae3_miroc585 <- readRDS(here('data', 'df_pklarvae3_miroc585.rds'))

pklarvae1_list <- list(df_pklarvae1_cesm126, df_pklarvae1_cesm585,
                    df_pklarvae1_gfdl126, df_pklarvae1_gfdl585,
                    df_pklarvae1_miroc126, df_pklarvae1_miroc585)

pklarvae2_list <- list(df_pklarvae2_cesm126, df_pklarvae2_cesm585,
                    df_pklarvae2_gfdl126, df_pklarvae2_gfdl585,
                    df_pklarvae2_miroc126, df_pklarvae2_miroc585)

pklarvae3_list <- list(df_pklarvae3_cesm126, df_pklarvae3_cesm585,
                    df_pklarvae3_gfdl126, df_pklarvae3_gfdl585,
                    df_pklarvae3_miroc126, df_pklarvae3_miroc585)

pklarvae_COG <- COG_time(pklarvae1_list, pklarvae2_list, pklarvae3_list)
saveRDS(pklarvae_COG, here('data', 'pklarvae_COG'))

rm(df_pklarvae1_cesm126, df_pklarvae1_cesm585, df_pklarvae1_gfdl126, df_pklarvae1_gfdl585,
   df_pklarvae1_miroc126, df_pklarvae1_miroc585, df_pklarvae2_cesm126, df_pklarvae2_cesm585, 
   df_pklarvae2_gfdl126, df_pklarvae2_gfdl585, df_pklarvae2_miroc126, df_pklarvae2_miroc585,
   df_pklarvae3_cesm126, df_pklarvae3_cesm585, df_pklarvae3_gfdl126, df_pklarvae3_gfdl585,
   df_pklarvae3_miroc126, df_pklarvae3_miroc585, pklarvae1_list, pklarvae2_list, pklarvae3_list)

## Flathead Sole ----
### Eggs ----
# Load data
df_fhsegg1_cesm126 <- readRDS(here('data', 'df_fhsegg1_cesm126.rds'))
df_fhsegg2_cesm126 <- readRDS(here('data', 'df_fhsegg2_cesm126.rds'))
df_fhsegg3_cesm126 <- readRDS(here('data', 'df_fhsegg3_cesm126.rds'))
df_fhsegg1_cesm585 <- readRDS(here('data', 'df_fhsegg1_cesm585.rds'))
df_fhsegg2_cesm585 <- readRDS(here('data', 'df_fhsegg2_cesm585.rds'))
df_fhsegg3_cesm585 <- readRDS(here('data', 'df_fhsegg3_cesm585.rds'))
df_fhsegg1_gfdl126 <- readRDS(here('data', 'df_fhsegg1_gfdl126.rds'))
df_fhsegg2_gfdl126 <- readRDS(here('data', 'df_fhsegg2_gfdl126.rds'))
df_fhsegg3_gfdl126 <- readRDS(here('data', 'df_fhsegg3_gfdl126.rds'))
df_fhsegg1_gfdl585 <- readRDS(here('data', 'df_fhsegg1_gfdl585.rds'))
df_fhsegg2_gfdl585 <- readRDS(here('data', 'df_fhsegg2_gfdl585.rds'))
df_fhsegg3_gfdl585 <- readRDS(here('data', 'df_fhsegg3_gfdl585.rds'))
df_fhsegg1_miroc126 <- readRDS(here('data', 'df_fhsegg1_miroc126.rds'))
df_fhsegg2_miroc126 <- readRDS(here('data', 'df_fhsegg2_miroc126.rds'))
df_fhsegg3_miroc126 <- readRDS(here('data', 'df_fhsegg3_miroc126.rds'))
df_fhsegg1_miroc585 <- readRDS(here('data', 'df_fhsegg1_miroc585.rds'))
df_fhsegg2_miroc585 <- readRDS(here('data', 'df_fhsegg2_miroc585.rds'))
df_fhsegg3_miroc585 <- readRDS(here('data', 'df_fhsegg3_miroc585.rds'))

fhsegg1_list <- list(df_fhsegg1_cesm126, df_fhsegg1_cesm585,
                    df_fhsegg1_gfdl126, df_fhsegg1_gfdl585,
                    df_fhsegg1_miroc126, df_fhsegg1_miroc585)

fhsegg2_list <- list(df_fhsegg2_cesm126, df_fhsegg2_cesm585,
                    df_fhsegg2_gfdl126, df_fhsegg2_gfdl585,
                    df_fhsegg2_miroc126, df_fhsegg2_miroc585)

fhsegg3_list <- list(df_fhsegg3_cesm126, df_fhsegg3_cesm585,
                    df_fhsegg3_gfdl126, df_fhsegg3_gfdl585,
                    df_fhsegg3_miroc126, df_fhsegg3_miroc585)

fhsegg_COG <- COG_time(fhsegg1_list, fhsegg2_list, fhsegg3_list)
saveRDS(fhsegg_COG, here('data', 'fhsegg_COG'))

rm(df_fhsegg1_cesm126, df_fhsegg1_cesm585, df_fhsegg1_gfdl126, df_fhsegg1_gfdl585,
   df_fhsegg1_miroc126, df_fhsegg1_miroc585, df_fhsegg2_cesm126, df_fhsegg2_cesm585, 
   df_fhsegg2_gfdl126, df_fhsegg2_gfdl585, df_fhsegg2_miroc126, df_fhsegg2_miroc585,
   df_fhsegg3_cesm126, df_fhsegg3_cesm585, df_fhsegg3_gfdl126, df_fhsegg3_gfdl585,
   df_fhsegg3_miroc126, df_fhsegg3_miroc585, fhsegg1_list, fhsegg2_list, fhsegg3_list)

### Larvae ----
# Load data
df_fhslarvae1_cesm126 <- readRDS(here('data', 'df_fhslarvae1_cesm126.rds'))
df_fhslarvae2_cesm126 <- readRDS(here('data', 'df_fhslarvae2_cesm126.rds'))
df_fhslarvae3_cesm126 <- readRDS(here('data', 'df_fhslarvae3_cesm126.rds'))
df_fhslarvae1_cesm585 <- readRDS(here('data', 'df_fhslarvae1_cesm585.rds'))
df_fhslarvae2_cesm585 <- readRDS(here('data', 'df_fhslarvae2_cesm585.rds'))
df_fhslarvae3_cesm585 <- readRDS(here('data', 'df_fhslarvae3_cesm585.rds'))
df_fhslarvae1_gfdl126 <- readRDS(here('data', 'df_fhslarvae1_gfdl126.rds'))
df_fhslarvae2_gfdl126 <- readRDS(here('data', 'df_fhslarvae2_gfdl126.rds'))
df_fhslarvae3_gfdl126 <- readRDS(here('data', 'df_fhslarvae3_gfdl126.rds'))
df_fhslarvae1_gfdl585 <- readRDS(here('data', 'df_fhslarvae1_gfdl585.rds'))
df_fhslarvae2_gfdl585 <- readRDS(here('data', 'df_fhslarvae2_gfdl585.rds'))
df_fhslarvae3_gfdl585 <- readRDS(here('data', 'df_fhslarvae3_gfdl585.rds'))
df_fhslarvae1_miroc126 <- readRDS(here('data', 'df_fhslarvae1_miroc126.rds'))
df_fhslarvae2_miroc126 <- readRDS(here('data', 'df_fhslarvae2_miroc126.rds'))
df_fhslarvae3_miroc126 <- readRDS(here('data', 'df_fhslarvae3_miroc126.rds'))
df_fhslarvae1_miroc585 <- readRDS(here('data', 'df_fhslarvae1_miroc585.rds'))
df_fhslarvae2_miroc585 <- readRDS(here('data', 'df_fhslarvae2_miroc585.rds'))
df_fhslarvae3_miroc585 <- readRDS(here('data', 'df_fhslarvae3_miroc585.rds'))

fhslarvae1_list <- list(df_fhslarvae1_cesm126, df_fhslarvae1_cesm585,
                       df_fhslarvae1_gfdl126, df_fhslarvae1_gfdl585,
                       df_fhslarvae1_miroc126, df_fhslarvae1_miroc585)

fhslarvae2_list <- list(df_fhslarvae2_cesm126, df_fhslarvae2_cesm585,
                       df_fhslarvae2_gfdl126, df_fhslarvae2_gfdl585,
                       df_fhslarvae2_miroc126, df_fhslarvae2_miroc585)

fhslarvae3_list <- list(df_fhslarvae3_cesm126, df_fhslarvae3_cesm585,
                       df_fhslarvae3_gfdl126, df_fhslarvae3_gfdl585,
                       df_fhslarvae3_miroc126, df_fhslarvae3_miroc585)

fhslarvae_COG <- COG_time(fhslarvae1_list, fhslarvae2_list, fhslarvae3_list)
saveRDS(fhslarvae_COG, here('data', 'fhslarvae_COG'))

rm(df_fhslarvae1_cesm126, df_fhslarvae1_cesm585, df_fhslarvae1_gfdl126, df_fhslarvae1_gfdl585,
   df_fhslarvae1_miroc126, df_fhslarvae1_miroc585, df_fhslarvae2_cesm126, df_fhslarvae2_cesm585, 
   df_fhslarvae2_gfdl126, df_fhslarvae2_gfdl585, df_fhslarvae2_miroc126, df_fhslarvae2_miroc585,
   df_fhslarvae3_cesm126, df_fhslarvae3_cesm585, df_fhslarvae3_gfdl126, df_fhslarvae3_gfdl585,
   df_fhslarvae3_miroc126, df_fhslarvae3_miroc585, fhslarvae1_list, fhslarvae2_list, fhslarvae3_list)

## Alaska Plaice ----
### Eggs ----
# Load data
df_akpegg1_cesm126 <- readRDS(here('data', 'df_akpegg1_cesm126.rds'))
df_akpegg2_cesm126 <- readRDS(here('data', 'df_akpegg2_cesm126.rds'))
df_akpegg3_cesm126 <- readRDS(here('data', 'df_akpegg3_cesm126.rds'))
df_akpegg1_cesm585 <- readRDS(here('data', 'df_akpegg1_cesm585.rds'))
df_akpegg2_cesm585 <- readRDS(here('data', 'df_akpegg2_cesm585.rds'))
df_akpegg3_cesm585 <- readRDS(here('data', 'df_akpegg3_cesm585.rds'))
df_akpegg1_gfdl126 <- readRDS(here('data', 'df_akpegg1_gfdl126.rds'))
df_akpegg2_gfdl126 <- readRDS(here('data', 'df_akpegg2_gfdl126.rds'))
df_akpegg3_gfdl126 <- readRDS(here('data', 'df_akpegg3_gfdl126.rds'))
df_akpegg1_gfdl585 <- readRDS(here('data', 'df_akpegg1_gfdl585.rds'))
df_akpegg2_gfdl585 <- readRDS(here('data', 'df_akpegg2_gfdl585.rds'))
df_akpegg3_gfdl585 <- readRDS(here('data', 'df_akpegg3_gfdl585.rds'))
df_akpegg1_miroc126 <- readRDS(here('data', 'df_akpegg1_miroc126.rds'))
df_akpegg2_miroc126 <- readRDS(here('data', 'df_akpegg2_miroc126.rds'))
df_akpegg3_miroc126 <- readRDS(here('data', 'df_akpegg3_miroc126.rds'))
df_akpegg1_miroc585 <- readRDS(here('data', 'df_akpegg1_miroc585.rds'))
df_akpegg2_miroc585 <- readRDS(here('data', 'df_akpegg2_miroc585.rds'))
df_akpegg3_miroc585 <- readRDS(here('data', 'df_akpegg3_miroc585.rds'))

akpegg1_list <- list(df_akpegg1_cesm126, df_akpegg1_cesm585,
                     df_akpegg1_gfdl126, df_akpegg1_gfdl585,
                     df_akpegg1_miroc126, df_akpegg1_miroc585)

akpegg2_list <- list(df_akpegg2_cesm126, df_akpegg2_cesm585,
                     df_akpegg2_gfdl126, df_akpegg2_gfdl585,
                     df_akpegg2_miroc126, df_akpegg2_miroc585)

akpegg3_list <- list(df_akpegg3_cesm126, df_akpegg3_cesm585,
                     df_akpegg3_gfdl126, df_akpegg3_gfdl585,
                     df_akpegg3_miroc126, df_akpegg3_miroc585)

akpegg_COG <- COG_time(akpegg1_list, akpegg2_list, akpegg3_list)
saveRDS(akpegg_COG, here('data', 'akpegg_COG'))

rm(df_akpegg1_cesm126, df_akpegg1_cesm585, df_akpegg1_gfdl126, df_akpegg1_gfdl585,
   df_akpegg1_miroc126, df_akpegg1_miroc585, df_akpegg2_cesm126, df_akpegg2_cesm585, 
   df_akpegg2_gfdl126, df_akpegg2_gfdl585, df_akpegg2_miroc126, df_akpegg2_miroc585,
   df_akpegg3_cesm126, df_akpegg3_cesm585, df_akpegg3_gfdl126, df_akpegg3_gfdl585,
   df_akpegg3_miroc126, df_akpegg3_miroc585, akpegg1_list, akpegg2_list, akpegg3_list)

### Larvae ----
# Load data
df_akplarvae1_cesm126 <- readRDS(here('data', 'df_akplarvae1_cesm126.rds'))
df_akplarvae2_cesm126 <- readRDS(here('data', 'df_akplarvae2_cesm126.rds'))
df_akplarvae3_cesm126 <- readRDS(here('data', 'df_akplarvae3_cesm126.rds'))
df_akplarvae1_cesm585 <- readRDS(here('data', 'df_akplarvae1_cesm585.rds'))
df_akplarvae2_cesm585 <- readRDS(here('data', 'df_akplarvae2_cesm585.rds'))
df_akplarvae3_cesm585 <- readRDS(here('data', 'df_akplarvae3_cesm585.rds'))
df_akplarvae1_gfdl126 <- readRDS(here('data', 'df_akplarvae1_gfdl126.rds'))
df_akplarvae2_gfdl126 <- readRDS(here('data', 'df_akplarvae2_gfdl126.rds'))
df_akplarvae3_gfdl126 <- readRDS(here('data', 'df_akplarvae3_gfdl126.rds'))
df_akplarvae1_gfdl585 <- readRDS(here('data', 'df_akplarvae1_gfdl585.rds'))
df_akplarvae2_gfdl585 <- readRDS(here('data', 'df_akplarvae2_gfdl585.rds'))
df_akplarvae3_gfdl585 <- readRDS(here('data', 'df_akplarvae3_gfdl585.rds'))
df_akplarvae1_miroc126 <- readRDS(here('data', 'df_akplarvae1_miroc126.rds'))
df_akplarvae2_miroc126 <- readRDS(here('data', 'df_akplarvae2_miroc126.rds'))
df_akplarvae3_miroc126 <- readRDS(here('data', 'df_akplarvae3_miroc126.rds'))
df_akplarvae1_miroc585 <- readRDS(here('data', 'df_akplarvae1_miroc585.rds'))
df_akplarvae2_miroc585 <- readRDS(here('data', 'df_akplarvae2_miroc585.rds'))
df_akplarvae3_miroc585 <- readRDS(here('data', 'df_akplarvae3_miroc585.rds'))

akplarvae1_list <- list(df_akplarvae1_cesm126, df_akplarvae1_cesm585,
                        df_akplarvae1_gfdl126, df_akplarvae1_gfdl585,
                        df_akplarvae1_miroc126, df_akplarvae1_miroc585)

akplarvae2_list <- list(df_akplarvae2_cesm126, df_akplarvae2_cesm585,
                        df_akplarvae2_gfdl126, df_akplarvae2_gfdl585,
                        df_akplarvae2_miroc126, df_akplarvae2_miroc585)

akplarvae3_list <- list(df_akplarvae3_cesm126, df_akplarvae3_cesm585,
                        df_akplarvae3_gfdl126, df_akplarvae3_gfdl585,
                        df_akplarvae3_miroc126, df_akplarvae3_miroc585)

akplarvae_COG <- COG_time(akplarvae1_list, akplarvae2_list, akplarvae3_list)
saveRDS(akplarvae_COG, here('data', 'akplarvae_COG'))

rm(df_akplarvae1_cesm126, df_akplarvae1_cesm585, df_akplarvae1_gfdl126, df_akplarvae1_gfdl585,
   df_akplarvae1_miroc126, df_akplarvae1_miroc585, df_akplarvae2_cesm126, df_akplarvae2_cesm585, 
   df_akplarvae2_gfdl126, df_akplarvae2_gfdl585, df_akplarvae2_miroc126, df_akplarvae2_miroc585,
   df_akplarvae3_cesm126, df_akplarvae3_cesm585, df_akplarvae3_gfdl126, df_akplarvae3_gfdl585,
   df_akplarvae3_miroc126, df_akplarvae3_miroc585, akplarvae1_list, akplarvae2_list, akplarvae3_list)

## Yellowfin Sole ----
### Larvae ----
# Load data
df_yfslarvae1_cesm126 <- readRDS(here('data', 'df_yfslarvae1_cesm126.rds'))
df_yfslarvae2_cesm126 <- readRDS(here('data', 'df_yfslarvae2_cesm126.rds'))
df_yfslarvae3_cesm126 <- readRDS(here('data', 'df_yfslarvae3_cesm126.rds'))
df_yfslarvae1_cesm585 <- readRDS(here('data', 'df_yfslarvae1_cesm585.rds'))
df_yfslarvae2_cesm585 <- readRDS(here('data', 'df_yfslarvae2_cesm585.rds'))
df_yfslarvae3_cesm585 <- readRDS(here('data', 'df_yfslarvae3_cesm585.rds'))
df_yfslarvae1_gfdl126 <- readRDS(here('data', 'df_yfslarvae1_gfdl126.rds'))
df_yfslarvae2_gfdl126 <- readRDS(here('data', 'df_yfslarvae2_gfdl126.rds'))
df_yfslarvae3_gfdl126 <- readRDS(here('data', 'df_yfslarvae3_gfdl126.rds'))
df_yfslarvae1_gfdl585 <- readRDS(here('data', 'df_yfslarvae1_gfdl585.rds'))
df_yfslarvae2_gfdl585 <- readRDS(here('data', 'df_yfslarvae2_gfdl585.rds'))
df_yfslarvae3_gfdl585 <- readRDS(here('data', 'df_yfslarvae3_gfdl585.rds'))
df_yfslarvae1_miroc126 <- readRDS(here('data', 'df_yfslarvae1_miroc126.rds'))
df_yfslarvae2_miroc126 <- readRDS(here('data', 'df_yfslarvae2_miroc126.rds'))
df_yfslarvae3_miroc126 <- readRDS(here('data', 'df_yfslarvae3_miroc126.rds'))
df_yfslarvae1_miroc585 <- readRDS(here('data', 'df_yfslarvae1_miroc585.rds'))
df_yfslarvae2_miroc585 <- readRDS(here('data', 'df_yfslarvae2_miroc585.rds'))
df_yfslarvae3_miroc585 <- readRDS(here('data', 'df_yfslarvae3_miroc585.rds'))

yfslarvae1_list <- list(df_yfslarvae1_cesm126, df_yfslarvae1_cesm585,
                        df_yfslarvae1_gfdl126, df_yfslarvae1_gfdl585,
                        df_yfslarvae1_miroc126, df_yfslarvae1_miroc585)

yfslarvae2_list <- list(df_yfslarvae2_cesm126, df_yfslarvae2_cesm585,
                        df_yfslarvae2_gfdl126, df_yfslarvae2_gfdl585,
                        df_yfslarvae2_miroc126, df_yfslarvae2_miroc585)

yfslarvae3_list <- list(df_yfslarvae3_cesm126, df_yfslarvae3_cesm585,
                        df_yfslarvae3_gfdl126, df_yfslarvae3_gfdl585,
                        df_yfslarvae3_miroc126, df_yfslarvae3_miroc585)

yfslarvae_COG <- COG_time(yfslarvae1_list, yfslarvae2_list, yfslarvae3_list)
saveRDS(yfslarvae_COG, here('data', 'yfslarvae_COG'))

rm(df_yfslarvae1_cesm126, df_yfslarvae1_cesm585, df_yfslarvae1_gfdl126, df_yfslarvae1_gfdl585,
   df_yfslarvae1_miroc126, df_yfslarvae1_miroc585, df_yfslarvae2_cesm126, df_yfslarvae2_cesm585, 
   df_yfslarvae2_gfdl126, df_yfslarvae2_gfdl585, df_yfslarvae2_miroc126, df_yfslarvae2_miroc585,
   df_yfslarvae3_cesm126, df_yfslarvae3_cesm585, df_yfslarvae3_gfdl126, df_yfslarvae3_gfdl585,
   df_yfslarvae3_miroc126, df_yfslarvae3_miroc585, yfslarvae1_list, yfslarvae2_list, yfslarvae3_list)

## Northern Rock Sole ----
### Larvae ----
# Load data
df_nrslarvae1_cesm126 <- readRDS(here('data', 'df_nrslarvae1_cesm126.rds'))
df_nrslarvae2_cesm126 <- readRDS(here('data', 'df_nrslarvae2_cesm126.rds'))
df_nrslarvae3_cesm126 <- readRDS(here('data', 'df_nrslarvae3_cesm126.rds'))
df_nrslarvae1_cesm585 <- readRDS(here('data', 'df_nrslarvae1_cesm585.rds'))
df_nrslarvae2_cesm585 <- readRDS(here('data', 'df_nrslarvae2_cesm585.rds'))
df_nrslarvae3_cesm585 <- readRDS(here('data', 'df_nrslarvae3_cesm585.rds'))
df_nrslarvae1_gfdl126 <- readRDS(here('data', 'df_nrslarvae1_gfdl126.rds'))
df_nrslarvae2_gfdl126 <- readRDS(here('data', 'df_nrslarvae2_gfdl126.rds'))
df_nrslarvae3_gfdl126 <- readRDS(here('data', 'df_nrslarvae3_gfdl126.rds'))
df_nrslarvae1_gfdl585 <- readRDS(here('data', 'df_nrslarvae1_gfdl585.rds'))
df_nrslarvae2_gfdl585 <- readRDS(here('data', 'df_nrslarvae2_gfdl585.rds'))
df_nrslarvae3_gfdl585 <- readRDS(here('data', 'df_nrslarvae3_gfdl585.rds'))
df_nrslarvae1_miroc126 <- readRDS(here('data', 'df_nrslarvae1_miroc126.rds'))
df_nrslarvae2_miroc126 <- readRDS(here('data', 'df_nrslarvae2_miroc126.rds'))
df_nrslarvae3_miroc126 <- readRDS(here('data', 'df_nrslarvae3_miroc126.rds'))
df_nrslarvae1_miroc585 <- readRDS(here('data', 'df_nrslarvae1_miroc585.rds'))
df_nrslarvae2_miroc585 <- readRDS(here('data', 'df_nrslarvae2_miroc585.rds'))
df_nrslarvae3_miroc585 <- readRDS(here('data', 'df_nrslarvae3_miroc585.rds'))

nrslarvae1_list <- list(df_nrslarvae1_cesm126, df_nrslarvae1_cesm585,
                        df_nrslarvae1_gfdl126, df_nrslarvae1_gfdl585,
                        df_nrslarvae1_miroc126, df_nrslarvae1_miroc585)

nrslarvae2_list <- list(df_nrslarvae2_cesm126, df_nrslarvae2_cesm585,
                        df_nrslarvae2_gfdl126, df_nrslarvae2_gfdl585,
                        df_nrslarvae2_miroc126, df_nrslarvae2_miroc585)

nrslarvae3_list <- list(df_nrslarvae3_cesm126, df_nrslarvae3_cesm585,
                        df_nrslarvae3_gfdl126, df_nrslarvae3_gfdl585,
                        df_nrslarvae3_miroc126, df_nrslarvae3_miroc585)

nrslarvae_COG <- COG_time(nrslarvae1_list, nrslarvae2_list, nrslarvae3_list)
saveRDS(nrslarvae_COG, here('data', 'nrslarvae_COG'))

rm(df_nrslarvae1_cesm126, df_nrslarvae1_cesm585, df_nrslarvae1_gfdl126, df_nrslarvae1_gfdl585,
   df_nrslarvae1_miroc126, df_nrslarvae1_miroc585, df_nrslarvae2_cesm126, df_nrslarvae2_cesm585, 
   df_nrslarvae2_gfdl126, df_nrslarvae2_gfdl585, df_nrslarvae2_miroc126, df_nrslarvae2_miroc585,
   df_nrslarvae3_cesm126, df_nrslarvae3_cesm585, df_nrslarvae3_gfdl126, df_nrslarvae3_gfdl585,
   df_nrslarvae3_miroc126, df_nrslarvae3_miroc585, nrslarvae1_list, nrslarvae2_list, nrslarvae3_list)

## Pacific Cod ----
### Larvae ----
# Load data
df_pcodlarvae1_cesm126 <- readRDS(here('data', 'df_pcodlarvae1_cesm126.rds'))
df_pcodlarvae2_cesm126 <- readRDS(here('data', 'df_pcodlarvae2_cesm126.rds'))
df_pcodlarvae3_cesm126 <- readRDS(here('data', 'df_pcodlarvae3_cesm126.rds'))
df_pcodlarvae1_cesm585 <- readRDS(here('data', 'df_pcodlarvae1_cesm585.rds'))
df_pcodlarvae2_cesm585 <- readRDS(here('data', 'df_pcodlarvae2_cesm585.rds'))
df_pcodlarvae3_cesm585 <- readRDS(here('data', 'df_pcodlarvae3_cesm585.rds'))
df_pcodlarvae1_gfdl126 <- readRDS(here('data', 'df_pcodlarvae1_gfdl126.rds'))
df_pcodlarvae2_gfdl126 <- readRDS(here('data', 'df_pcodlarvae2_gfdl126.rds'))
df_pcodlarvae3_gfdl126 <- readRDS(here('data', 'df_pcodlarvae3_gfdl126.rds'))
df_pcodlarvae1_gfdl585 <- readRDS(here('data', 'df_pcodlarvae1_gfdl585.rds'))
df_pcodlarvae2_gfdl585 <- readRDS(here('data', 'df_pcodlarvae2_gfdl585.rds'))
df_pcodlarvae3_gfdl585 <- readRDS(here('data', 'df_pcodlarvae3_gfdl585.rds'))
df_pcodlarvae1_miroc126 <- readRDS(here('data', 'df_pcodlarvae1_miroc126.rds'))
df_pcodlarvae2_miroc126 <- readRDS(here('data', 'df_pcodlarvae2_miroc126.rds'))
df_pcodlarvae3_miroc126 <- readRDS(here('data', 'df_pcodlarvae3_miroc126.rds'))
df_pcodlarvae1_miroc585 <- readRDS(here('data', 'df_pcodlarvae1_miroc585.rds'))
df_pcodlarvae2_miroc585 <- readRDS(here('data', 'df_pcodlarvae2_miroc585.rds'))
df_pcodlarvae3_miroc585 <- readRDS(here('data', 'df_pcodlarvae3_miroc585.rds'))

pcodlarvae1_list <- list(df_pcodlarvae1_cesm126, df_pcodlarvae1_cesm585,
                        df_pcodlarvae1_gfdl126, df_pcodlarvae1_gfdl585,
                        df_pcodlarvae1_miroc126, df_pcodlarvae1_miroc585)

pcodlarvae2_list <- list(df_pcodlarvae2_cesm126, df_pcodlarvae2_cesm585,
                        df_pcodlarvae2_gfdl126, df_pcodlarvae2_gfdl585,
                        df_pcodlarvae2_miroc126, df_pcodlarvae2_miroc585)

pcodlarvae3_list <- list(df_pcodlarvae3_cesm126, df_pcodlarvae3_cesm585,
                        df_pcodlarvae3_gfdl126, df_pcodlarvae3_gfdl585,
                        df_pcodlarvae3_miroc126, df_pcodlarvae3_miroc585)

pcodlarvae_COG <- COG_time(pcodlarvae1_list, pcodlarvae2_list, pcodlarvae3_list)
saveRDS(pcodlarvae_COG, here('data', 'pcodlarvae_COG'))

rm(df_pcodlarvae1_cesm126, df_pcodlarvae1_cesm585, df_pcodlarvae1_gfdl126, df_pcodlarvae1_gfdl585,
   df_pcodlarvae1_miroc126, df_pcodlarvae1_miroc585, df_pcodlarvae2_cesm126, df_pcodlarvae2_cesm585, 
   df_pcodlarvae2_gfdl126, df_pcodlarvae2_gfdl585, df_pcodlarvae2_miroc126, df_pcodlarvae2_miroc585,
   df_pcodlarvae3_cesm126, df_pcodlarvae3_cesm585, df_pcodlarvae3_gfdl126, df_pcodlarvae3_gfdl585,
   df_pcodlarvae3_miroc126, df_pcodlarvae3_miroc585, pcodlarvae1_list, pcodlarvae2_list, pcodlarvae3_list)