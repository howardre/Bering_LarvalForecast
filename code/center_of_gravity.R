# Title: Examine primary attributes of Bering 10K model.
# Authors: Lorenzo Ciannelli, Rebecca Howard
# Date: 07/17/2023

# Load data
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

time1_list <- list(df_pkegg1_cesm126, df_pkegg1_cesm585,
                   df_pkegg1_gfdl126, df_pkegg1_gfdl585,
                   df_pkegg1_miroc126, df_pkegg1_miroc585)

time2_list <- list(df_pkegg2_cesm126, df_pkegg2_cesm585,
                   df_pkegg2_gfdl126, df_pkegg2_gfdl585,
                   df_pkegg2_miroc126, df_pkegg2_miroc585)

time3_list <- list(df_pkegg3_cesm126, df_pkegg3_cesm585,
                   df_pkegg3_gfdl126, df_pkegg3_gfdl585,
                   df_pkegg3_miroc126, df_pkegg3_miroc585)

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
COG_time <- function(list1, list2, list3){
  time1_COG <- list()
  time2_COG <- list()
  time3_COG <- list()
  
  time1_COG <- COG_loop(list1, time1_COG)
  time2_COG <- COG_loop(list2, time2_COG)
  time3_COG <- COG_loop(list3, time3_COG)
  
  avg1_lat <- rowMeans(do.call(cbind, lapply(time1_COG, "[", "value")))
  avg1_lon <- rowMeans(do.call(cbind, lapply(time1_COG, "[", "value1")))
  
  avg2_lat <- rowMeans(do.call(cbind, lapply(time2_COG, "[", "value")))
  avg2_lon <- rowMeans(do.call(cbind, lapply(time2_COG, "[", "value1")))
  
  avg3_lat <- rowMeans(do.call(cbind, lapply(time3_COG, "[", "value")))
  avg3_lon <- rowMeans(do.call(cbind, lapply(time3_COG, "[", "value1")))
  
  avg_list <- list(avg1_lat, avg1_lon, avg2_lat, avg2_lon, avg3_lat, avg3_lon)
  return(avg_list)
}

pkegg_COG <- COG_time(time1_list, time2_list, time3_list)
