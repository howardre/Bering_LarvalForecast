# Calculate COG for predictions
COG_pred <- function(data){
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
  return(COG)
}

# Calculate COG for temperature
COG_temp <- function(data){
  temp1 <- data %>% dplyr::mutate(pred = temperature * lat) %>%
    dplyr::summarize(value = sum(pred, na.rm = TRUE)) %>%
    dplyr::select(value)
  
  temp2 <- data %>%
    dplyr::summarize(value = sum(temperature, na.rm = TRUE)) %>%
    dplyr::select(value)
  
  COG_y <- temp1 / temp2
  
  temp3 <- data %>% dplyr::mutate(pred = temperature * lon) %>%
    dplyr::summarize(value1 = sum(pred, na.rm = TRUE)) %>%
    dplyr::select(value1)
  
  temp4 <- data %>%
    dplyr::summarize(value1 = sum(temperature, na.rm = TRUE)) %>%
    dplyr::select(value1)
  
  COG_x <- temp3 / temp4
  
  COG <- as.data.frame(cbind(COG_y, COG_x))
  return(COG)
}

# Make loop to calculate per time frame
COG_loop <- function(data, the_list){
  for (i in 1:length(data)){
    the_list[[i]] <- COG(data[[i]])
  }
  print(the_list)
}

# Get the COG for each time period
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
  avg_df <- as.data.frame(avg_df)
  avg_data <- as.data.frame(tibble::rownames_to_column(avg_df, "period"))
  avg_dist <- mutate(avg_data,
                     distance = distHaversine(cbind(lon, lat),
                                              cbind(lag(lon), lag(lat))) / 1000)
  start_dist <- distHaversine(c(avg_data$lon[1], avg_data$lat[1]),
                              c(avg_data$lon[4], avg_data$lat[4])) / 1000
  final_df <- list(avg_dist, start_dist)
  return(final_df)
}

# Scatterplot of COGs
plot_COG <- function(COG){
  ggplot(data = COG[[1]], 
         aes(x = lon, y = lat)) +
    geom_point(aes(color = period), size = 4) +
    geom_path() +
    theme_classic() +
    theme(axis.ticks = element_blank(),
          axis.text = element_text(family = "serif", size = 16),
          axis.title = element_text(family = "serif", size = 20),
          axis.text.x = element_text(angle = 45, vjust = 0.7),
          strip.text = element_text(family = "serif", size = 20),
          legend.title = element_text(family = "serif", size = 18),
          legend.text = element_text(family = "serif", size = 16)) +
    labs(x = "Longitude",
         y = "Latitude")
}

# Calculate distance between life stage COG for each time period
lifestage_dist <- function(data1, data2){
  d1 <- distHaversine(c(data1[[1]]$lon[1], data1[[1]]$lat[1]),
                      c(data2[[1]]$lon[1], data2[[1]]$lat[1])) / 1000
  d2 <- distHaversine(c(data1[[1]]$lon[2], data1[[1]]$lat[2]),
                      c(data2[[1]]$lon[2], data2[[1]]$lat[2])) / 1000
  d3 <- distHaversine(c(data1[[1]]$lon[3], data1[[1]]$lat[3]),
                      c(data2[[1]]$lon[3], data2[[1]]$lat[3])) / 1000
  d4 <- distHaversine(c(data1[[1]]$lon[4], data1[[1]]$lat[4]),
                      c(data2[[1]]$lon[4], data2[[1]]$lat[4])) / 1000
  distances <- list(d1, d2, d3, d4)
  return(distances)
}