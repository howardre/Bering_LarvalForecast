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
  temp1 <- data %>% dplyr::mutate(temp = temperature * lat) %>%
    dplyr::summarize(value = sum(temp, na.rm = TRUE)) %>%
    dplyr::select(value)
  
  temp2 <- data %>%
    dplyr::summarize(value = sum(temperature, na.rm = TRUE)) %>%
    dplyr::select(value)
  
  COG_y <- temp1 / temp2
  
  temp3 <- data %>% dplyr::mutate(temp = temperature * lon) %>%
    dplyr::summarize(value1 = sum(temp, na.rm = TRUE)) %>%
    dplyr::select(value1)
  
  temp4 <- data %>%
    dplyr::summarize(value1 = sum(temperature, na.rm = TRUE)) %>%
    dplyr::select(value1)
  
  COG_x <- temp3 / temp4
  
  COG <- as.data.frame(cbind(COG_y, COG_x))
  return(COG)
}

# Make loop to calculate per time frame
COG_loop1 <- function(data, the_list){
  for (i in 1:length(data)){
    the_list[[i]] <- COG_pred(data[[i]])
  }
  print(the_list)
}

COG_loop2 <- function(data, the_list){
  for (i in 1:length(data)){
    the_list[[i]] <- COG_temp(data[[i]])
  }
  print(the_list)
}


# Get the COG for each time period
COG_time <- function(hindcast, list1, list2, list3){
  time1_COG <- list() # create empty lists
  time2_COG <- list()
  time3_COG <- list()
  temp1_COG <- list()
  temp2_COG <- list()
  temp3_COG <- list()
  
  hindcast_COG <- COG_pred(hindcast) # calculate hindcast prediction
  time1_COG <- COG_loop1(list1, time1_COG) # calculate predictions for each within period
  time2_COG <- COG_loop1(list2, time2_COG)
  time3_COG <- COG_loop1(list3, time3_COG)
  hindcast_temp <- COG_temp(hindcast)
  temp1_COG <- COG_loop2(list1, temp1_COG) # calculate temperature for each within period
  temp2_COG <- COG_loop2(list2, temp2_COG)
  temp3_COG <- COG_loop2(list3, temp3_COG)
  
  avg1_lat <- rowMeans(do.call(cbind, lapply(time1_COG, "[", "value")), na.rm = TRUE)
  avg1_lon <- rowMeans(do.call(cbind, lapply(time1_COG, "[", "value1")), na.rm = TRUE)
  
  avg2_lat <- rowMeans(do.call(cbind, lapply(time2_COG, "[", "value")), na.rm = TRUE)
  avg2_lon <- rowMeans(do.call(cbind, lapply(time2_COG, "[", "value1")), na.rm = TRUE)
  
  avg3_lat <- rowMeans(do.call(cbind, lapply(time3_COG, "[", "value")), na.rm = TRUE)
  avg3_lon <- rowMeans(do.call(cbind, lapply(time3_COG, "[", "value1")), na.rm = TRUE)
  
  lats <- rbind(hindcast_COG$value, avg1_lat, avg2_lat, avg3_lat)
  lons <- rbind(hindcast_COG$value1, avg1_lon, avg2_lon, avg3_lon)
  
  avg4_lat <- rowMeans(do.call(cbind, lapply(temp1_COG, "[", "value")))
  avg4_lon <- rowMeans(do.call(cbind, lapply(temp1_COG, "[", "value1")))
  
  avg5_lat <- rowMeans(do.call(cbind, lapply(temp2_COG, "[", "value")))
  avg5_lon <- rowMeans(do.call(cbind, lapply(temp2_COG, "[", "value1")))
  
  avg6_lat <- rowMeans(do.call(cbind, lapply(temp3_COG, "[", "value")))
  avg6_lon <- rowMeans(do.call(cbind, lapply(temp3_COG, "[", "value1")))
  
  lats1 <- rbind(hindcast_COG$value, avg4_lat, avg5_lat, avg6_lat)
  lons1 <- rbind(hindcast_COG$value1, avg4_lon, avg5_lon, avg6_lon)
  
  avg_df <- cbind(lats, lons)
  avg_df1 <- cbind(lats1, lons1)
  colnames(avg_df) <- c("lat", "lon")
  rownames(avg_df) <- c("hindcast", "time1", "time2", "time3")
  colnames(avg_df1) <- c("lat", "lon")
  rownames(avg_df1) <- c("hindcast", "time1", "time2", "time3")
  avg_df <- as.data.frame(avg_df)
  avg_df1 <- as.data.frame(avg_df1)
  avg_data <- as.data.frame(tibble::rownames_to_column(avg_df, "period"))
  avg_temp <- as.data.frame(tibble::rownames_to_column(avg_df1, "period"))
  avg_dist <- mutate(avg_data,
                     distance = distHaversine(cbind(lon, lat),
                                              cbind(lag(lon), lag(lat))) / 1000)
  start_dist <- distHaversine(c(avg_data$lon[1], avg_data$lat[1]),
                              c(avg_data$lon[4], avg_data$lat[4])) / 1000
  final_df <- list(avg_temp, avg_dist, start_dist)
  return(final_df)
}

# Scatterplot of COGs
plot_COG <- function(COG){
  ggplot(data = COG[[2]], 
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
  d1 <- distHaversine(c(data1[[2]]$lon[1], data1[[2]]$lat[1]),
                      c(data2[[2]]$lon[1], data2[[2]]$lat[1])) / 1000
  d2 <- distHaversine(c(data1[[2]]$lon[2], data1[[2]]$lat[2]),
                      c(data2[[2]]$lon[2], data2[[2]]$lat[2])) / 1000
  d3 <- distHaversine(c(data1[[2]]$lon[3], data1[[2]]$lat[3]),
                      c(data2[[2]]$lon[3], data2[[2]]$lat[3])) / 1000
  d4 <- distHaversine(c(data1[[2]]$lon[4], data1[[2]]$lat[4]),
                      c(data2[[2]]$lon[4], data2[[2]]$lat[4])) / 1000
  distances <- list(d1, d2, d3, d4)
  return(distances)
}

temp_dist <- function(data){
  d1 <- distHaversine(c(data[[1]]$lon[1], data[[1]]$lat[1]),
                      c(data[[2]]$lon[1], data[[2]]$lat[1])) / 1000
  d2 <- distHaversine(c(data[[1]]$lon[2], data[[1]]$lat[2]),
                      c(data[[2]]$lon[2], data[[2]]$lat[2])) / 1000
  d3 <- distHaversine(c(data[[1]]$lon[3], data[[1]]$lat[3]),
                      c(data[[2]]$lon[3], data[[2]]$lat[3])) / 1000
  d4 <- distHaversine(c(data[[1]]$lon[4], data[[1]]$lat[4]),
                      c(data[[2]]$lon[4], data[[2]]$lat[4])) / 1000
  distances <- list(d1, d2, d3, d4)
  return(distances)
}
