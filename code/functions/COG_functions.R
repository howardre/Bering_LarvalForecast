# Calculate COG for predictions
COG_pred <- function(data){
  temp1 <- data %>% dplyr::mutate(pred = pred_scaled * lat) %>%
    dplyr::summarize(value = sum(pred, na.rm = TRUE)) %>%
    dplyr::select(value)
  
  temp2 <- data %>%
    dplyr::summarize(value = sum(pred_scaled, na.rm = TRUE)) %>%
    dplyr::select(value)
  
  COG_x <- temp1 / temp2
  
  temp3 <- data %>% dplyr::mutate(pred = pred_scaled * lon) %>%
    dplyr::summarize(value1 = sum(pred, na.rm = TRUE)) %>%
    dplyr::select(value1)
  
  temp4 <- data %>%
    dplyr::summarize(value1 = sum(pred_scaled, na.rm = TRUE)) %>%
    dplyr::select(value1)
  
  COG_y <- temp3 / temp4
  
  COG <- as.data.frame(cbind(COG_x, COG_y))
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
  
  COG_x <- temp1 / temp2
  
  temp3 <- data %>% dplyr::mutate(temp = temperature * lon) %>%
    dplyr::summarize(value1 = sum(temp, na.rm = TRUE)) %>%
    dplyr::select(value1)
  
  temp4 <- data %>%
    dplyr::summarize(value1 = sum(temperature, na.rm = TRUE)) %>%
    dplyr::select(value1)
  
  COG_y <- temp3 / temp4
  
  COG <- as.data.frame(cbind(COG_x, COG_y))
  return(COG)
}

COG_calc <- function(hindcast, my_list1, my_list2, my_list3, my_list4, my_list5, my_list6){
  keys1 <- unique(names(my_list1[[1]])) # create list of names of year data frames
  keys2 <- unique(names(my_list2[[1]]))
  keys3 <- unique(names(my_list3[[1]]))
  abun_years1 <- lapply(keys1, function(x) colMeans(do.call(rbind, lapply(my_list1, "[[", x)), na.rm = TRUE)) # average per year
  abun_years2 <- lapply(keys2, function(x) colMeans(do.call(rbind, lapply(my_list2, "[[", x)), na.rm = TRUE))
  abun_years3 <- lapply(keys3, function(x) colMeans(do.call(rbind, lapply(my_list3, "[[", x)), na.rm = TRUE))
  temp_years1 <- lapply(keys1, function(x) colMeans(do.call(rbind, lapply(my_list4, "[[", x)), na.rm = TRUE))
  temp_years2 <- lapply(keys2, function(x) colMeans(do.call(rbind, lapply(my_list5, "[[", x)), na.rm = TRUE))
  temp_years3 <- lapply(keys3, function(x) colMeans(do.call(rbind, lapply(my_list6, "[[", x)), na.rm = TRUE))
  hindcast_COG <- data.frame(avg_lat = rowMeans(do.call(cbind, lapply(hindcast[[2]], "[", "value")), na.rm = TRUE),
                             avg_lon = rowMeans(do.call(cbind, lapply(hindcast[[2]], "[", "value1")), na.rm = TRUE))
  historic_COG <- data.frame(avg_lat = rowMeans(do.call(cbind, lapply(hindcast[[3]], "[", "value")), na.rm = TRUE),
                             avg_lon = rowMeans(do.call(cbind, lapply(hindcast[[3]], "[", "value1")), na.rm = TRUE))
  overall_COG1 <- data.frame(avg_lat = rowMeans(do.call(cbind, lapply(abun_years1, "[", "value")), na.rm = TRUE),
                             avg_lon = rowMeans(do.call(cbind, lapply(abun_years1, "[", "value1")), na.rm = TRUE)) # average per period
  overall_COG2 <- data.frame(avg_lat = rowMeans(do.call(cbind, lapply(abun_years2, "[", "value")), na.rm = TRUE),
                             avg_lon = rowMeans(do.call(cbind, lapply(abun_years2, "[", "value1")), na.rm = TRUE))
  overall_COG3 <- data.frame(avg_lat = rowMeans(do.call(cbind, lapply(abun_years3, "[", "value")), na.rm = TRUE),
                             avg_lon = rowMeans(do.call(cbind, lapply(abun_years3, "[", "value1")), na.rm = TRUE))
  climate_COG1 <- data.frame(avg_lat = rowMeans(do.call(cbind, lapply(temp_years1, "[", "value")), na.rm = TRUE),
                             avg_lon = rowMeans(do.call(cbind, lapply(temp_years1, "[", "value1")), na.rm = TRUE))
  climate_COG2 <- data.frame(avg_lat = rowMeans(do.call(cbind, lapply(temp_years2, "[", "value")), na.rm = TRUE),
                             avg_lon = rowMeans(do.call(cbind, lapply(temp_years2, "[", "value1")), na.rm = TRUE))
  climate_COG3 <- data.frame(avg_lat = rowMeans(do.call(cbind, lapply(temp_years3, "[", "value")), na.rm = TRUE),
                             avg_lon = rowMeans(do.call(cbind, lapply(temp_years3, "[", "value1")), na.rm = TRUE))
  lat_abun <- rbind(hindcast_COG$avg_lat, overall_COG1$avg_lat, overall_COG2$avg_lat, overall_COG3$avg_lat) # combine into data frame
  lat_clim <- rbind(historic_COG$avg_lat, climate_COG1$avg_lat, climate_COG2$avg_lat, climate_COG3$avg_lat)
  lon_abun <- rbind(hindcast_COG$avg_lon, overall_COG1$avg_lon, overall_COG2$avg_lon, overall_COG3$avg_lon)
  lon_clim <- rbind(historic_COG$avg_lon, climate_COG1$avg_lon, climate_COG2$avg_lon, climate_COG3$avg_lon)
  avg_abun <- cbind(lat_abun, lon_abun)
  avg_clim <- cbind(lat_clim, lon_clim)
  colnames(avg_abun) <- c("lat", "lon")
  rownames(avg_abun) <- c("hindcast", "time1", "time2", "time3")
  colnames(avg_clim) <- c("lat", "lon")
  rownames(avg_clim) <- c("hindcast", "time1", "time2", "time3")
  avg_abun <- as.data.frame(avg_abun)
  avg_clim <- as.data.frame(avg_clim)
  avg_data <- as.data.frame(tibble::rownames_to_column(avg_abun, "period"))
  avg_temp <- as.data.frame(tibble::rownames_to_column(avg_clim, "period"))
  avg_dist <- mutate(avg_data,
                     distance = distHaversine(cbind(lon, lat),
                                              cbind(lag(lon), 
                                                    lag(lat))) / 1000)
  start_dist <- distHaversine(c(avg_data$lon[1], avg_data$lat[1]),
                              c(avg_data$lon[4], avg_data$lat[4])) / 1000
  COG_list <- list(hindcast,
                   abun_years1, abun_years2, abun_years3, 
                   temp_years1, temp_years2, temp_years3, 
                   avg_temp, avg_dist, start_dist)
  return(COG_list)
}

# Make loop to calculate per time frame
# COG_loop1 <- function(data, the_list){
#   for (i in 1:length(data)){
#     the_list[[i]] <- COG_pred(data[[i]])
#   }
#   print(the_list)
# }
# 
# COG_loop2 <- function(data, the_list){
#   for (i in 1:length(data)){
#     the_list[[i]] <- COG_temp(data[[i]])
#   }
#   print(the_list)
# }


# Get the COG for each time period
# COG_time <- function(hindcast, list1, list2, list3){
#   time1_COG <- list() # create empty lists
#   time2_COG <- list()
#   time3_COG <- list()
#   temp1_COG <- list()
#   temp2_COG <- list()
#   temp3_COG <- list()
#   
#   hindcast_COG <- COG_pred(hindcast) # calculate hindcast prediction
#   time1_COG <- COG_loop1(list1, time1_COG) # calculate predictions for each within period
#   time2_COG <- COG_loop1(list2, time2_COG)
#   time3_COG <- COG_loop1(list3, time3_COG)
#   hindcast_temp <- COG_temp(hindcast)
#   temp1_COG <- COG_loop2(list1, temp1_COG) # calculate temperature for each within period
#   temp2_COG <- COG_loop2(list2, temp2_COG)
#   temp3_COG <- COG_loop2(list3, temp3_COG)
#   
#   avg1_lat <- rowMeans(do.call(cbind, lapply(time1_COG, "[", "value")), na.rm = TRUE)
#   avg1_lon <- rowMeans(do.call(cbind, lapply(time1_COG, "[", "value1")), na.rm = TRUE)
#   
#   avg2_lat <- rowMeans(do.call(cbind, lapply(time2_COG, "[", "value")), na.rm = TRUE)
#   avg2_lon <- rowMeans(do.call(cbind, lapply(time2_COG, "[", "value1")), na.rm = TRUE)
#   
#   avg3_lat <- rowMeans(do.call(cbind, lapply(time3_COG, "[", "value")), na.rm = TRUE)
#   avg3_lon <- rowMeans(do.call(cbind, lapply(time3_COG, "[", "value1")), na.rm = TRUE)
#   
#   lats <- rbind(hindcast_COG$value, avg1_lat, avg2_lat, avg3_lat)
#   lons <- rbind(hindcast_COG$value1, avg1_lon, avg2_lon, avg3_lon)
#   
#   avg4_lat <- rowMeans(do.call(cbind, lapply(temp1_COG, "[", "value")))
#   avg4_lon <- rowMeans(do.call(cbind, lapply(temp1_COG, "[", "value1")))
#   
#   avg5_lat <- rowMeans(do.call(cbind, lapply(temp2_COG, "[", "value")))
#   avg5_lon <- rowMeans(do.call(cbind, lapply(temp2_COG, "[", "value1")))
#   
#   avg6_lat <- rowMeans(do.call(cbind, lapply(temp3_COG, "[", "value")))
#   avg6_lon <- rowMeans(do.call(cbind, lapply(temp3_COG, "[", "value1")))
#   
#   lats1 <- rbind(hindcast_COG$value, avg4_lat, avg5_lat, avg6_lat)
#   lons1 <- rbind(hindcast_COG$value1, avg4_lon, avg5_lon, avg6_lon)
#   
#   avg_df <- cbind(lats, lons)
#   avg_df1 <- cbind(lats1, lons1)
#   colnames(avg_df) <- c("lat", "lon")
#   rownames(avg_df) <- c("hindcast", "time1", "time2", "time3")
#   colnames(avg_df1) <- c("lat", "lon")
#   rownames(avg_df1) <- c("hindcast", "time1", "time2", "time3")
#   avg_df <- as.data.frame(avg_df)
#   avg_df1 <- as.data.frame(avg_df1)
#   avg_data <- as.data.frame(tibble::rownames_to_column(avg_df, "period"))
#   avg_temp <- as.data.frame(tibble::rownames_to_column(avg_df1, "period"))
#   avg_dist <- mutate(avg_data,
#                      distance = distHaversine(cbind(lon, lat),
#                                               cbind(lag(lon), lag(lat))) / 1000)
#   start_dist <- distHaversine(c(avg_data$lon[1], avg_data$lat[1]),
#                               c(avg_data$lon[4], avg_data$lat[4])) / 1000
#   final_df <- list(avg_temp, avg_dist, start_dist)
#   return(final_df)
# }
# 
# # Scatterplot of COGs
plot_COG <- function(COG, species){
  ggplot() +
    geom_path(data = COG[[9]],
              aes(x = lon, y = lat, linetype = "Species"),
              size = 1) +
    geom_point(data = COG[[9]],
               aes(x = lon, y = lat,
                   color = period), 
               size = 5) +
    geom_path(data = COG[[8]],
              aes(x = lon, y = lat, linetype = "Temperature"),
              size = 1) +
    geom_point(data = COG[[8]],
               aes(x = lon, y = lat,
                   color = period), 
               size = 5) +
    scale_linetype_manual(values = c("solid", "dotted")) +
    scale_color_manual(values = c("goldenrod3", "coral2", "darkturquoise", "darkslateblue"),
                       labels = c("Hindcast", "2015-2039", "2040-2069", "2070-2099")) +
    theme_classic() +
    theme(axis.ticks = element_blank(),
          axis.text = element_text(family = "serif", size = 18),
          plot.title = element_text(size = 24, family = "serif", 
                                    face = "bold", hjust = 0.5),
          axis.title = element_text(family = "serif", size = 22),
          axis.text.x = element_text(angle = 45, vjust = 0.7),
          legend.title = element_text(family = "serif", size = 22,
                                      face = "bold"),
          legend.text = element_text(family = "serif", size = 20)) +
    labs(x = "Longitude \u00B0W",
         y = "Latitude \u00B0N",
         color = "Time Period",
         linetype = "COG",
         title = species)
}
# 
# Calculate distance between life stage COG for each time period
lifestage_dist <- function(data1, data2){
  d1 <- distHaversine(c(data1[[9]]$lon[1], data1[[9]]$lat[1]),
                      c(data2[[9]]$lon[1], data2[[9]]$lat[1])) / 1000
  d2 <- distHaversine(c(data1[[9]]$lon[2], data1[[9]]$lat[2]),
                      c(data2[[9]]$lon[2], data2[[9]]$lat[2])) / 1000
  d3 <- distHaversine(c(data1[[9]]$lon[3], data1[[9]]$lat[3]),
                      c(data2[[9]]$lon[3], data2[[9]]$lat[3])) / 1000
  d4 <- distHaversine(c(data1[[9]]$lon[4], data1[[9]]$lat[4]),
                      c(data2[[9]]$lon[4], data2[[9]]$lat[4])) / 1000
  distances <- list(d1, d2, d3, d4)
  return(distances)
}

temp_dist <- function(data){
  d1 <- distHaversine(c(data[[8]]$lon[1], data[[8]]$lat[1]),
                      c(data[[9]]$lon[1], data[[9]]$lat[1])) / 1000
  d2 <- distHaversine(c(data[[8]]$lon[2], data[[8]]$lat[2]),
                      c(data[[9]]$lon[2], data[[9]]$lat[2])) / 1000
  d3 <- distHaversine(c(data[[8]]$lon[3], data[[8]]$lat[3]),
                      c(data[[9]]$lon[3], data[[9]]$lat[3])) / 1000
  d4 <- distHaversine(c(data[[8]]$lon[4], data[[8]]$lat[4]),
                      c(data[[9]]$lon[4], data[[9]]$lat[4])) / 1000
  distances <- list(d1, d2, d3, d4)
  return(distances)
}

# period_dist <- function(data){
#   p1 <- (data[[1]][[2]][[length(data[[1]][[2]])]]$value - # hindcast has different # of years depending on species
#            data[[1]][[2]][[1]]$value) / length(data[[1]][[2]])
#   p2 <- as.vector((data[[2]][[25]][1] - data[[2]][[1]][1]) / 25)
#   p3 <- as.vector((data[[3]][[30]][1] - data[[3]][[1]][1]) / 30)
#   p4 <- as.vector((data[[4]][[30]][1] - data[[4]][[1]][1]) / 30)
#   t1 <- (data[[1]][[3]][[length(data[[1]][[3]])]]$value - 
#            data[[1]][[3]][[1]]$value) / length(data[[1]][[3]])
#   t2 <- as.vector((data[[5]][[25]][1] - data[[5]][[1]][1]) / 25)
#   t3 <- as.vector((data[[6]][[30]][1] - data[[6]][[1]][1]) / 30)
#   t4 <- as.vector((data[[7]][[30]][1] - data[[7]][[1]][1]) / 30)
#   abundance <- c(p1, p2, p3, p4)
#   temperature <- c(t1, t2, t3, t4)
#   final_df <- data.frame(abundance, temperature)
#   return(final_df)
# }

velocity_calc <- function(data_COG){
  data_hindcast <- as.data.frame(cbind(lapply(data_COG[[1]][[2]], "[[", 1), deparse.level = 1))
  data_hindcast$latitude <- as.numeric(data_hindcast$V1)
  data_hindcast <- tibble::rowid_to_column(data_hindcast, "id")
  data_period1 <- as.data.frame(cbind(lapply(data_COG[[2]], "[[", 1), deparse.level = 1))
  data_period1$latitude <- as.numeric(data_period1$V1)
  data_period1 <- tibble::rowid_to_column(data_period1, "id")
  data_period2 <- as.data.frame(cbind(lapply(data_COG[[3]], "[[", 1), deparse.level = 1))
  data_period2$latitude <- as.numeric(data_period2$V1)
  data_period2 <- tibble::rowid_to_column(data_period2, "id")
  data_period3 <- as.data.frame(cbind(lapply(data_COG[[4]], "[[", 1), deparse.level = 1))
  data_period3$latitude <- as.numeric(data_period3$V1)
  data_period3 <- tibble::rowid_to_column(data_period3, "id")
  temp_hindcast <- as.data.frame(cbind(lapply(data_COG[[1]][[3]], "[[", 1), deparse.level = 1))
  temp_hindcast$latitude <- as.numeric(temp_hindcast$V1)
  temp_hindcast <- tibble::rowid_to_column(temp_hindcast, "id")
  temp_period1 <- as.data.frame(cbind(lapply(data_COG[[5]], "[[", 1), deparse.level = 1))
  temp_period1$latitude <- as.numeric(temp_period1$V1)
  temp_period1 <- tibble::rowid_to_column(temp_period1, "id")
  temp_period2 <- as.data.frame(cbind(lapply(data_COG[[6]], "[[", 1), deparse.level = 1))
  temp_period2$latitude <- as.numeric(temp_period2$V1)
  temp_period2 <- tibble::rowid_to_column(temp_period2, "id")
  temp_period3 <- as.data.frame(cbind(lapply(data_COG[[7]], "[[", 1), deparse.level = 1))
  temp_period3$latitude <- as.numeric(temp_period3$V1)
  temp_period3 <- tibble::rowid_to_column(temp_period3, "id")
  hindcast_lm1 <- lm(formula = latitude ~ id, data = data_hindcast, na.action = na.omit)
  hindcast_lm2 <- lm(formula = latitude ~ id, data = temp_hindcast, na.action = na.omit)
  data_lm1 <- lm(formula = latitude ~ id, data = data_period1, na.action = na.omit)
  data_lm2 <- lm(formula = latitude ~ id, data = data_period2, na.action = na.omit)
  data_lm3 <- lm(formula = latitude ~ id, data = data_period3, na.action = na.omit)
  temp_lm1 <- lm(formula = latitude ~ id, data = temp_period1, na.action = na.omit)
  temp_lm2 <- lm(formula = latitude ~ id, data = temp_period2, na.action = na.omit)
  temp_lm3 <- lm(formula = latitude ~ id, data = temp_period3, na.action = na.omit)
  hindcast_slope1 <- hindcast_lm1$coefficients[2]
  hindcast_slope2 <- hindcast_lm2$coefficients[2]
  data_slope1 <- data_lm1$coefficients[2]
  data_slope2 <- data_lm2$coefficients[2]
  data_slope3 <- data_lm3$coefficients[2]
  temp_slope1 <- temp_lm1$coefficients[2]
  temp_slope2 <- temp_lm2$coefficients[2]
  temp_slope3 <- temp_lm3$coefficients[2]
  return(list(data_slope1, data_slope2, data_slope3,
              temp_slope1, temp_slope2, temp_slope3,
              hindcast_slope1, hindcast_slope2))
}
