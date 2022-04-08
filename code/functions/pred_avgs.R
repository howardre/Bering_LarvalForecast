# Prediction function for the first time period (2015 - 2039)
predict_cells <- function(range, data, day,
                          month, scenario, temps,
                          salts, formula){
  preds <- pred_loop(range, data, day,
                     month, scenario, temps,
                     salts, formula)
  df <- data.frame(lat = preds[[1]]$lat,
                   lon = preds[[1]]$lon,
                   avg_pred = rowMedians(as.matrix(do.call(cbind, lapply(preds > quantile(preds, 0.95, na.rm = T), "[", "pred")))))
  df$pred_scaled <- rescale(df$avg_pred)
  return(df)
}


d[[x]][d[[x]]$pred > quantile(d[[x]]$pred, 0.95, na.rm = T), ]

# Prediction function for the second and third time periods (2040 - 2069, 2070 - 2099)
pred_list2 <- function(range, data, day, 
                       month, scenario, temps, 
                       salts, formula){
  preds <- pred_loop(range, data, day,
                     month, scenario, temps,
                     salts, formula)
  # Combine into one data frame
  df <- list(preds[[1]], preds[[2]],
             preds[[3]], preds[[4]],
             preds[[5]], preds[[6]],
             preds[[7]], preds[[8]],
             preds[[9]], preds[[10]],
             preds[[11]], preds[[12]],
             preds[[13]], preds[[14]],
             preds[[15]], preds[[16]],
             preds[[17]], preds[[18]],
             preds[[19]], preds[[20]],
             preds[[21]], preds[[22]],
             preds[[23]], preds[[24]],
             preds[[25]], preds[[26]],
             preds[[27]], preds[[28]],
             preds[[29]], preds[[30]]) %>%
    reduce(inner_join, by = c("lon", "lat", "dist", "doy")) 
  return(df)
}

# Function to generate average prediction from all predictions for first time period
pred_avgs1 <- function(df){
x <- grepl("pred", names(df), fixed = T)
avgs <- data.frame(lat = df$lat, 
                   lon = df$lon, 
                   dist = df$dist,
                   avg_pred = rowSums(df[, x][df[, x] < quantile(df[, x], 0.95, na.rm = T), ])/25)
avgs$pred_scaled <- rescale(avgs$avg_pred)
return(avgs)
}

# Function to generate average prediction from all predictions for second and third time periods
pred_avgs2 <- function(df){
  x <- grepl("pred", names(df), fixed = T)
  avgs <- data.frame(lat = df$lat, 
                     lon = df$lon, 
                     dist = df$dist,
                     avg_pred = rowSums(df[, x])/30)
  avgs$pred_scaled <- rescale(avgs$avg_pred)
  return(avgs)
}

# Function to calculate and plot mean predicted abundance per year
avg_plot <- function(df, range){
preds <- sapply(df, function(x) colMeans(select(x, pred), na.rm = T))
avgs <- data.frame(year = c(range), 
                   avg_pred = preds)
avgs$avg_scaled <- rescale(avgs$avg_pred)
plot <- ggplot(avgs) +
  geom_line(aes(x = year,
                y = avg_scaled))
return(plot)
}
